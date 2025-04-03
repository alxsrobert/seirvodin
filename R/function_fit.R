#' Run odin.dust model
#'
#' @param list_specs list containing the specification of the model fitting process, output from seirvodin::specs_run().
#' @param list_data list containing all data needed for model fitting, expecting 8 elements: `ref_m`: contact matrix between age groups (matrix); `ref_d`: degree of connectivity between regions (matrix); `new_birth`: number of daily births per region (matrix); `mean_import_per_reg`: number of annual import by year and region (matrix);  `year_per_age`: duration of each age group (vector); `dt_vacc`: vaccine coverage by year and age (data table with at least 5 columns: years, region, coverage, yob (for year of birth), and dose); `N`: number of inhabitants at t0 by region and age (matrix);`dt_case`: data table containing the number of cases by date, age and region (three columns: date (numeric starting at 0); cases; and population (character containing unique id for a strata of vaccine status, age, and region))
#' @param model dust model.
#' @param proposal_matrix Proposal matrix to run the MCMC, leave to NULL to start with default matrix
#' @param init Vector of initial value for each parameter
#' @param prior List of prior function for each parameter (for parameters without prior, the element should be  NULL)
#' @param list_min_max list containing two vectors: `min` and `max`, each list contains the minimum (and maximum) value for each parameter
#'
#' @return An `mcstate_pmcmc` object, containing all information on the model fit.
#' @export
run_model <- function(list_specs, list_data, init, list_prior, list_min_max, model, proposal_matrix = NULL){
  
  ## Create list containing all data
  all_data <- compute_from_data(
    list_data = list_data, list_specs = list_specs
  )
  
  ## Initialise list of parameters
  mcmc_pars <- create_mcmc_pars(
    list_specs = list_specs, list_data = all_data, init =init, list_prior = list_prior, 
    list_min_max = list_min_max, proposal_matrix = proposal_matrix
  )

  data_wide <- pivot_wider(list_data$dt_case, names_from = "population", values_from = "cases")
  
  ## Import data
  particle_data <- mcstate::particle_filter_data(
    data_wide, time = "date", rate = 1, initial_time = 0
  )
  
  ## Initialise the model using the parameters previously defined
  filter <- mcstate::particle_deterministic$new(
    particle_data, model, case_compare, index
  )
  
  ## Create control object
  control <- mcstate::pmcmc_control(
    n_steps = list_specs$n_steps, save_state = FALSE,  adaptive_proposal = TRUE, 
    save_trajectories = FALSE, progress = TRUE)
  
  ## Run model
  pmcmc_run <- mcstate::pmcmc(
    mcmc_pars, filter, control = control)
  
  return(pmcmc_run)
}

#'  Set specifications of the model run
#'
#' @param N_year Number of years the simulations run for.
#' @param year_start starting year of the simulations.
#' @param n_steps Number of MCMC steps.
#' @param waning Character values, corresponds to whether waning is included in the model, expect one of three values: "no", "since_vax", or "since_eli"
#' @param alpha Duration of incubation period.
#' @param gamma Duration of infectious period.
#'
#' @return List of specifications
#' @export
#'
specs_run <- function(N_year, year_start, n_steps, waning, alpha = 11, gamma = 8){
  if(!(waning %in% c("no", "since_vax", "since_eli"))) stop("waning should be `no`, `since_vax`, or `since_eli`")
  ## Duration of the run (in days)
  N_time <- 365 * N_year
  
  ## Define time step
  dt <- 1
  
  state_names <- c("new_IS", "new_IV1", "new_IV2")
  
  return(list(state_names = state_names, year_start = year_start, N_year = N_year, 
              N_time = N_time, n_steps = n_steps, waning = waning, alpha = alpha, 
              gamma = gamma))
}

#' Set filter function  
#'
#' @noRd
#' @keywords internal
case_compare <- function(state, observed, pars = NULL) {
  incidence_modelled <- state[,, drop = TRUE]
  ## Remove the first four entries of observed (date_start, date_end, time_start, time_end)
  incidence_observed <- unlist(observed)[-seq_len(4)]
  lambda <- incidence_modelled + 1e-10
  if(length(incidence_modelled) > length(incidence_observed)){
    colSums(dpois(x = incidence_observed, lambda = lambda, log = TRUE))
  } else sum(dpois(x = incidence_observed, lambda = lambda, log = TRUE))
}

#' Index function to return only states IS, IV1, and IV2
#'
#'
#'
#' @noRd
#' @keywords internal
index <- function(info) {
  list(run = c(new_IS = info$index$new_IS
               , new_IV1 = info$index$new_IV1,
               new_IV2 = info$index$new_IV2
  ), state = c())
}

