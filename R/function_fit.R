#' Run odin.dust model
#'
#' @param list_specs list containing the specification of the model fitting 
#' process, output from seirvodin::specs_run().
#' @param list_data list containing all data needed for model fitting, 
#' expecting 8 elements: `ref_m`: contact matrix between age groups (matrix);
#' `ref_d`: degree of connectivity between regions (matrix);
#' `new_birth`: number of daily births per region (matrix); 
#' `mean_import_per_reg`: number of annual import by year and region (matrix);  
#' `year_per_age`: duration of each age group (vector); 
#' `dt_vacc`: vaccine coverage by year and age (data table with at least 5 
#' columns: years, region, coverage, yob (for year of birth), and dose); 
#' `N`: number of inhabitants at t0 by region and age (matrix);
#' `dt_case`: data table containing the number of cases by date, age and region
#' (three columns: date (numeric starting at 0); cases; and population
#' (character containing unique id for each strata of vaccine status, age, and
#' region))
#' @param model dust model.
#' @param proposal_matrix Proposal matrix to run the MCMC, leave to NULL to 
#' start with default matrix
#' @param init Vector of initial value for each parameter
#' @param list_prior List of prior function for each parameter (for parameters 
#' without prior, the element should be  NULL)
#' @param list_min_max list containing two vectors: `min` and `max`, each list
#' contains the minimum (and maximum) value for each parameter
#'
#' @importFrom tidyr pivot_wider
#' @return An `mcstate_pmcmc` object, containing all information on the model fit.
#' @export
#' @examples 
#' N_year = 5
#' # Define the specifications
#' specs <- specs_run(N_year = N_year, year_start = 2020, n_steps = 20, 
#'                    waning = "no", alpha = 11, gamma = 8)
#' 
#' # Contact matrix between age groups
#' polymod <- socialmixr::polymod
#' ref_m <- 1e6 * socialmixr::contact_matrix(
#'   survey = polymod,
#'   countries = "United Kingdom",
#'   age.limits = c(0, 2, 5, 10),
#'   symmetric = TRUE, per.capita = TRUE)$matrix.per.capita
#' 
#' age <- c("[0-2)", "[2-5)", "[5-10)", "10+")
#' 
#' # Duration of each age group
#' year_per_age <- c(2, 3, 5, 70)
#' 
#' # Degree of connectivity between regions
#' ref_d <- matrix(c(1,2,2,1), nrow = 2, ncol = 2)
#' 
#' # Number of daily births per region
#' new_birth <- rbind(A = rexp(365 * N_year, 1/(500/365)), 
#'                    B = rexp(365 * N_year, 1/(1000/365)))
#' 
#' # Number of inhabitants at t0 by region and age
#' N <- cbind(A = c(1000, 1500, 3000, 50000), B = c(2000, 3000, 75000, 90000))
#' rownames(N) <- age
#' 
#' # Number of annual import by year and region
#' mean_import_per_reg <- cbind(A = round(rexp(N_year, 1/5)),
#'                              B = round(rexp(N_year, 1/10)))
#' 
#' 
#' # Vaccine coverage by year and age
#' # first: vaccine coverage for the first dose
#' dt_vacc_dose1 <- data.frame(
#'   years = 2010:2024, region = rep(c("A", "B"), each = N_year + 10), 
#'   coverage = runif(30, .85, .97), yob = rep(2010:2024 - 2, 2), dose = 1)
#' # second: vaccine coverage for the second dose
#' dt_vacc_dose2 <- data.frame(
#'   years = 2010:2024, region = rep(c("A", "B"), each = N_year + 10), 
#'   coverage = runif(30, .75, .85), yob = rep(2010:2024 - 5, 2), dose = 2)
#' # Merge both vaccine coverage data
#' dt_vacc <- data.table::as.data.table(rbind.data.frame(dt_vacc_dose1, dt_vacc_dose2))
#' 
#' 
#' # Import daily number of cases by region, age group and vaccination status
#' data("simulated_outbreak")
#' 
#' ## Create unique row id
#' simulated_outbreak[, population := paste(vaccinated, region, age_groups, 
#'                                          sep = "_")]
#' simulated_outbreak[, population := factor(population)]
#' ## Select columns dates, cases, and population
#' simulated_outbreak <- simulated_outbreak[, .(date, cases, population)]
#' 
#' 
#' # Prepare the data
#' list_data <- list(
#'   ref_m = ref_m,
#'   year_per_age = year_per_age,
#'   ref_d = ref_d,
#'   new_birth = new_birth,
#'   N = N,
#'   mean_import_per_reg = mean_import_per_reg,
#'   dt_vacc = dt_vacc,
#'   dt_case = simulated_outbreak
#' )
#' 
#' 
#' # Define initial parameter values
#' init <- c(beta = 0.5, delta = 0.1, X = 0.1, Y = 0.1, X_import = 0.1, 
#'           Y_import = 0.1, v_fail = 0.02, vacc = 0.5, report_import = 0.5, 
#'           recov_3 = .5, recov_4 = .9, catchup2_3 = .2, b = .5, c = .5, 
#'           theta = .5)
#' # Define minimum and maximum values for each parameter
#' list_min_max <- list(
#'   min = c(beta = 0, delta = 0.003, X = 0, Y = 0, X_import = 0, Y_import = 0, 
#'           v_fail = 0, vacc = 0, report_import = 0, recov_4 = 0, recov_3 = 0, 
#'           catchup2_3 = 0, b = 0, c = 0, theta = 0),
#'   max = c(beta = 50, delta = 0.2, X = 1, Y = 7, X_import = 1, Y_import = 7, 
#'           v_fail = 1, vacc = 1, report_import = 1, recov_4 = 1, recov_3 = 1, 
#'           catchup2_3 = 1, b = 1, c = 1, theta = 1)
#' )
#' # Define prior functions
#' list_prior <- list(
#'   delta = function(x) dnorm(1/x, mean = 100, sd = 15, log = TRUE),
#'   v_fail = function(x) dnorm(x, mean = 0.02, sd = 0.01, log = TRUE),
#'   report_import = function(x) dnorm(x, mean = 0.5, sd = 0.1, log = TRUE)
#' )
#' # Run the model
#' model_run <- run_model(
#'   specs, list_data, init, list_prior, list_min_max, seirv_age_region)
#' 
run_model <- function(list_specs, list_data, init, list_prior, list_min_max, 
                      model, proposal_matrix = NULL){
  
  ## Create list containing all data
  all_data <- compute_from_data(
    list_data = list_data, list_specs = list_specs
  )
  
  ## Initialise list of parameters
  mcmc_pars <- create_mcmc_pars(
    list_specs = list_specs, list_data = all_data, init =init, 
    list_prior = list_prior, list_min_max = list_min_max, 
    proposal_matrix = proposal_matrix
  )

  data_wide <- pivot_wider(list_data$dt_case, names_from = "population", 
                           values_from = "cases")
  
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
#' @param waning Character values, corresponds to whether waning is included in
#'  the model, expect one of three values: "no", "since_vax", or "since_eli"
#' @param alpha Duration of incubation period.
#' @param gamma Duration of infectious period.
#'
#' @return List of specifications
#' @export
#'
specs_run <- function(
    N_year, year_start, n_steps, waning, alpha = 11, gamma = 8){
  if(!(waning %in% c("no", "since_vax", "since_eli"))) 
    stop("waning should be `no`, `since_vax`, or `since_eli`")
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
#' @importFrom stats dpois
case_compare <- function(state, observed, pars = NULL) {
  incidence_modelled <- state[,, drop = TRUE]
  ## Remove the first four entries of observed (date_start, date_end, 
  ## time_start, time_end)
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

