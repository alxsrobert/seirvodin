#' Function to generate n_samples * n_part stochastic outbreaks for a given model output
#'
#' @param model_run mcstate_pmcmc object, output from seirvodin::run_model().
#' @param model dust model.
#' @param list_specs list containing the specification of the simulations, output from seirvodin::specs_simulations().
#' @param list_data list containing all data needed for the simulations, expecting 7 elements: `ref_m`: contact matrix between age groups (matrix); `ref_d`: degree of connectivity between regions (matrix); `new_birth`: number of daily births per region (matrix); `mean_import_per_reg`: number of annual import by year and region (matrix);  `year_per_age`: duration of each age group (vector); `dt_vacc`: vaccine coverage by year and age (data table with at least 5 columns: years, region, coverage, yob (for year of birth), and dose); `N`: number of inhabitants at t0 by region and age (matrix).
#' @param n_part Number of simulation per sample.
#' @param aggreg_year boolean: whether to aggregate the number of cases simulated by year.
#' @param verbose boolean: print the simulation step for each sample.
#'
#' @return 3-dimensional array containing the number of cases by compartment (stratified by vaccine status, age and region), simulation, and day (or year if aggreg_year is TRUE)
#' @export
#'
generate_outbreaks <- function(
    model_run, model, list_specs, list_data, n_part, aggreg_year = TRUE, 
    verbose = FALSE
){
  samples <- extract_sample(
    model_run = model_run, burnin = list_specs$burnin, 
    n_samples = list_specs$n_samples
  )
  
  ## Create list containing all data
  all_data <- compute_from_data(
    list_data = list_data, list_specs = list_specs
  )
    
  if(list_specs$nowane) samples[, "v_leak"] <- 0
  
  # The rownames of all_output is a combination of vaccination status, age, and region
  names_age <- rownames(list_data$N)
  names_regions <- colnames(list_data$N)
  
  # Initialise empty 3-d array, which will contain the number of cases per 
  # vaccination status, age, region, simulation, and day
  all_output <- array(NA, dim = c(
    length(list_specs$state_names) * length(names_age) * 
      length(names_regions), nrow(samples) * n_part, 
    if(aggreg_year) list_specs$N_year else list_specs$N_time))
  
  rownames(all_output) <- 
    c(paste0(rep(paste0(
      rep(list_specs$state_names, each = length(names_regions)), 
      "_reg", rep(seq_along(names_regions), length(list_specs$state_names))), 
      each = length(names_age)), "_age", seq_along(names_age))
    )
  # for each element of sample, generate n_part stochastic outbreaks
  for(k in seq_len(nrow(samples))){
    all_output[, seq((k-1) * n_part + 1, (k) * n_part),] <- 
      generate_outbreaks_1sample(
        sample = samples[k,], model = model, data = all_data, aggreg_year = aggreg_year,
        states = rownames(all_output), n_part = n_part, waning = list_specs$waning, 
        list_specs = list_specs, year_start = list_specs$year_start, N_time = list_specs$N_time, 
        deterministic = list_specs$deterministic)
    if(verbose) print(k)
  }
  return(all_output)
}

#' Set specifications of the simulation run
#'
#' @param N_year Number of years the simulations run for.
#' @param year_start starting year of the simulations.
#' @param waning does the model include waning of immunity? Should be one of "no", "since_eli", "since_vax", or "early".
#' @param nowane boolean: if set to TRUE, v_leak will be set to 0, used to quantify the impact of waning on the overall number of cases
#' @param deterministic boolean: set to TRUE to run a deterministic version of the model
#' @param burnin Numeric: duration of the burnin period
#' @param n_samples Numeric: Number of samples used to generate the simulations
#' @param alpha Duration of incubation period.
#' @param gamma Duration of infectious period.

#' @return List of specifications
#' @export
#'
specs_simulations <- function(N_year, year_start, waning, burnin, n_samples, 
                              nowane = FALSE, deterministic = FALSE, alpha = 11, 
                              gamma = 8){
  if(!(waning %in% c("no", "since_vax", "since_eli", "early"))) stop("waning should be `no`, `since_vax`, `since_eli`, or `early`")
  ## Duration of the run (in days)
  ## Duration of the run (in days)
  N_time <- 365 * N_year
  
  ## Define time step
  dt <- 1
  
  state_names <- c("new_IS", "new_IV1", "new_IV2")
  
  return(list(state_names = state_names, year_start = year_start, N_year = N_year, 
              N_time = N_time, waning = waning, burnin = burnin, n_samples = n_samples, 
              nowane = nowane, deterministic = deterministic, alpha = alpha, gamma = gamma))
}

#' Extract samples to generate simulations
#'
#' @param model_run mcstate_pmcmc object, output from seirvodin::run_model().
#' @param burnin Numeric: duration of the burnin period
#' @param n_samples Numeric: Number of samples used to generate the simulations
#'
#' @return Matrix contaning n_samples row
#' @noRd
#'
#' @keywords internal
extract_sample <- function(model_run, burnin, n_samples){
  samples <- model_run$pars[-c(1:burnin),]
  samples <- (samples[seq(1, nrow(samples), nrow(samples)/n_samples),])
  
  return(samples)
}

#' Generate n_part simulations from a given sample
#'
#' @param sample vector containing the value of the different parameters
#' @param data list containing all datasets, output from compute_from_data().
#' @param states vector containing the states that will be returned from the model run
#' @param waning does the model include waning of immunity? "no", "since_eli", "since_vax", or "early".
#' @param year_start starting year of the simulations.
#' @param deterministic boolean: set to yes to run a deterministic version of the model
#' @param N_time Number of days the simulations run for.
#' @param aggreg_year boolean: whether to aggregate the number of cases simulated by year.
#' @param model dust model.
#' @param list_specs list containing the specification of the simulations, output from seirvodin::specs_simulations().
#' @param n_part Number of simulations run per parameter set.
#'
#' @return 3-dimensional array containing the number of cases by compartment (stratified by vaccine status, age and region), simulation, and day (or year if aggreg_year is TRUE)
#' @noRd
#'
#' @keywords internal
generate_outbreaks_1sample <- function(sample, model, data, states, list_specs,
                                       waning, year_start, deterministic, N_time,
                                       n_part = 1, aggreg_year = TRUE){
  
  S <- data$S
  R <- data$R
  N_age <- nrow(data$N)
  N_reg <- ncol(data$N)
  
  ## Define number of contacts
  beta <- sample["beta"]
  
  ## Define vaccine effectiveness against infection
  # Primary vaccine failure
  v_fail <- sample["v_fail"]
  # secondary vaccine failure (increasing with age)
  if(any(names(sample) == "v_leak")) v_leak <- sample["v_leak"]
  # Secondary vaccine failure (constant)
  if(any(names(sample) == "v_sec")) v_sec <- sample["v_sec"]
  # Vaccine coverage in individuals born in the 1970s
  if(any(names(sample) == "v_70s")) v_70s <- sample["v_70s"]
  # Against onwards infection
  vacc <- sample["vacc"]
  
  ## Define spatial kernel parameter
  b <- if(any(names(sample) == "b")) sample["b"] else 1
  c <- if(any(names(sample) == "c")) sample["c"] else 1
  theta <- if(any(names(sample) == "theta")) sample["theta"] else 1
  
  if(any(names(sample) == "alpha")) alpha <- sample["alpha"] else alpha <- list_specs[["alpha"]]
  if(any(names(sample) == "gamma")) gamma <- sample["gamma"] else gamma <- list_specs[["gamma"]]
  
  ## Define the seasonality parameters
  X <- sample["X"]
  Y <- sample["Y"]
  X_import <- sample["X_import"]
  Y_import <- sample["Y_import"]
  report_import <- sample["report_import"]
  
  # Duration of maternal antibodies
  delta <- sample["delta"]
  
  catchup <- rep(0, nrow(data$V1))
  catchup2 <- rep(0, nrow(data$V1))
  recov <- rep(0, nrow(data$V1))
  v <- rep(0, nrow(data$V1))
  
  # Set proportion of recovered in 2010 per age group
  recov <- data$R * 0
  
  for(i in seq_along(catchup)){
    if(any(names(sample) == paste0("catchup_", i))){
      # Individuals in V1 who were vaccinated during an MMR2 catchup are set as V2
      catchup[i] <- sample[[paste0("catchup_", i)]]
      data$V2[i,] <- round(data$V2[i, ] + data$V1[i, ] * (catchup[i]))
      data$V1[i,] <- round(data$V1[i, ] * (1 - catchup[i]))
    } 
    if(any(names(sample) == paste0("catchup2_", i))){
      catchup2[i] <- sample[[paste0("catchup2_", i)]]
      # Individuals who were vaccinated during a catchup campaign are set as V2
      S[i,] <- round(S[i,] * (1 - catchup2[i]))
      data$V1[i,] <- round(data$V1[i,] * (1 - catchup2[i]))
      data$V2[i,] <- round(data$V2[i,] + S[i,] * catchup2[i] + 
                             data$V1[i,] * catchup2[i])
    } 
    # Individuals who started recovered
    if(any(names(sample) == paste0("recov_", i))) recov[i,] <- sample[[paste0("recov_", i)]]
    # Individuals who started vaccinated
    if(any(names(sample) == paste0("v_", i))){
      v[i] <- sample[[paste0("v_", i)]]
      S_vax <- round(S[i,] * v[i])
      data$V1[i,] <- S_vax + data$V1[i,]
      S[i,] <- (S[i,] - S_vax)
    }
  }
  
  
  # Compute the number of importations per region 
  import_per_reg <- data$mean_import_per_reg
  
  mean_import <- import_per_reg / report_import
  
  array_cov1 <- array(data$array_cov1[-1,,], 
                      dim = c(nrow(data$array_cov1) - 1, 
                              ncol(data$array_cov1), 
                              dim(data$array_cov1)[3]),
                      dimnames = list(rownames(data$array_cov1[-1,,]), 
                                      colnames(data$array_cov1))
  )
  array_cov2 <- array(data$array_cov2[-1,,], 
                      dim = c(nrow(data$array_cov2) - 1, 
                              ncol(data$array_cov2), 
                              dim(data$array_cov2)[3]),
                      dimnames = list(rownames(data$array_cov2[-1,,]), 
                                      colnames(data$array_cov2))
  )
  
  ## Create the age and region-stratified model
  seir_model <- model$new(pars = list(
    m = data$ref_m, d = data$ref_d, b = b, c = c, theta = theta,
    beta = beta, X = X, Y = Y, delta = delta, X_import = X_import, 
    Y_import = Y_import, v_fail = v_fail, vacc = vacc,
    v_leak = if(exists("v_leak")) v_leak else 0, 
    v_sec = if(exists("v_sec")) v_sec else 0, 
    v_70s = if(exists("v_70s")) v_70s else 0, 
    mean_import = mean_import, 
    N_time = N_time, N_age = N_age, N_reg = N_reg, V1_init = data$V1, 
    V2_init = data$V2, 
    Es_init = S * 0, Ev1_init = S * 0, Ev2_init = S * 0,
    Is_init = S * 0, Iv1_init = S * 0, Iv2_init = S * 0, 
    S_init = S + R, recov = recov, 
    R_init = R * 0, RV1_init = R * 0, RV2_init = R * 0,
    
    array_cov1 = array_cov1, array_cov2 = array_cov2,
    array_new = data$new_birth, 
    dt = 1, 
    alpha = 1 / alpha, gamma = 1 / gamma, 
    waning = if(waning == "no") 0 else if(waning == "since_vax") 1 else
      if (waning == "since_eli") 2 else if (waning == "early") 3,
    year_start = year_start, year_per_age = data$year_per_age
  ),
  time = 1, n_particles = n_part, n_threads = 1L, seed = 1L, deterministic = deterministic)
  
  ## Define array output_sim, containing the number of individuals in each compartment per day
  ## in each simulation
  output_sim <- array(NA, dim = c(seir_model$info()$len, n_part, N_time))
  
  ## Rename the rows of output_sim, as the name of the state + region + age group
  all_states <- c("S", "M", "V1", "V2", "V1p", "V2p", "R", "RV1", "RV2", 
                  "Es", "Ev1", "Ev2", "Is", "Iv1", "Iv2", "new_IS", "new_IV1", "new_IV2")
  ## The first row of output_sim contains the time
  rownames(output_sim) <- c(
    "Time", "iter", "new_IV1_tot", "new_IV2_tot", 
    paste0(rep(paste0(rep(all_states, each = nrow(data$ref_d)), "_reg", 
                      rep(seq_len(nrow(data$ref_d)), length(all_states))), 
               each = nrow(data$ref_m)), "_age", seq_len(nrow(data$ref_m)))
  )
  
  ## Run the model for each time step
  if(aggreg_year){
    output_year <- array(NA, dim = c(length(states), ncol(output_sim), dim(output_sim)[3]/365),
                         dimnames = list(states))
  }
  iter <- 0
  for (t in seq_len(N_time)) {
    ## Generate stochastic simulations
    output_sim[ , , t] <- seir_model$run(t)
    if(t%%365.25 < 1 & aggreg_year){
      # aggregate by year
      iter <- t%/%365.25
      output_year[,,iter] <- 
        apply(output_sim[states, , seq_len(N_time) %/% 365.25 == (iter - 1)], 
              c(1,2), sum)
    } else if(t == N_time & aggreg_year){
      # aggregate last year of simulation
      iter <- iter + 1
      output_year[,,iter] <- 
        apply(output_sim[states, , seq_len(N_time) %/% 365.25 == (iter - 1)], 
              c(1,2), sum)
    }
  }
  if(aggreg_year) return(output_year) else return(output_sim[states, , ])
}

