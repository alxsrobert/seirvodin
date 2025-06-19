#' Function to create mcstate::pmcmc_parameters
#'
#' @param list_data list containing all data, output from `compute_from_data()`.
#' @param list_specs list containing the specification, expecting at least 2 elements: `year_start`: starting year of the simulations; `N_year`: Number of years the simulations run for.
#' @param proposal_matrix Proposal matrix to run the MCMC, if NULL it will be set to default version.
#' @param init Vector of initial value for each parameter
#' @param list_prior List of prior function for each parameter (for parameters without prior, the element should be  NULL)
#' @param list_min_max list containing two vectors: `min` and `max`, each list contains the minimum (and maximum) value for each parameter
#'
#' @return A `pmcmc_parameters` object, containing the initial parameters, the proposal matrix and the transform function
#'
#' @noRd
#' @keywords internal
create_mcmc_pars <- function(list_data, list_specs, init, list_prior, list_min_max, proposal_matrix){
  ## Extract elements from list_data (generated with the function compute_from_data)
  ref_d <- list_data[["ref_d"]]
  ref_m <- list_data[["ref_m"]]
  new_birth <- list_data[["new_birth"]]
  mean_import_per_reg <- list_data[["mean_import_per_reg"]]
  N_age <- nrow(list_data[["N"]])
  N_reg <- ncol(list_data[["N"]])
  N_time <- list_specs[["N_time"]]
  S <- list_data[["S"]]
  V1 <- list_data[["V1"]]
  V2 <- list_data[["V2"]]
  RV1 <- RV2 <- R <- list_data[["R"]]
  array_cov1 <- list_data[["array_cov1"]]
  array_cov2 <- list_data[["array_cov2"]]
  year_per_age <- list_data[["year_per_age"]]
  
  if(!all(names(init) %in% names(list_min_max[["min"]])) ||
     !all(names(init) %in% names(list_min_max[["max"]]))){
    stop("init, list_min_max[[`min`]], and list_min_max[[`max`]] should all contain the same elements")
  }
  
  # Create mcstate::pmcmc_parameter from initial, prior, and min_max
  list_pars <- create_pars(init = init, prior = list_prior, min = list_min_max[["min"]], 
                           max = list_min_max[["max"]])
  
  # Create proposal matrix
  if(is.null(proposal_matrix)){
    proposal_matrix <- matrix(0, nrow = length(list_pars), ncol = length(list_pars))
    colnames(proposal_matrix) <- rownames(proposal_matrix) <- names(list_pars)
    diag(proposal_matrix) <- .01
  }
  
  # Function for parameter transformation
  transform <- make_transform(
    list_pars, m = ref_m, d = ref_d, import = mean_import_per_reg, 
    N_time = N_time, N_age = N_age, N_reg = N_reg, 
    
    S_init = S, V1_init = V1, V2_init = V2, R_init = R, RV1_init = RV1, RV2_init = RV2, 
    # All E and I compartments start empty
    Es_init = S * 0, Ev1_init = S * 0, Ev2_init = S * 0, 
    Is_init = S * 0, Iv1_init = S * 0, Iv2_init = S * 0,
    
    array_cov1 = array(array_cov1[-1,,], 
                       dim = c(nrow(array_cov1) - 1, ncol(array_cov1),
                               dim(array_cov1)[3]),
                       dimnames = list(rownames(array_cov1[-1,,]), 
                                       colnames(array_cov1))), 
    array_cov2 = array(array_cov2[-1,,], 
                       dim = c(nrow(array_cov2) - 1, ncol(array_cov2), 
                               dim(array_cov2)[3]),
                       dimnames = list(rownames(array_cov2[-1,,]), 
                                       colnames(array_cov2))), 
    array_new = new_birth, 
    
    waning = if(list_specs$waning == "no") 0 else if(list_specs$waning == "since_vax") 1 else
      if (list_specs$waning == "since_eli") 2, 
    
    alpha = list_specs[["alpha"]], gamma = list_specs[["gamma"]], 
    
    dt = 1, year_start = list_specs$year_start, year_per_age = year_per_age,
    # If the following parameters are included in pars, this reference level will be ignored
    b = 1, c = 1, theta = 1, v_leak = 0, v_sec = 0
  )
  
  ## create pmcmc_parameters object
  mcmc_pars <- mcstate::pmcmc_parameters$new(
    parameters = list_pars, proposal = proposal_matrix, transform = transform
  )
  
  return(mcmc_pars)
}

#' Transformation function to move between inferred parameters and parameters 
#' needed to implement the model (see https://mrc-ide.github.io/mcstate/reference/pmcmc_parameters.html)
#'
#' @param pars list with each element corresponding to a fitted parameter of the model. Each element must be a `pmcmc_parameter` object (generated with mcstate::pmcmc_parameter).
#' @param m Contact matrix between age groups
#' @param d Connectivity matrix between regions
#' @param import Average number of importations per year and region
#' @param N_time Number of days the model runs for
#' @param N_age Number of age groups
#' @param N_reg Number of regions
#' @param b Spatial parameter 1
#' @param c Spatial parameter 2
#' @param theta Spatial parameter 3
#' @param v_leak Waning of immunity
#' @param v_sec Rate of secondary vaccine failure (without waning)
#' @param S_init Number of individuals in Susceptible compartment (by age and region)
#' @param V1_init Number of individuals vaccinated once (by age and region)
#' @param V2_init Number of individuals vaccinated twice (by age and region)
#' @param Es_init Number of unvaccinated individuals in the Exposed compartment (by age and region)
#' @param Ev1_init Number of individuals vaccinated once in the Exposed compartment (by age and region)
#' @param Ev2_init Number of individuals vaccinated twice in the Exposed compartment (by age and region)
#' @param Is_init Number of unvaccinated individuals in the Infected compartment (by age and region)
#' @param Iv1_init Number of individuals vaccinated once in the Infected compartment (by age and region)
#' @param Iv2_init Number of individuals vaccinated twice in the Infected compartment (by age and region)
#' @param R_init Number of unvaccinated individuals in the Recovered compartment (by age and region)
#' @param RV1_init Number of individuals vaccinated once in the Recovered compartment (by age and region)
#' @param RV2_init Number of individuals vaccinated twice in the Recovered compartment (by age and region)
#' @param array_cov1 Daily rate of vaccination (from S to V1)
#' @param array_cov2 Daily rate of vaccination (from V1 to V2)
#' @param array_new Daily number of births
#' @param alpha Duration of incubation period.
#' @param gamma Duration of infectious period.
#' @param dt time step
#' @param year_per_age Duration (in year) spent in each age group.
#' @param waning Character values, corresponds to whether waning is included in the model, expect one of three values: "no", "since_vax", or "since_eli"
#' @param year_start starting year of the simulations.
#'
#' @noRd
#' @keywords internal
make_transform <- function(pars,m, d, import, N_time, N_age, N_reg, 
                           b, c, theta, v_leak, v_sec, 
                           S_init, V1_init, V2_init, Es_init, 
                           Ev1_init, Ev2_init, Is_init, Iv1_init, Iv2_init, 
                           R_init, RV1_init, RV2_init, array_cov1, array_cov2, 
                           array_new, dt, year_per_age, waning, alpha, gamma,
                           year_start){
  function(pars){
    # Extract parameters from the vector pars
    beta <- pars[["beta"]]
    delta <- pars[["delta"]]
    X <- pars[["X"]]
    Y <- pars[["Y"]]
    X_import <- pars[["X_import"]]
    Y_import <- pars[["Y_import"]]
    v_fail <- pars[["v_fail"]]
    vacc <- pars[["vacc"]]
    report_import <- pars[["report_import"]]

    # Extract optional parameters
    if(any(names(pars) == "v_leak")) v_leak <- pars[["v_leak"]]
    if(any(names(pars) == "v_sec")) v_sec <- pars[["v_sec"]]
    if(any(names(pars) == "b")) b <- pars[["b"]]
    if(any(names(pars) == "c")) c <- pars[["c"]]
    if(any(names(pars) == "theta")) theta <- pars[["theta"]]
    if(any(names(pars) == "alpha")) alpha <- pars[["alpha"]]
    if(any(names(pars) == "gamma")) gamma <- pars[["gamma"]]

    catchup <- rep(0, nrow(V1_init))
    catchup2 <- rep(0, nrow(V1_init))
    recov <- rep(0, nrow(V1_init))
    v <- rep(0, nrow(V1_init))
    
    V_tot <- V1_init + V2_init
    # Set proportion of recovered at year_start per age group
    recov <- R_init * 0
    
    for(i in seq_along(catchup)){
      if(any(names(pars) == paste0("catchup_", i))){
        # Individuals in V1 who were vaccinated during an MMR2 catchup in 1996 are set as V2
        catchup[i] <- pars[[paste0("catchup_", i)]]
        V1_init[i,] <- round(V_tot[i, ] * (1 - catchup[i]))
        V2_init[i,] <- round(V_tot[i, ] * (catchup[i]))
      } 
      if(any(names(pars) == paste0("catchup2_", i))){
        catchup2[i] <- pars[[paste0("catchup2_", i)]]
        # Individuals who were vaccinated during a catchup campaign are set as V2
        S_init[i,] <- round(S_init[i,] * (1 - catchup2[i]))
        V1_init[i,] <- round(V1_init[i,] * (1 - catchup2[i]))
        V2_init[i,] <- round(V2_init[i,] + S_init[i,] * catchup2[i] + 
                               V1_init[i,] * catchup2[i])
      } 
      # Individuals who started recovered
      if(any(names(pars) == paste0("recov_", i))) recov[i,] <- pars[[paste0("recov_", i)]]
      # Individuals who started vaccinated
      if(any(names(pars) == paste0("v_", i))){
        V1_init[i,] <- round(S_init[i,] * v[i])
        S_init[i,] <- (S_init[i,] - V1_init[i,])
      }
    }
    list(beta = beta, X = X, Y = Y, delta = delta, X_import = X_import, Y_import = Y_import,
         b = b, c = c, theta = theta, m = m, d = d, vacc = vacc, 
         mean_import = import / report_import, N_time = N_time, 
         N_age = N_age, N_reg = N_reg, v_fail = v_fail, v_leak = v_leak, v_sec = v_sec, 
         V1_init = V1_init, V2_init = V2_init, Es_init = Es_init, 
         Ev1_init = Ev1_init, Ev2_init = Ev2_init, Is_init = Is_init, 
         Iv1_init = Iv1_init, Iv2_init = Iv2_init, 
         
         alpha = 1 / alpha, gamma = 1 / gamma, 
         
         recov = recov, R_init = R_init * 0, RV1_init = R_init * 0, RV2_init = R_init * 0,
         S_init = S_init + R_init,
         
         array_cov1 = array_cov1, array_cov2 = array_cov2, 
         array_new = array_new, dt = dt, year_per_age = year_per_age,
         waning = waning, year_start = year_start)
  }
}

#' Create list of pmcmc_parameter object
#'
#' @param init Vector of initial value for each parameter
#' @param prior List of prior function for each parameter
#' @param min List of minimum value for each parameter
#' @param max List of maximum value for each parameter
#'
#' @return List of pmcmc_parameter
#'
#' @noRd
#' @keywords internal
create_pars <- function(init, prior, min, max){
  pars <- list()
  for(i in seq_along(init)){
    name_i <- names(init)[i]
    if(!any(names(prior) == name_i)) p <- NULL else p <- prior[[name_i]]
    pars <- append(pars, list(mcstate::pmcmc_parameter(
      name = names(init)[i], initial = init[[name_i]], min = min[[name_i]], 
      max = max[[name_i]], prior = p)))
  }
  return(pars)
}
