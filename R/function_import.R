#' Compute initial distribution, and daily vaccination rate 
#'
#' @param list_data list containing all data needed for model fitting, expecting at least 3 elements: `year_per_age`: duration of each age group (vector); `dt_vacc`: vaccine coverage by year and age (data table with at least 5 columns: years, region, coverage, yob (for year of birth), and dose); `N`: number of inhabitants at t0 by region and age (matrix)
#' @param list_specs list containing the specification, expecting at least 2 elements: `year_start`: starting year of the simulations; `N_year`: Number of years the simulations run for.
#'
#' @return List of data sets: all elements from list_data plus initial distribution and daily vaccination rate.
#'
#' @noRd
#' @keywords internal
compute_from_data <- function(list_data, list_specs){
  ## Compute initial state
  list_state <- compute_initial_state(
    dt_vacc = list_data$dt_vacc, 
    N = list_data$N, 
    year_per_age = list_data$year_per_age,
    year_start = list_specs$year_start)
  
  ## Compute vaccine coverage through time
  list_array_cov <- compute_vax_cov(
    dt_vacc = list_data$dt_vacc, 
    N = list_data$N,
    year_start = list_specs$year_start, 
    N_year = list_specs$N_year)
  
  list_data <- append(list_data, append(list_state, list_array_cov))
  
  return(list_data)
}

#' Compute the distribution of individuals in each compartment at t = 0
#'
#' @param dt_vacc data table containing vaccine coverage by year and age
#' @param N Matrix containing the number of inhabitants per age group and region at t=0
#' @param year_per_age Duration (in year) spent in each age group.
#' @param year_start starting year.
#'
#' @noRd
#' @keywords internal
compute_initial_state <- function(dt_vacc, N, year_per_age, year_start){
  #### Define initial number of cases in compartment
  regions <- colnames(N)
  age_names <- rownames(N)
  N_reg <- length(regions)
  N_age <- length(age)
  
  # avoid data.table NOTE
  yob <- dose <- years <- region <- NULL

  ## Define vaccine uptake and proportion of recovered at t = 0
  ## Initialise empty matrices for proportion of recovery and coverage
  recov <- matrix(0, ncol = N_reg, nrow = N_age, byrow = F)
  cov1 <- matrix(0, ncol = N_reg, nrow = N_age)
  cov2 <- matrix(0, ncol = N_reg, nrow = N_age)
  
  colnames(cov1) <- colnames(cov2) <- toupper(regions)
  
  mean_vacc1 <- numeric(N_reg)
  mean_vacc2 <- numeric(N_reg)
  nb_years <- 0
  yob_transition <- c(year_start - cumsum(year_per_age))
  
  for(j in seq_along(yob_transition[-1])){
    age <- j + 1
    max_year <- yob_transition[j]
    min_year <- ifelse(j == length(yob_transition), -Inf, yob_transition[j+1] + 1)
    
    if(max_year >= min(dt_vacc$yob)){
      if(length(unique(dt_vacc[yob %in% seq(min_year, max_year) & dose == 1, years])) == 1){
        last_vacc <- unique(dt_vacc[yob %in% seq(min_year, max_year) & dose == 1, years])
      } else{
        last_vacc <- max(dt_vacc[yob %in% seq(min_year, max_year) & dose == 1 & years <= year_start, years])
      }

      cov1[age, unique(dt_vacc[yob <= max_year & yob >= min_year & dose == 1, region])] <- 
        dt_vacc[yob <= max_year & yob >= min_year & dose == 1 & years == last_vacc, 
                lapply(.SD, mean), by = region, .SDcols = "coverage"]$coverage
      
      if(any(dt_vacc$yob <= max_year & dt_vacc$yob >= min_year & dt_vacc$dose == 2 & 
             dt_vacc$years == last_vacc)){
        cov2[age, unique(dt_vacc[yob <= max_year & yob >= min_year & dose == 2, region])] <- 
          dt_vacc[yob <= max_year & yob >= min_year & dose == 2 & years == last_vacc, 
                  lapply(.SD, mean), by = region, .SDcols = "coverage"]$coverage
      }
    }
  }
  
  
  cov1 <- cov1 - cov2
  cov1[cov1 < 0] <- 0
  
  unvax <- 1 - cov1 - cov2
  
  if(any(unvax > 1 | unvax < 0)) stop("Proportion should be between 0 and 1")
  
  ## Number of inhabitants in each compartment
  V1 <- round(N * (1 - recov) * cov1, 0)
  V2 <- round(N * (1 - recov) * cov2, 0)
  S <- round(N * (1 - recov) * unvax, 0)
  R <- N * 0
  
  return(list(S = S, V1 = V1, V2 = V2, R = R))
}

#' Compute vaccine coverage at each year
#'
#' @param dt_vacc data table contaning vaccine coverage by year and age
#' @param year_start starting year.
#' @param N Matrix containing the number of inhabitants per age group and region at t=0
#' @param N_year Number of years the simulations run for.
#'
#'
#' @noRd
#' @keywords internal
compute_vax_cov <- function(dt_vacc, N, year_start, N_year){
  # Compute the number of regions and age
  regions <- colnames(N)
  age_names <- rownames(N)
  N_reg <- length(regions)
  N_age <- length(age_names)
  
  # avoid data.table NOTE
  agegroup <- age <- v1 <- v2 <- coverage <- year <- NULL
  
  # Create empty data table to compute the vaccine coverage per region / age / year
  vacc_per_age <- data.table(regions = rep(toupper(regions), N_age * N_year), 
                             age = rep(rep(age_names, each = N_reg), N_year),
                             year = rep(year_start + 0:(N_year - 1), each = N_age * N_reg))
  
  vacc_per_age[, agegroup := substr(age, 4, 5)]
  vacc_per_age[age == "age0", v1 := 0]
  vacc_per_age[age == "age0", v2 := 0]
  vacc_per_age[age %in% c("age1", "age2", "age3", "age4", "age5"),
               v1 := dt_vacc[paste(year, tolower(regions), 1, 
                                   agegroup, sep = "_"), coverage]]
  vacc_per_age[age %in% c("age1", "age2", "age3", "age4", "age5"),
               v2 := dt_vacc[paste(year, tolower(regions), 2, 
                                   agegroup, sep = "_"), coverage]]
  vacc_per_age[, agegroup := NULL]
  
  vacc_per_age[is.na(v1), v1 := 0]
  vacc_per_age[is.na(v2), v2 := 0]
  
  # Proportion of the population vaccinated once
  vacc_per_age$v1 <- vacc_per_age$v1 - vacc_per_age$v2
  vacc_per_age$v1[vacc_per_age$v1 < 0] <- 0
  
  vacc_per_age$regions <- factor(vacc_per_age$regions, levels = toupper(regions))
  
  vacc_per_age$age <- factor(vacc_per_age$age, levels = unique(vacc_per_age$age))
  
  # re-order rows
  vacc_per_age <- vacc_per_age[order(vacc_per_age$year, vacc_per_age$regions, 
                                     vacc_per_age$age), ]
  
  ## Create 3-d arrays containing the vaccine coverage (v1 then v2) per age / region / DAY
  array_cov1 <- array(data = vacc_per_age$v1, 
                      dim = c(N_age, N_reg, N_year),
                      dimnames = list(unique(vacc_per_age$age), 
                                      unique(vacc_per_age$regions))
  )[,, rep(seq_len(N_year), each  = 365)]
  
  array_cov2 <- array(data = vacc_per_age$v2, 
                      dim = c(N_age, N_reg, N_year),
                      dimnames = list(unique(vacc_per_age$age), 
                                      unique(vacc_per_age$regions))
  )[,, rep(seq_len(N_year), each  = 365)]
  return(list(array_cov1 = array_cov1, array_cov2 = array_cov2))
}

