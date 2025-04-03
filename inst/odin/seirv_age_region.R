#### Include transmission between regions (using a spatial kernel) and age groups
#### Include outside importations

#### Structure of the script:
#### 1- Core equations
#### 2- Draw number of individuals changing compartments
#### 3- Compute the probability of movement between compartments
#### 4- compute the force of infection
#### 5- Compute the number of importations
#### 6- Compute ageing
#### 7- Define initial status / dimension of the compartment matrices

#### 1- Core equations for transitions between compartments: ####

### All compartments are stratified by age (ROWS) and regions (COLUMNS)
## Susceptibles compartments: 
# New value of S, V1, and V2 is: 
#   (number of individuals at the previous time step) +
#   (nb of individuals who aged into the compartment) -
#   (nb of individuals from the compartment who aged) - 
#   (nb of susceptible individuals who were exposed)

# First age groups: Add births
# M is the number of children protected by maternal immunity
# All values at age groups higher than 1 is 0
# n_MS: number of people moving from M to S (losing maternal immunity at t)
update(M[1, ]) <- M[1, j] +
  new_birth[j] -             # Births
  n_MS[j]                    # Losing maternal immunity
# S number of people in the Susceptible compartment
# S_Ses: number of people from S to "Exposed and unvaccinated"
update(S[1, ]) <- S[1, j] +  # Number of susceptible age 0 to 1 at the previous time step
  n_MS[j] - n_S[1, j] -      # from maternal immun + ageing out
  n_SEs[1, j]                # Exposed
# In the first age group, V1. V2. V1p, and V2p should be 0
# V1 number of people in the single-vaccinated compartment (with vaccine failure)
update(V1[1,]) <- V1[1, j] - # Initial
  n_V1[1, j] -               # Ageing out
  n_v1E[1, j]                # Exposed
# V2 number of people in the double-vaccinated compartment (double vaccine failure)
update(V2[1,]) <- V2[1, j] - # Initial
  n_V2V2[1, j] -             # Ageing out
  n_v2E[1, j]                # Exposed

# V1 number of people in the single-vaccinated and protected compartment
update(V1p[1,]) <- V1p[1, j] -# Initial
  n_V1p[1, j] -               # Ageing out
  n_v1pE[1, j]                # Exposed
# V1 number of people in the double-vaccinated and protected compartment
update(V2p[1,]) <- V2p[1, j] -# Initial
  n_V2pV2p[1, j] -            # Ageing out
  n_v2pE[1, j]                # Exposed


# Following age groups
update(S[2 : len_ageing,]) <- S[i, j] +   # Initial
  n_SS[i - 1, j] -                        # Ageing in
  n_S[i, j] -                             # Ageing out
  n_SEs[i, j]                             # Exposed

update(V1[2 : len_ageing,]) <- V1[i, j] + # Initial
  n_V1V1[i - 1, j] + n_SV1[i - 1, j] -    # Ageing in (from vaccinated and susceptible)
  n_V1[i, j] -                            # Ageing out
  n_v1E[i, j]                             # Exposed

update(V2[2 : len_ageing,]) <- V2[i, j] + # Initial
  n_V2V2[i - 1, j] + n_V1V2[i - 1, j] -   # Ageing in (from already double vaccinated, or single vaccinated)
  n_V2V2[i, j] -                          # Ageing out
  n_v2E[i, j]                             # Exposed

update(V1p[2 : len_ageing,]) <- V1p[i, j] + # Initial
  n_V1pV1p[i - 1, j] + n_SV1p[i - 1, j] -   # Ageing in
  n_V1p[i, j] -                             # Ageing out
  n_v1pE[i, j]                              # Exposed

update(V2p[2 : len_ageing,]) <- V2p[i, j] +                      # Initial
  n_V2pV2p[i - 1, j] + n_V1V2p[i - 1, j] + n_V1pV2p[i - 1, j] -  # Ageing in (from double vaccinated, single vaccinated by not protected, and single vaccinated and protected)
  n_V2pV2p[i, j] -                                               # Ageing out
  n_v2pE[i, j]                                                   # Exposed



# Last age groups: Add ageing population
update(S[N_age,]) <- S[N_age, j] +                # Initial
  n_SS[len_ageing, j] -                           # Ageing in
  n_SEs[N_age, j]                                 # Exposed

update(V1[N_age,]) <- V1[N_age, j] +              # Initial
  n_V1V1[len_ageing, j] + n_SV1[len_ageing, j] -  # Ageing in
  n_v1E[N_age, j]                                 # Exposed

update(V2[N_age,]) <- V2[N_age, j] +              # Initial
  n_V2V2[len_ageing, j] + n_V1V2[len_ageing, j] - # Ageing in
  n_v2E[N_age, j]                                 # Exposed

update(V1p[N_age,]) <- V1p[N_age, j] +              # Initial
  n_V1pV1p[len_ageing, j] + n_SV1p[len_ageing, j] - # Ageing in
  n_v1pE[N_age, j]                                  # Exposed

update(V2p[N_age,]) <- V2p[N_age, j] +                                        # Initial
  n_V2pV2p[len_ageing, j] + n_V1V2p[len_ageing, j] + n_V1pV2p[len_ageing, j] -# Ageing in
  n_v2pE[N_age, j]                                                            # Exposed


## Exposed compartments: No ageing (simplification since course of infection for measles is short)
# New value of ES, EV1, and EV2 is: 
#   (old value) + (new_exposed) - (nb of exposed individuals moving to infectious)
update(Es[,]) <-  Es[i, j] +   # Initial
  n_SEs[i, j] -                # New exposed
  n_EsIs[i, j]                 # Moving to infected
update(Ev1[,]) <- Ev1[i, j] +  # Initial
  n_v1E[i, j] + n_v1pE[i,j] -  # New exposed
  n_Ev1Iv1[i, j]               # Moving to infected
update(Ev2[,]) <- Ev2[i, j] +  # Initial
  n_v2E[i, j] + n_v2pE[i,j] -  # New exposed
  n_Ev12v2[i, j]               # Moving to infected

## Infected compartments: No ageing (simplification since course of infection for measles is short)
# New value of IS is:
#   (old value) + (new_infected) + (new_imports) - (nb of infectious moving to recovered)
# import_t corresponds to the number of new imports in this region / age group
# Assumption: All imports are unvaccinated
update(Is[,]) <- Is[i, j] +       # Initial
  n_EsIs[i, j] + import_t[i,j] -  # New infected
  n_IsR[i, j]                     # Moving to recovered
update(Iv1[,]) <- Iv1[i, j] +     # Initial
  n_Ev1Iv1[i, j] -                # New infected
  n_Iv1R[i, j]                    # Moving to recovered
update(Iv2[,]) <- Iv2[i, j] +     # Initial
  n_Ev12v2[i, j] -                # New infected
  n_Iv2R[i, j]                    # Moving to recovered

## Recovered compartments: similar to Susceptible compartments
# First age group:
update(R[1, ]) <- R[1, j] +      # Initial
  n_IsR[1, j] -                  # New recovered
  n_R[1, j]                      # Ageing out
update(RV1[1, ]) <- RV1[1, j] +  # Initial
  n_Iv1R[1, j] -                 # New recovered
  n_RV1[1, j]                    # Ageing out
update(RV2[1, ]) <- RV2[1, j] +  # Initial
  n_Iv2R[1, j] -                 # New recovered
  n_RV2RV2[1, j]                 # Ageing out

# Following age groups:
update(R[2 : len_ageing,]) <- R[i, j] +         # Initial
  n_RR[i - 1, j] -                              # Ageing in
  n_R[i, j] +                                   # Ageing out
  n_IsR[i, j]                                   # Recovered
update(RV1[2 : len_ageing,]) <- RV1[i, j] +     # Initial
  n_RV1RV1[i - 1, j] + n_RRV1[i - 1, j] -       # Ageing in (from single vaccinated or unvaccinated who gained vaccination)
  n_RV1[i, j] +                                 # Ageing out
  n_Iv1R[i, j]                                  # Recovered
update(RV2[2 : len_ageing,]) <- RV2[i, j] +     # Initial
  n_RV2RV2[i - 1, j] + n_RV1RV2[i - 1, j] -     # Ageing in (from double vaccinated or single vaccinated who gained vaccination)
  n_RV2RV2[i, j] +                              # Ageing out
  n_Iv2R[i, j]                                  # Recovered

# Last age group (Remove death here):
update(R[N_age,]) <- R[N_age, j] +                     # Initial
  n_RR[len_ageing, j] +                                # Ageing in
  n_IsR[N_age, j] -                                    # Recovered
  new_birth[j]                                         # Deaths
update(RV1[N_age,]) <- RV1[N_age, j] +                 # Initial
  n_RV1RV1[len_ageing, j] + n_RRV1[len_ageing, j] +    # Ageing in
  n_Iv1R[N_age, j]                                     # Recovered
update(RV2[N_age,]) <- RV2[N_age, j] +                 # Initial
  n_RV2RV2[len_ageing, j] + n_RV1RV2[len_ageing, j] +  # Ageing in
  n_Iv2R[N_age, j]                                     # Recovered

## Track number of new infections each day per age / region
# Unvaccinated
update(new_IS[, ]) <- n_EsIs[i, j]
# Single vaccinated
update(new_IV1[, ]) <- n_Ev1Iv1[i, j]
# Double vaccinated
update(new_IV2[, ]) <- n_Ev12v2[i, j]
# Total number of single vaccinated
update(new_Iv1_tot) <- sum(n_Ev1Iv1[, ])
# Total number of double vaccinated
update(new_Iv2_tot) <- sum(n_Ev12v2[, ])

#### 2- Compute the number of individuals changing between compartments via infection ####
## Draw from binomial distributions: 
# Depends on the number of individuals in the compartment and the probability of transition
# p_XY represents the probability of moving from compartment X to Y (compute in section 3) 
## From maternal immunity to susceptible
n_MS[] <- rbinom(M[1, i], p_MS[i])

## From susceptible (S, V1, V2, V1p, V2p) to exposed
n_SEs[,] <- rbinom(S[i, j], p_SE[i, j])
n_v1E[,] <- rbinom(V1[i, j], p_SE[i, j])
n_v2E[,] <- rbinom(V2[i, j], p_SE[i, j])
n_v1pE[,] <- rbinom(V1p[i, j], p_SEv1[i, j])
n_v2pE[,] <- rbinom(V2p[i, j], p_SEv2[i, j])

## From Exposed to Infected
n_EsIs[,] <- rbinom(Es[i, j], p_EI)
n_Ev1Iv1[,] <- rbinom(Ev1[i, j], p_EI)
n_Ev12v2[,] <- rbinom(Ev2[i, j], p_EI)

## From Infected to Recovered
n_IsR[,] <- rbinom(Is[i, j], p_IR)
n_Iv1R[,] <- rbinom(Iv1[i, j], p_IR)
n_Iv2R[,] <- rbinom(Iv2[i, j], p_IR)



#### 3- Compute the probabilities of transition between each compartment: ####
## From S to E depends on lambda_t, the force of infection (per age / region / date)
p_SE[,] <- 1 - exp(-lambda_t[i, j] * dt) # S to E
## From v1 and v2 to E also depends on lambda, and v2, the protection from infection 
## brought by the vaccine 
p_SEv1[,] <- 1 - exp(-lambda_t[i, j] * v2[i] * dt) # V1 to Ev1
p_SEv2[,] <- 1 - exp(-lambda_t[i, j] * v2[i] * dt) # V2 to Ev2
## From M to S depends on the duration of the maternal immunity
p_MS[] <- 1 - exp(-delta * dt)
## From E to I depends on the duration of the latent period
p_EI <- 1 - exp(-alpha * dt) # E to I
## From I to R depends on the duration of the infectious period
p_IR <- 1 - exp(-gamma * dt) # I to R

## Define the duration of  the latent / infectious periods
alpha <- user(1/11)
gamma <- user(1/8)
delta <- user(1/120)

## Define the default protection from infection brought by the vaccine
## v_fail => Primary vaccine failure (move from S to V1 rather than V1p)
v_fail <- user(.02)
## v_leak => Loss of protection with age
v_leak <- user(0)
## v_sec: baseline risk of secondary vaccine failure
v_sec <- user(0)
## curr_year: current year
curr_year <- year_start + trunc(step / 365)

## v2: risk of secondary vaccine failure by age group
dim(v2) <- N_age
# ## If no waning, v2 = v_sec (if v_sec = 0: impossible to be infected in V1p and V2p)
# ## If waning, depends on age and v_leak
dim(prop_after_change) <- N_age
prop_after_change[] <- if(waning == 3 && (sum(year_per_age[1:(i-1)]) + 1) > (curr_year - 2006))
  # if the youngest year of birth in age group i is before 2006, then set 0 (waning starts at 5)
  0 else if(waning == 3 && sum(year_per_age[1:i]) <= (curr_year - 2006))
    # if the oldest year of birth in age group i is after 2006, then set 1 (waning starts at 3)
    2 else if(waning == 3)
      # if the oldest year of birth in age group i is after 2006, and the youngest is after, compute the number of years before 2006
      (curr_year - 2006 - sum(year_per_age[1:(i-1)])) / year_per_age[i] else
        0

v2[] <-
  if(v_leak == 0 || (sum(year_per_age[1:i]) <= (5 - prop_after_change[i])))
    v_sec else if(waning == 1)
      v_sec + v_leak * ((sum(year_per_age[1:(i-1)]) + year_per_age[i] / 2) - 5) else if(waning == 2)
        v_sec + v_leak * (min((sum(year_per_age[1:(i-1)]) + year_per_age[i] / 2) - 5,
                              curr_year - 2000)) else if(waning == 3)
                                v_sec + v_leak * ((sum(year_per_age[1:(i-1)]) + year_per_age[i] / 2) -
                                                    (5 - 2*prop_after_change[i])) else
                                                      v_sec
# v2[] <- if(v_leak == 0 || sum(year_per_age[1:i]) <= 5) v_sec else if(waning == 1)
#   v_sec + v_leak * ((sum(year_per_age[1:(i-1)]) + year_per_age[i] / 2) - 5) else if(waning == 2)
#     v_sec + v_leak * (min((sum(year_per_age[1:(i-1)]) + year_per_age[i] / 2) - 5, 
#                           curr_year - 2000)) else
#                             v_sec


year_start <- user()
waning <- user()

#### 4- Compute lambda_t, the force of infection including the impact of seasonality ####

# lambda_t depends on lambda (per age group / region) and seasonality
lambda_t[,] <- lambda[i,j] * (1 + X * cos(2 * 3.14159 * time / 365.25 + Y))
# Default values of X and Y (seasonality parameters)
X <- user(1)
Y <- user(1)

### Compute the force of infection from beta (the average nb of contacts), 
### N (the total number of inhabitants) and case_ij (the potential of infection
### stemming from the currently infectious cases, depending on the age / region)
lambda[,] <- beta * cases_ij[i, j] / N

## Default values of beta (average number of contacts)
beta <- user()
## Total population size
N <- sum(S[,]) + sum(M[,]) + sum(V1[,]) + sum(V2[,]) + sum(V1p[,]) + sum(V2p[,]) +  sum(Es[,]) + 
  sum(Ev1[,]) + sum(Ev2[,]) + sum(Is[,]) + sum(Iv1[,]) + sum(Iv2[,]) + sum(R[,]) + 
  sum(RV1[,]) + sum(RV2[,])
## Population size stratified by region / age group
N_strat[,] <- S[i,j] + M[i,j] + V1[i,j] + V2[i,j] + V1p[i,j] + V2p[i,j] + Es[i,j] + Ev1[i,j] +
  Ev2[i,j] + Is[i,j] + Iv1[i,j] + Iv2[i,j] + R[i,j] + RV1[i,j] + RV2[i,j]
dim(N_strat) <- c(N_age, N_reg)

## Compute the potential of infection given the current number of infectious cases
## in each age group / region. Transmission from cases A (age i, region j) to 
## B (age k, region l) is m[i, k] * d[j, l]**-a
## with m the age contact matrix, and d the distance matrix, 
## Therefore, we compute cases_ijkl (the number of connection from [i,j] to [k,l])
cases_ijkl[, , , ] <- m[i, k] * d_a[j, l] *
  (Is[i, j] + vacc * Iv1[i, j] + vacc * Iv2[i, j])
## vacc: protection from onwards transmission brought by vaccination 
## d_a: the distance matrix after applying the spatial kernel (as a rate)

d_a[,] <- if(i == j) 1 else 
  (sum(N_strat[,j])^b) * (theta * ((d[i, j] - 1 )^(-1)) * sum(N_strat[,i])^c) / 
  sum(N_strat[,i])
theta <- user() # No default value, has to be defined by the user
c <- user() # No default value, has to be defined by the user

## We sum over k and l to get the total potential for transmission towards 
## individuals in (i,j) 
cases_ij[, ] <- sum(cases_ijkl[, , i, j])

# Dimension of the transmission matrices
dim(cases_ijkl) <- c(N_age, N_reg, N_age, N_reg)
dim(cases_ij) <- c(N_age, N_reg)
# Dimension of the force of infection matrices
dim(lambda) <- c(N_age, N_reg)
dim(lambda_t) <- c(N_age, N_reg)
# Dimension of the distance matrices and age contact matrix
dim(d_a) <- c(N_reg, N_reg)
dim(m) <- c(N_age, N_age)
dim(d) <- c(N_reg, N_reg)

# Default value of the matrices and parameters
m[, ] <- user() # No default value, has to be defined by the user
d[, ] <- user() # No default value, has to be defined by the user
b <- user(0)
N_age <- user(2)
N_reg <- user(2)
vacc <- user(1)



#### 5- Compute the number of importations ####

## find the current day 
row_number <- trunc(1 + (time / 365.25))
## Draw the number of importations from a poisson distribution using import 
## seasonality and average number of importation
import_t[,] <- rpois(
  mean_import[row_number, j] * (12 * N_strat[i,j])/sum(N_strat[,j]) *
    (1 + X_import * cos(2 * 3.14159 * time / 365.25 + Y_import)))

X_import <- user()
Y_import <- user()

## Define time step
dt <- user(1)
initial(time) <- 1
update(time) <- (step + 1) * dt
initial(iter) <- 1
update(iter) <- step + 1

## Initialise the number of importations / dimensions of the matrices
mean_import[,] <- user()
dim(mean_import) <- c(N_year, N_reg)
dim(import_t) <- c(N_age, N_reg)
N_time <- user(2)
N_year <- N_time/ 365


#### 6- Compute the number of individuals ageing per region / age ####

## Draw number of new births and deaths per region
new_birth[] <- rpois(array_new[i, iter])
# array_new: number of births per region / day
array_new[,] <- user()
dim(array_new) <- c(N_reg, N_time)
dim(new_birth) <- c(N_reg)

## Compute number of individuals ageing in each compartment
# First, compute the number of individuals per vaccination status (S, V1, and V2)
pop_per_age_s[,] <- (S[i, j] + Es[i, j] + Is[i, j] + R[i, j])
pop_per_age_v1[,] <- (V1[i, j] + Ev1[i, j] + Iv1[i, j] + RV1[i, j] + V1p[i, j])
pop_per_age_v2[,] <- (V2[i, j] + Ev2[i, j] + Iv2[i, j] + RV2[i, j] + V2p[i, j])

dim(pop_per_age_s) <- c(N_age, N_reg)
dim(pop_per_age_v1) <- c(N_age, N_reg)
dim(pop_per_age_v2) <- c(N_age, N_reg)
year_per_age[] <- user()
dim(year_per_age) <- c(N_age)

# Draw the number of individuals ageing each day (depending on the population per vaccination status)
# It changes when iter %% 365 == 1, ie when there's the beginning of a new year
N_ageing_S[,] <- if(iter %% 365 == 1 && i != 1) 
  (pop_per_age_s[i, j])/(year_per_age[i] * 365) else 
    if(iter %% 365 == 1 && i == 1)  # in 0-1 yo children: take into account the duration of maternal immunity
      (pop_per_age_s[i, j])/(year_per_age[i] * 365) / (1 - 1/(delta*365)) else 
        N_ageing_S[i, j]
N_ageing_V1[,] <- if(iter %% 365 == 1) (pop_per_age_v1[i, j])/(year_per_age[i] * 365) else N_ageing_V1[i, j]
N_ageing_V2[,] <- if(iter %% 365 == 1) (pop_per_age_v2[i, j])/(year_per_age[i] * 365) else N_ageing_V2[i, j]

dim(N_ageing_S) <- c(N_age, N_reg)
dim(N_ageing_V1) <- c(N_age, N_reg)
dim(N_ageing_V2) <- c(N_age, N_reg)

# Compute the mean number of individuals ageing per compartment using the proportion
# of individuals from pop_per_age that belong to the compartment of interest
# i.e: number of people ageing in compartment 
# X = number of people ageing in this vaccination status *
#         proportion of the population in this vaccination status that belongs to X
mean_S[,] <- if(pop_per_age_s[i, j] > 0) 
  N_ageing_S[i, j] * (S[i, j]/(pop_per_age_s[i, j])) else 0
mean_V1[,] <- if(pop_per_age_v1[i, j] > 0) 
  N_ageing_V1[i, j] * (V1[i, j]/(pop_per_age_v1[i, j])) else 0
mean_V2[,] <- if(pop_per_age_v2[i, j] > 0) 
  N_ageing_V2[i, j] * (V2[i, j]/(pop_per_age_v2[i, j])) else 0
mean_R[,] <- if(pop_per_age_s[i, j] > 0) 
  N_ageing_S[i, j] * (R[i, j]/(pop_per_age_s[i, j])) else 0
mean_RV1[,] <- if(pop_per_age_v1[i, j] > 0) 
  N_ageing_V1[i, j] * (RV1[i, j]/(pop_per_age_v1[i, j])) else 0
mean_RV2[,] <- if(pop_per_age_v2[i, j] > 0) 
  N_ageing_V2[i, j] * (RV2[i, j]/(pop_per_age_v2[i, j])) else 0

mean_V1p[,] <- if(pop_per_age_v1[i, j] > 0) 
  N_ageing_V1[i, j] * (V1p[i, j]/(pop_per_age_v1[i, j])) else 0
mean_V2p[,] <- if(pop_per_age_v2[i, j] > 0) 
  N_ageing_V2[i, j] * (V2p[i, j]/(pop_per_age_v2[i, j])) else 0

## Compute overall number of movements between compartments using a poisson distribution
n_S[,] <- if(mean_S[i, j] < 0 ) 0  else if(mean_S[i, j] < -1) rpois(-1) else  
  rpois(mean_S[i, j])
n_V1[,] <- if(mean_V1[i, j] < 0 ) 0  else if(mean_V1[i, j] < -1) rpois(-1) else  
  rpois(mean_V1[i, j])
n_V2V2[,] <- if(mean_V2[i, j] < 0) 0 else if(mean_V2[i, j] < -1)  rpois(-1) else 
  rpois(mean_V2[i, j])
n_R[,] <- if(mean_R[i, j] < 0) 0 else if(mean_R[i, j] < -1)  rpois(-1) else 
  rpois(mean_R[i, j])
n_RV1[,] <- if(mean_RV1[i, j] < 0) 0 else if(mean_RV1[i, j] < -1)  rpois(-1) else 
  rpois(mean_RV1[i, j])
n_RV2RV2[,] <- if(mean_RV2[i, j] < 0) 0 else if(mean_RV2[i, j] < -1)  rpois(-1) else 
  rpois(mean_RV2[i, j])
n_V1p[,] <- if(mean_V1p[i, j] < 0) 0 else if(mean_V1p[i, j] < -1) rpois(-1) else 
  rpois(mean_V1p[i, j])
n_V2pV2p[,] <- if(mean_V2p[i, j] < 0) 0 else if(mean_V2p[i, j] < -1) rpois(-1) else 
  rpois(mean_V2p[i, j])
len_ageing <- N_age - 1

# If number of people ageing is above the number of people left in the compartment, set it 
# to the current number of individuals left in the compartments
# The number of people left in the compartment is computed from the initial number of people and 
# the number of new exposed / recovered
n_S[,] <- if(n_S[i,j] >= (S[i, j] - n_SEs[i, j])) (S[i, j]- n_SEs[i, j]) else n_S[i,j]
n_V1[,] <- if(n_V1[i,j] >= (V1[i, j] - n_v1E[i, j])) (V1[i, j] - n_v1E[i, j]) else n_V1[i,j]
n_V2V2[,] <- if(n_V2V2[i,j] >= (V2[i, j] - n_v2E[i, j])) (V2[i, j] - n_v2E[i, j]) else n_V2V2[i,j]
n_R[,] <- if(n_R[i,j] >= (R[i, j])) R[i, j] else n_R[i,j]
n_RV1[,] <- if(n_RV1[i,j] >= RV1[i, j]) RV1[i, j] else n_RV1[i,j]
n_RV2RV2[,] <- if(n_RV2RV2[i,j] >= RV2[i, j]) RV2[i, j] else n_RV2RV2[i,j]

n_V1p[,] <- if(n_V1p[i,j] >= (V1p[i, j] - n_v1pE[i, j])) (V1p[i, j] - n_v1pE[i, j]) else n_V1p[i,j]
n_V2pV2p[,] <- if(n_V2pV2p[i,j] >= (V2p[i, j] - n_v2pE[i, j])) (V2p[i, j] - n_v2pE[i, j]) else n_V2pV2p[i,j]

## Compute the proportion of individuals that get vaccinated as they age
# number of people in S who move to V1 = (1 - susceptible who do not get vaccinated)
#         = 1 - ((1 - v1 - v2)[at age i] / (1 - v1 - v2)[at age i - 1])
#  i.e (1 - proportion of unvaccinated at i-1 who are not getting vaccinated at i) 
prop_v1[1:len_ageing,] <- if(iter == 1 || i == 1) array_cov1[i, j, iter] else 
  if(iter <= 365)
    (1 - (1 - array_cov1[i, j, iter] - array_cov2[i, j, iter]) / 
       (1 - array_cov1[i - 1, j, 1] - array_cov2[i - 1, j, 1])) else
         (1 - (1 - array_cov1[i, j, iter] - array_cov2[i, j, iter]) / 
            (1 - array_cov1[i - 1, j, iter - 365] - array_cov2[i - 1, j, iter - 365]))
prop_v1[,] <- if(prop_v1[i,j] < 0)  0 else 
  if(prop_v1[i,j] > 1) 1 else 
    prop_v1[i, j]

dim(prop_v1) <- c(len_ageing, N_reg)
# Draw the number of individuals gaining vaccination as they age
n_SVtot[1 : len_ageing,] <- rbinom(n_S[i, j], prop_v1[i, j])
# Compute the number of primary vaccine failure
n_SV1[1 : len_ageing,] <- rbinom(n_SVtot[i, j], v_fail)
n_SV1p[1 : len_ageing,] <- n_SVtot[i, j] - n_SV1[i, j]

n_RRV1[1 : len_ageing,] <- rbinom(n_R[i, j], prop_v1[i, j])
# Draw the number of individuals who remain susceptible
n_SS[1 : len_ageing,] <-  n_S[i, j] - n_SVtot[i, j]
n_RR[1 : len_ageing,] <- n_R[i, j] - n_RRV1[i, j]

# number of people in V1 who move to V2 = V1 who get vaccinated
#         = (v2[at age i] - V2[at age i-1]) / (v1[at age i - 1])
#  i.e (1 - proportion of V1 at i-1 who are getting vaccinated at i) 
prop_v1v2[1:len_ageing,] <- if(iter == 1 || i == 1) 0 else if(iter <= 365)
  (array_cov2[i, j, iter] - array_cov2[i - 1, j, 1]) / 
  (array_cov1[i - 1, j, 1]) else
    (array_cov2[i, j, iter] - array_cov2[i - 1, j, iter - 365]) / 
  (array_cov1[i - 1, j, iter - 365])

prop_v1v2[1:len_ageing,] <- 
  if(iter <= 365 && i >= 1 && array_cov1[i - 1, j, 1] == 0) 0 else 
    if(iter > 365 && i >= 1 && array_cov1[i - 1, j, iter - 365] == 0) 0 else
      prop_v1v2[i,j]

prop_v1v2[,] <- if(prop_v1v2[i,j] < 0)  0 else if(prop_v1v2[i,j] > 1) 1 else 
  prop_v1v2[i, j]
dim(prop_v1v2) <- c(len_ageing, N_reg)
# Draw the number of individuals gaining vaccination as they age 
# number of individuals going from single vaccinated & protected to double vaccinated & protected
n_V1pV2p[1 : len_ageing,] <- rbinom(n_V1p[i, j], prop_v1v2[i, j])

# number of individuals going from primary vaccine failure to double vaccinated
n_V1V2tot[1 : len_ageing,] <- rbinom(n_V1[i, j], prop_v1v2[i, j])
# number of double primary vaccine failure
n_V1V2[1 : len_ageing,] <- rbinom(n_V1V2tot[i, j], v_fail)
n_V1V2p[1 : len_ageing,] <- n_V1V2tot[i, j] - n_V1V2[i, j]

n_RV1RV2[1 : len_ageing,] <- rbinom(n_RV1[i, j], prop_v1v2[i, j])
# Draw the number of individuals who remain in V1
n_V1V1[1 : len_ageing,] <- n_V1[i, j] - n_V1V2tot[i, j]
n_V1pV1p[1 : len_ageing,] <- n_V1p[i, j] - n_V1pV2p[i, j]
n_RV1RV1[1 : len_ageing,] <- n_RV1[i, j] - n_RV1RV2[i, j]

array_cov1[,,] <- user()
dim(array_cov1) <- c(len_ageing, N_reg, N_time)
array_cov2[,,] <- user()
dim(array_cov2) <- c(len_ageing, N_reg, N_time)


#### 7- Initial conditions ####

## Compute the initial number of susceptibles from S_init and the duration of maternal immunity
initial(S[,]) <- if(i == 1) round(S_init[i, j] * (1 - 1/(delta * 365))) else
  round(S_init[i, j] * (1 - recov[i,j]))
initial(M[,]) <- if(i == 1) round(S_init[i, j] * 1/(delta * 365)) else 0

## Use v_fail to compute the number of primary vaccine failure at t=0
initial(V1[,]) <- round(V1_init[i, j] * v_fail * (1 - recov[i,j]))
initial(V2[,]) <- round(V2_init[i, j] * (v_fail ** 2) * (1 - recov[i,j]))
initial(V1p[,]) <- round(V1_init[i, j] * (1 - v_fail))
initial(V2p[,]) <- round(V2_init[i, j] * (1 - v_fail**2))

## Use v_fail and recov to compute the number of primary vaccine failure and recovered at t=0
initial(R[,]) <- round(S_init[i,j] * recov[i,j])
initial(RV1[,]) <- round(V1_init[i, j] * v_fail * (recov[i,j]))
initial(RV2[,]) <- round(V2_init[i, j] * (v_fail ** 2) * (recov[i,j]))


initial(Es[,]) <- Es_init[i, j]
initial(Ev1[,]) <- Ev1_init[i, j]
initial(Ev2[,]) <- Ev2_init[i, j]
initial(Is[,]) <- Is_init[i, j]
initial(Iv1[,]) <- Iv1_init[i, j]
initial(Iv2[,]) <- Iv2_init[i, j]
initial(new_IS[, ]) <- 0
initial(new_IV1[, ]) <- 0
initial(new_IV2[, ]) <- 0
initial(new_Iv1_tot) <- 0
initial(new_Iv2_tot) <- 0
recov[,] <- user()
dim(recov) <- c(N_age, N_reg)

## Default values
S_init[,] <- user()
V1_init[,] <- user()
V2_init[,] <- user()
Es_init[,] <- user()
Ev1_init[,] <- user()
Ev2_init[,] <- user()
Is_init[,] <- user()
Iv1_init[,] <- user()
Iv2_init[,] <- user()

### Dimensions of each compartment
dim(S) <- c(N_age, N_reg)
dim(M) <- c(N_age, N_reg)
dim(V1) <- c(N_age, N_reg)
dim(V2) <- c(N_age, N_reg)
dim(V1p) <- c(N_age, N_reg)
dim(V2p) <- c(N_age, N_reg)
dim(Es) <- c(N_age, N_reg)
dim(Ev1) <- c(N_age, N_reg)
dim(Ev2) <- c(N_age, N_reg)
dim(Is) <- c(N_age, N_reg)
dim(Iv1) <- c(N_age, N_reg)
dim(Iv2) <- c(N_age, N_reg)
dim(R) <- c(N_age, N_reg)
dim(RV1) <- c(N_age, N_reg)
dim(RV2) <- c(N_age, N_reg)
dim(new_IS) <- c(N_age, N_reg)
dim(new_IV1) <- c(N_age, N_reg)
dim(new_IV2) <- c(N_age, N_reg)

### Initialise the dimensions of the compartments' initial value
dim(S_init) <- c(N_age, N_reg)
dim(V1_init) <- c(N_age, N_reg)
dim(V2_init) <- c(N_age, N_reg)
dim(Es_init) <- c(N_age, N_reg)
dim(Ev1_init) <- c(N_age, N_reg)
dim(Ev2_init) <- c(N_age, N_reg)
dim(Is_init) <- c(N_age, N_reg)
dim(Iv1_init) <- c(N_age, N_reg)
dim(Iv2_init) <- c(N_age, N_reg)


### Initialise the dimensions of the number of individual moving between the compartments
dim(n_SEs) <- c(N_age, N_reg)
dim(n_v1E) <- c(N_age, N_reg)
dim(n_v2E) <- c(N_age, N_reg)
dim(n_v1pE) <- c(N_age, N_reg)
dim(n_v2pE) <- c(N_age, N_reg)
dim(n_MS) <- c(N_reg)
dim(n_EsIs) <- c(N_age, N_reg)
dim(n_Ev1Iv1) <- c(N_age, N_reg)
dim(n_Ev12v2) <- c(N_age, N_reg)
dim(n_IsR) <- c(N_age, N_reg)
dim(n_Iv1R) <- c(N_age, N_reg)
dim(n_Iv2R) <- c(N_age, N_reg)
dim(p_SE) <- c(N_age, N_reg)
dim(p_MS) <- c(N_reg)
dim(p_SEv1) <- c(N_age, N_reg)
dim(p_SEv2) <- c(N_age, N_reg)

dim(n_SS) <- c(N_age, N_reg)
dim(n_SV1) <- c(N_age, N_reg)
dim(n_SVtot) <- c(N_age, N_reg)
dim(n_V1V2tot) <- c(N_age, N_reg)
dim(n_V1V1) <- c(N_age, N_reg)
dim(n_V1V2) <- c(N_age, N_reg)
dim(n_V2V2) <- c(N_age, N_reg)
dim(n_SV1p) <- c(N_age, N_reg)
dim(n_V1pV1p) <- c(N_age, N_reg)
dim(n_V1V2p) <- c(N_age, N_reg)
dim(n_V1pV2p) <- c(N_age, N_reg)
dim(n_RR) <- c(N_age, N_reg)
dim(n_RRV1) <- c(N_age, N_reg)
dim(n_RV1RV2) <- c(N_age, N_reg)
dim(n_RV1RV1) <- c(N_age, N_reg)
dim(n_RV2RV2) <- c(N_age, N_reg)
dim(n_S) <- c(N_age, N_reg)
dim(n_V1) <- c(N_age, N_reg)
dim(n_V1p) <- c(N_age, N_reg)
dim(n_V2pV2p) <- c(N_age, N_reg)
dim(n_R) <- c(N_age, N_reg)
dim(n_RV1) <- c(N_age, N_reg)
dim(mean_S) <- c(N_age, N_reg)
dim(mean_V1) <- c(N_age, N_reg)
dim(mean_V2) <- c(N_age, N_reg)
dim(mean_V1p) <- c(N_age, N_reg)
dim(mean_V2p) <- c(N_age, N_reg)
dim(mean_R) <- c(N_age, N_reg)
dim(mean_RV1) <- c(N_age, N_reg)
dim(mean_RV2) <- c(N_age, N_reg)
