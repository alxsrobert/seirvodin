##' seirv_age_region model. This is a dust model.
##' @name seirv_age_region
##' @title Age and region structured SEIRV model
##' @param beta: Infection rate.
##' @param delta: Duration of maternal immunity.
##' @param X: Seasonality of infection rate.
##' @param Y: Seasonality of infection rate.
##' @param X_import: Seasonality of importations.
##' @param Y_import: Seasonality of importations.
##' @param v_fail: Proportion of primary vaccine failure.
##' @param vacc: Risk of onward transmission from vaccinated cases, compared to unvaccinated cases.
##' @param report_import: Proportion of importations reported.
##' @param b: (optional) Spatial parameter.
##' @param c: (optional) Spatial parameter.
##' @param theta: (optional) Spatial parameter.
##' @param v_leak: (optional) Waning of immunity per year (set to 0 if not set by the user).
##' @param v_sec: (optional) Baseline risk of secondary vaccine failure (set to 0 if not set by the user).
##' @param recov_X: (optional) Proportion of individuals starting as recovered in a given age group, X corresponds to the age group targeted (e.g. recov_7 corresponds to the proportion of individuals in the seventh age group who will start as recovered).
##' @param catchup_X: (optional) Proportion of individuals vaccinated during a catch-up campaign (moving from V1 to V2) in a given age group, X corresponds to the age group targeted (e.g. catchup_7 corresponds to the proportion of individuals (in V1) in the seventh age group targeted by the catch-up campaign who will move to V2).
##' @param catchup2_X: (optional) Proportion of individuals vaccinated during a catch-up campaign (moving from S and V1 to V2) in a given age group, X corresponds to the age group targeted (e.g. catchup2_7 corresponds to the proportion of individuals (in V1 or S) in the seventh age group targeted by the catch-up campaign who will move to V2).
##' @param v_X: (optional) Proportion of individuals vaccinated (moving from S to V1) in a given age group, X corresponds to the age group targeted (e.g. v_7 corresponds to the proportion of individuals (in S) in the seventh age group targeted by the catch-up campaign who will move to V1).##' @export seirv_age_region
##' @param alpha: (optional) Duration of incubation period (default to 11 days).
##' @param gamma: (optional) Duration of infectious period (default to 8 days).
##' 
##' @import odin.dust
##' 
##' @export seirv_age_region
NULL
