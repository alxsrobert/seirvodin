#' Simulated outbreaks
#'
#' We generated a simulated dataset  to illustrate \code{seirvodin}. 
#' It contains 1,550 cases, from simulated outbreaks, stratified
#' by age group, region, onset date and vaccination status. 
#' The dataset contains the following column:
#'
#' \itemize{
#'
#' \item \code{date}: The onset date (from 1 to 1825).
#'
#' \item \code{region}: The region of origin ("A" or "B").
#' 
#' \item \code{age_groups}: The age group ("0-2", "2-5", "5-10", or ">10"). 
#' 
#' \item \code{vaccinated}: The vaccinated status (new_IS = unvaccinated, new_IV1 = one-dose recipient, new_IV2 = two-dose recipient)
#' 
#' \item \code{N}: The number of cases reported on that date, in this region / age group/ vaccine status.
#' 
#'
#' }
#'
#' @aliases simulated_outbreak
#' @docType data
#' @author Alexis Robert \email{alexis.robert@lshtm.ac.uk}
#' @keywords datasets
#'
#' @examples
#' data("simulated_outbreak")
#' simulated_outbreak
#' 
"simulated_outbreak"
