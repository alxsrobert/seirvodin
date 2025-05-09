% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/function_fit.R
\name{specs_run}
\alias{specs_run}
\title{Set specifications of the model run}
\usage{
specs_run(N_year, year_start, n_steps, waning, alpha = 11, gamma = 8)
}
\arguments{
\item{N_year}{Number of years the simulations run for.}

\item{year_start}{starting year of the simulations.}

\item{n_steps}{Number of MCMC steps.}

\item{waning}{Character values, corresponds to whether waning is included in
the model, expect one of three values: "no", "since_vax", or "since_eli"}

\item{alpha}{Duration of incubation period.}

\item{gamma}{Duration of infectious period.}
}
\value{
List of specifications
}
\description{
Set specifications of the model run
}
