# seirvodin

`seirvodin` presents a compartmental moodel implemented in `odin_dust`, this model is stratified by age, region and vaccination status. The model was used in [Robert et al, Lancet Public Health, 2024](https://www.thelancet.com/journals/lanpub/article/PIIS2468-2667(24)00181-6/fulltext). The structure of the model is presented in the figure below:

![Structure of the model for a given age and region](https://github.com/user-attachments/assets/b3ac8898-1ef3-4c6e-9763-a18004967216) (From [Robert et al, Lancet Public Health, 2024](https://www.thelancet.com/journals/lanpub/article/PIIS2468-2667(24)00181-6/fulltext))

`seirvodin` can be used to fit the model to case data, or generate stochastic (or deterministic) simulations. It was designed to analyse measles outbreaks, and therefore expects routine coverage data (for first and second dose). It can estimate the coverage of catchup campaigns prior the start of fitting (or simulations), but does not explicitly allow catchup campaigns during the fitting period (but the coverage of catchup campaigns could be integrated in the routine coverage data).

Installation
-------------

To install the development version from github (requires Rtools on windows and GSL headers on all platforms):

```{r, eval = FALSE}
devtools::install_github("alxsrobert/seirvodin")
```

Documentation
-------------

Five functions from `seirvodin` are exported:
* `seirv_age_region`: `dust` compartmental model stratified by age and region.
* `run_model()`: Function to fit the model to case data.
* `generate_outbreaks()`: Function to generate stochastic or deterministic simulations from the parameter set estimated during the model fitting process.
* `specs_run()`: Function to generate the specifications of the model fitting.
* `specs_simulations()`: Function to generate the specifications of the simulations.

The `seirv_age_region` model expects at least 9 parameters, and has several optional parameters:
* `beta`: Infection rate.
* `delta`: Duration of maternal immunity.
* `X`: Seasonality of infection rate.
* `Y`: Seasonality of infection rate.
* `X_import`: Seasonality of importations.
* `Y_import`: Seasonality of importations.
* `v_fail`: Proportion of primary vaccine failure.
* `vacc`: Risk of onward transmission from vaccinated cases, compared to unvaccinated cases.
* `report_import`: Proportion of importations reported
* `b`: (optional) Spatial parameter.
* `c`: (optional) Spatial parameter.
* `theta`: (optional) Spatial parameter.
* `v_leak`: (optional) Waning of immunity per year (set to 0 if not set by the user).
* `v_sec`: (optional) Baseline risk of secondary vaccine failure (set to 0 if not set by the user).
* `recov_X`: (optional) Proportion of individuals starting as recovered in a given age group, X corresponds to the age group targeted (e.g. `recov_7` corresponds to the proportion of individuals in the seventh age group who will start as recovered).
* `catchup_X`: (optional) Proportion of individuals vaccinated during a catch-up campaign (moving from V1 to V2) in a given age group, X corresponds to the age group targeted (e.g. `catchup_7` corresponds to the proportion of individuals (in V1) in the seventh age group targeted by the catch-up campaign who will move to V2).
* `catchup2_X`: (optional) Proportion of individuals vaccinated during a catch-up campaign (moving from S and V1 to V2) in a given age group, X corresponds to the age group targeted (e.g. `catchup2_7` corresponds to the proportion of individuals (in V1 or S) in the seventh age group targeted by the catch-up campaign who will move to V2).
* `v_X`: (optional) Proportion of individuals vaccinated (moving from S to V1) in a given age group, X corresponds to the age group targeted (e.g. `v_7` corresponds to the proportion of individuals (in S) in the seventh age group targeted by the catch-up campaign who will move to V1).



