# Population-Forecast
Probabilistic population forecasts for small regions using Bayesian hierarchical models.
Regional population is forecasted by combining separate Bayesian models for mortality, fertility, and migration. Each component is estimated using Stan and integrated into a cohort-component (Leslie matrix) projection framework. 

## Models
There are a collection of Bayesian hierarchical models found in `StanCode` that were used for the analysis.

All models are fitted using **Stan**. Information on how to install Stan can be found below.


## Data
Within the `Data` folder the following can be found: Deaths, Population, Birth, and Migration data by Age, Sex, Time and Region which was downloaded from <https://www.statistikdaten.bayern.de/genesis/online/logon> the database of the Bavarian statistical institute.
Also within the subfolder `Maps` are shapefiles which are needed for the creation of maps.

## R Code
There are multiple R files each requiring different packages to run.

* `01_LoadingData.R` Loads all raw data from Excel, transforms to long format, calculates person-years of exposure and imputes missing migration values. Finished data for analysis is saved as *TotalData.RData* which can also be found in the `Data` folder.
* `02_Functions.R` Includes all self-written functions needed for demographic calculations (Leslie matrices, life tables, cohort indices), spatial model preparation (BYM2 scaling factors, CAR data formatting), and visualization (population pyramids, Stan output summaries).
* `03_Mortality.R` Fits the RH_BYM2 Stan model separately for males and females and saves posterior mortality rate forecasts.
* `04_Fertility.R` Fits the Lee-Carter Poisson Stan model to age-specific fertility rates and saves posterior forecasts. Includes out-of-sample evaluation.
* `05_Migration.R` Fits the Azoze-Raftery skew-normal model for net migration counts and the Dirichlet model for age-specific migration schedules. Saves posterior forecasts for both components.
* `06_PopulationProjection.R` Combines mortality, fertility, and migration forecasts within a Leslie matrix projection to generate probabilistic population trajectories for 2024–2044. Includes out-of-sample validation for 2020–2024.
* `07_Graphics.R` Script for all visualization plots, including regional population projections with credible intervals, population pyramids, maps of projected change, and migration schedules.

The modelling part requires the following packages:
`rstan`, `tidyverse`, `sf`, `spdep`, `lme4`, `mice`

The analysis and projection part needs the following packages (in addition to the above):
`reshape2`

The visualization part needs:
`ggplot2`, `ggdist`, `cowplot`, `paletteer`

For **Stan** to run on the computer, one needs an interface to R called `rstan`. To install `rstan` see the following link: <https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started>
