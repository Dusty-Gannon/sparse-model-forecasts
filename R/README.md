# User-defined functions

**All functions are documented using Roxygen2. For documentation, see `/man/`.**

* `fitting_functions_lnorm_ricker.R`: Wrapper functions used to take in simulated data, fit Ricker models with log-normal errors, and summarize the results.

* `model_checking.R`: Functions used to verify that STAN and JAGS models are written correctly.

* `stat_model_sim_functions.R`: Functions for simulating data from statistical time series models.

* `terrestrial_data_processing_functions.R`: Functions for processing the simulated data from competitive communities (Ricker models).

* `terrestrial_lottery_sim_functions.R`: Functions to simulate a spatial lottery model with competition, dispersal, and mortality across a lattice.

* `terrestrial_ricker_sim_functions.R`: Functions used to simulate time series for competitive communities using Ricker models. These include functions for Poisson-distributed errors and log-normally distributed errors.
