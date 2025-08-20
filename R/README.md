# User-defined functions

**All functions are documented using Roxygen2. For documentation, see `/man/`.**

* `AR-p_sim_functions.R`: Functions to simulate AR-p datasets and extract results from stan model objects.

* `AR-p_fit_functions.R`: Functions to fit sparse models to simulated data.

* `continue_fourier.R`: Function to continue a set of Fourier basis functions into the future.

* `lag_covariates.R`: Function to create a set of covariates lagged by $\ell$ time steps and combine them into a dataframe.

* `model_checking.R`: Functions used to check model fit and predictions.

* `sim1_temporal_drivers`: Functions to simulate realistic temporal covariates with seasonal cycles and trends.

* `simulate_seasonal_ARp.R`: Simulate a seasonal AR timeseries by translating a seasonal cycle into a stationary vector of phis. Includes the option to add a sparse covariate matrix. 

* `simulation_processing_functions.R`: Functions to summarize the results from fitted AR-p simulation models.

* `tau0_from_data.R`: Function to provide an initial $\tau_0$ based on Piironen and Vehtari, 2017.

