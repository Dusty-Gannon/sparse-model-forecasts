# Analysis code

* `AICvsSTANfigures.R`: This R script generates figures comparing the performance of stepwise AIC model selection and Bayesian variable selection using the horseshoe prior (via STAN) in time series regression settings. It includes:
    - Figure 1: A simulated example illustrating prediction intervals and coefficient recovery for both AIC and STAN-based models.

    - Figure 2: Violin plots comparing true positive rate (TPR), true negative rate (TNR), and RMSE across varying levels of predictor correlation.

    - Figure 3: Violin plots comparing prediction RMSE under covariate shift and under assumptions of stationarity.

The script relies on precomputed simulation results that can be replicated using scripts in the `Simulations/` directory. Outputs are saved as PDF figures in the `Figures/` directory.

* `data_ex.R`: Exploratory plots for Prism and tree ring data.

* `datatrends.R`: Exploratory plots and analyses for Prism and tree ring data.

* `plot_individual_ARp_simulation.R`: Script to construct Figure 1 of the manuscript.

* `PrismDat_Agg.R`: Script to aggregate the prism data into different parts of the water-year.

* `tree_dat_forecasts_env_covariates.R`: Script to fit Bayesian sparse models and conduct AIC-based model selection to predict tree growth based on relationships with historical Prism data. This script is used to create **Figures 7 and 8** of the main text.

* `tree_dat_forecasts_fourier.R`: Script to compare the proposed method to the Auto-Arima method of Hyndman and Khandakar, 2008. 




