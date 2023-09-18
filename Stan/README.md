# Directory contents

-   `AR-p_err3_FHS_DG.stan`: Stan model to fit a Gaussian dynamic regression model with AR error structure. 'Finnish Horseshoe priors' are placed on the AR coefficients and the regression coefficients, each with separate global shrinkage priors so that model selection of the regression coefficient and AR coefficients can happen separately. Weakly-informative priors can be specified for a number of covariate coefficients (including the intercept term) by defining the data input `P_0` $\ge$ 1. These must come as the first columns in the model matrix supplied as data to Stan.

    -   **New feature**: Forecasts are done internally by supplying a matrix of "future" covariates as well as the number of steps to forecast.

-   `AR-p_err3_Gauss_DG.stan`: Stan model to fit a Gaussian dynamic regression model with AR error structure. Weakly-informative priors are used for all parameters, but flat priors can be specified for specific covariate coefficients (including the intercept term) by defining the data input `P_0` $\ge$ 1. These must come as the first columns in the model matrix supplied as data to Stan.

-   `AR-p_err3_Flat_DG.stan`: Stan model to fit a Gaussian dynamic regression model with AR error structure. Improper flat priors are used for all parameters.

**Note**: A script demonstrating how to use and compare the above three models using Fourier components to model seasonality can be found in `/Model_evals/testing_AR-p-err_mods_DG.R`.

-   `sparse_reg_FHS.stan`: Stan model to fit a Bayesian sparse regression using Regularized Horseshoe priors. Weakly-informative priors can be specified for a number of covariate coefficients (including the intercept term) by defining the data input `P_0` $\ge$ 1. These must come as the first columns in the model matrix supplied as data to Stan. It is recommended that `P_0` is at least 1 to leave the intercept term unregularized. Predictions can be constructed by supplying a model matrix, `X_new`, of `N_new` observations. If no predictions are required, set `N_new = 0` and `X_new = matrix(data = 1, nrow = 0, ncol = P)` as a null or dummy matrix in the data input to Stan.

-   `AR-p_FHS-p-beta.stan`: Stan model to fit a Gaussian AR(p) process with arbitrary p. 'Finnish Horseshoe Priors' are placed on the AR coefficients and the regression coefficients to perform Bayesian model selection/regularization on both the phi and beta vectors.

-   `AR-p_FHS.stan`: Stan model to fit a Gaussian AR(p) process with arbitrary p. 'Finnish Horseshoe Priors' are placed on the regression coefficients to perform Bayesian model selection/regularization on the beta vector.

-   `AR-p.stan`: Stan model to fit a Gaussian AR(p) process with arbitrary p. Commonly used weakly-informative priors are placed on the AR coefficients and the regression coefficients.

-   `ARMA-p-q_FHS.stan`: Stan model to fit a Gaussian ARMA(p, q) process with arbitrary p and q. 'Finnish Horseshoe Priors' are placed on the regression coefficients to perform Bayesian model selection/regularization on the beta vector.

# Archived

-   `Gamma_LatAR1_FHS_EBtest.stan`: Gamma-distributed response with 'Finnish Horseshoe Priors' on the regression coefficients to perform model selection. The autocorrelation is introduced through a latent AR process on the link scale. Because of identifiability problems with the dispersion parameter and error of the latent process, an `empirical Bayes` approach is implemented, supplying the MLE for the dispersion parameter in terms of the other parameters in the model. **Note that these Latent variable models are computationally costly with long time series.**

-   `Pois_ACM-p-q_FHS.stan`: Auto-regressive conditional mean model with arbitrary p and q and conditionally Poisson-distributed response. 'Finnish Horseshoe Priors' are placed on the regression coefficients to perform model selection.

-   `Pois_GLARMA-1-1_FHS.stan`: Generalized linear ARMA(1, 1) model with conditionally Poisson-distributed response. 'Finnish Horseshoe Priors' are placed on the regression coefficients to perform model selection.

-   `Pois_LatAR1_FHS.stan`: Conditionally Poisson-distributed response with 'Finnish Horseshoe Priors' on the regression coefficients to perform model selection. The autocorrelation is introduced through a latent AR process on the link scale. **Note that these Latent variable models are computationally costly with long time series.**

-   `Pois_ricker_fixed_lambda_FHS.stan`: Special case of the `Pois_ACM-p-q_FHS.stan` model above, catered to a Ricker community model in which intrinsic growth of the focal species is unknown but assumed fixed through time. 'Finnish Horseshoe Priors' are placed on the competition coefficients to help select important competing species.

-   `Pois_ricker_known_lambda_FHS.stan`: Special case of the `Pois_ACM-p-q_FHS.stan` model above, catered to a Ricker community model in which intrinsic growth of the focal species is assumed known and supplied as data. 'Finnish Horseshoe Priors' are placed on the competition coefficients to help select important competing species.

-   `Pois_ricker_partial_lambda_FHS.stan`: Special case of the `Pois_ACM-p-q_FHS.stan` model above, catered to a Ricker community model in which *maximum intrinsic growth* of the focal species under ideal growing conditions is assumed known, but intrinsic growth in each time point is determined by an environmental covariate. 'Finnish Horseshoe Priors' are placed on the competition coefficients to help select important competing species.
