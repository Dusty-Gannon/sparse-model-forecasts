# Directory contents

-   `AR-p_FHS-p-beta.stan`: Stan model to fit a Gaussian AR(p) process with arbitrary p. 'Finnish Horseshoe Priors' are placed on the AR coefficients and the regression coefficients to perform Bayesian model selection/regularization on both the phi and beta vectors.

-   `AR-p_FHS.stan`: Stan model to fit a Gaussian AR(p) process with arbitrary p. 'Finnish Horseshoe Priors' are placed on the regression coefficients to perform Bayesian model selection/regularization on the beta vector.

-   `AR-p.stan`: Stan model to fit a Gaussian AR(p) process with arbitrary p. Commonly used weakly-informative priors are placed on the AR coefficients and the regression coefficients.

-   `ARMA-p-q_FHS.stan`: Stan model to fit a Gaussian ARMA(p, q) process with arbitrary p and q. 'Finnish Horseshoe Priors' are placed on the regression coefficients to perform Bayesian model selection/regularization on the beta vector.

-   `Gamma_LatAR1_FHS_EBtest.stan`: Gamma-distributed response with 'Finnish Horseshoe Priors' on the regression coefficients to perform model selection. The autocorrelation is introduced through a latent AR process on the link scale. Because of identifiability problems with the dispersion parameter and error of the latent process, an `empirical Bayes` approach is implemented, supplying the MLE for the dispersion parameter in terms of the other parameters in the model. **Note that these Latent variable models are computationally costly with long time series.**

-   `Pois_ACM-p-q_FHS.stan`: Auto-regressive conditional mean model with arbitrary p and q and conditionally Poisson-distributed response. 'Finnish Horseshoe Priors' are placed on the regression coefficients to perform model selection.

-   `Pois_GLARMA-1-1_FHS.stan`: Generalized linear ARMA(1, 1) model with conditionally Poisson-distributed response. 'Finnish Horseshoe Priors' are placed on the regression coefficients to perform model selection.

-   `Pois_LatAR1_FHS.stan`: Conditionally Poisson-distributed response with 'Finnish Horseshoe Priors' on the regression coefficients to perform model selection. The autocorrelation is introduced through a latent AR process on the link scale. **Note that these Latent variable models are computationally costly with long time series.**

-   `Pois_ricker_fixed_lambda_FHS.stan`: Special case of the `Pois_ACM-p-q_FHS.stan` model above, catered to a Ricker community model in which intrinsic growth of the focal species is unknown but assumed fixed through time. 'Finnish Horseshoe Priors' are placed on the competition coefficients to help select important competing species.

-   `Pois_ricker_known_lambda_FHS.stan`: Special case of the `Pois_ACM-p-q_FHS.stan` model above, catered to a Ricker community model in which intrinsic growth of the focal species is assumed known and supplied as data. 'Finnish Horseshoe Priors' are placed on the competition coefficients to help select important competing species.

-   `Pois_ricker_partial_lambda_FHS.stan`: Special case of the `Pois_ACM-p-q_FHS.stan` model above, catered to a Ricker community model in which *maximum intrinsic growth* of the focal species under ideal growing conditions is assumed known, but intrinsic growth in each time point is determined by an environmental covariate. 'Finnish Horseshoe Priors' are placed on the competition coefficients to help select important competing species.
