######
# This script simulates conditional Poisson time series with
# a latent AR(1) process to induce autocorrelation and 3 non-zero
# regression coefficients, then fits them using Finnish-Horseshoe
# priors for the regression coefficients to induce sparsity
######

# libraries
  library(rstan)
  library(here)

# load user-defined functions (teton way, not devtools way)
  src_files <- list.files(here("R/"), pattern = "*.R", full.names = T)
  sapply(src_files, source, .GlobalEnv)

# load command-line arguments
  args <- commandArgs(trailingOnly = T)

## lengths of each time series to simulate and fit the model
#  ns <- seq(50, 500, by = 50)

# define parameters
  sim_params <- list(
    phi = 0.6,
    b = c(1, -1, 0.5), # non-zero effects
    P = 50,            # total number of candidate variables (minus the intercept)
    sigma_e = 0.5,     # sd of random innovations
    a = 2              # gamma shape parameter
  )


# function to simulate the data, fit the model, and report stats on model fit and predictive performance
  latent_ar1_FHS <- function(n, sim_params, mod_file, q_cutoff = 0.9, n_test = 30){

    # useful variables
    n_train <- n - n_test

    # create latent AR(1) variables
    lambda <- with(sim_params, {
      lambda <- vector(mode = "double", length = n)
      lambda[1] <- rnorm(1) * (sigma_e/sqrt(1 - phi^2))
      for(t in 2:n){
        lambda[t] <- phi * lambda[t-1] + rnorm(1, sd = sigma_e)
      }
      lambda
    })

    # create coefficient vector
    beta <- c(
      1, sim_params$b,
      rep(0, sim_params$P - length(sim_params$b))
    )

    # create model matrix
    X <- with(sim_params, {
      matrix(
        data = c(rep(1, n), rnorm(n * P)),
        nrow = n,
        ncol = P + 1
      )
    })

    # create means
    mu <- as.double(exp(X %*% beta + lambda))

    # convert back to gamma parameters
    b <- with(sim_params, {a / mu})

    # simulate response vector
    y <- rgamma(n, shape = sim_params$a, rate = b)

    # split data into training and testing data
    y_train <- y[1:n_train]
    y_test <- y[(n_train + 1):n]
    X_train <- X[1:n_train, ]
    X_test <- X[(n_train + 1):n, ]

    # estimate for tau0
    tau0 <- tau0(
      y = y_train,
      m0 = 5,
      M = sim_params$P,
      N = n_test,
      fam = "gamma"
    )

    # data list for model fit
    datlist <- with(sim_params,{
      list(
        N = n_train,
        P0 = 1,
        P = ncol(X) - 1,
        y = y_train,
        X_alpha = as.matrix(X_train[,1]),
        X_beta = X_train[,-1],
        tau0 = tau0,
        slab_scl = 1,
        slab_df = 6
      )
    })

    # compile stan model
    latar1_FHS <- stan_model(mod_file)

    # fit the model
    mfit2 <- sampling(
      latar1_FHS2,
      data = datlist,
      chains = 2,
      iter = 2000,
      #cores = 3,
      control = list(adapt_delta = 0.99, max_treedepth = 15)
    )

    # compute out-of-sample RMSE
    # first extract parameters
    beta_all_post <- cbind(
      rstan::extract(mfit, pars = "alpha")$alpha,
      rstan::extract(mfit, pars = "beta")$beta
    )
    phi_post <- rstan::extract(mfit, pars = "phi")$phi
    sigma_e_post <- rstan::extract(mfit, pars = "sigma_e")$sigma_e
    a_post <- rstan::extract(mfit, pars = "a")$a

    # extract the last lambda, then use it for predicting next y
    lambda_tm1_post <- rstan::extract(mfit, pars = "lambda")$lambda[, n_train]

    # for each draw from the posterior, predict the next n_test time points
    #  for the AR process variables
    lambda_hat <- sapply(
      1:length(sigma_e_post),
      FUN = function(i, phi, sigma_e, lambda_tm1, n_test){
        lambda_pred <- vector(mode = "double", length = n_test)
        # initial value based on lambda_tm1
        lambda_pred[1] <- phi[i] * lambda_tm1[i] + rnorm(1, sd = sigma_e[i])
        for(t in 2:n_test){
          lambda_pred[t] <- phi[i] * lambda_pred[t-1] + rnorm(1, sd = sigma_e[i])
        }
        return(lambda_pred)
      },
      phi = phi_post, sigma_e = sigma_e_post,
      lambda_tm1 = lambda_tm1_post,
      n_test = n_test
    )

    # create the linear predictor for each draw at the next n_test time points
    eta_hat <- X_test %*% t(beta_all_post)
    mu_hat <- eta_hat + lambda_hat

    # convert back to gamma parameters
    b_post <- apply(
      mu_hat, 2,
      function(x, a){
        a / exp(x)
      },
      a = as.double(sample(a_post, size = n_test, replace = T))
    )

    # for each draw from the predictive dist of mu_hat, draw a value for y
    y_hat <- apply(
      b_post, MARGIN = 1,
      FUN = function(b, n, a){
        rgamma(n = n, shape = a, rate = b)
      },
      n = ncol(mu_hat),
      a = sample(a_post, size = ncol(mu_hat), replace = T)
    )

    # RMSE
    RMSE <- RMSE_bayes(y_test, y_hat)

    # How many non-zero params
    beta_post <- rstan::extract(mfit, pars = "beta")$beta
    qdf <- data.frame(
      low = apply(beta_post, 2, quantile, probs = (1 - q_cutoff)/2),
      high = apply(beta_post, 2, quantile, probs = 1 - (1 - q_cutoff)/2)
    )
    nzs_estim <- sum(qdf$low > 0 | qdf$high < 0)

    # KL-divergence between prior and posterior phi
    KLD <- kl_divergence(
      prior_samps = runif(length(phi_post)),
      post_samps = phi_post
    )

    return(list(
      y_test = y_test,
      y_hat = y_hat,
      RMSE = RMSE,
      nzs_true = length(sim_params$b),
      nzs_estim = nzs_estim,
      KLD_phi = KLD
    ))


  }


# # make cluster and run in parallel
#   nProc <- length(ns)
#   cl <- makeCluster(nProc)

#  # model fits
#  mfits <- foreach(n = ns, .packages = c('here', 'rstan')) %dopar%
#    latent_ar1_FHS(
#      n = n,
#      sim_params = sim_params,
#      mod_file = here("Stan/Pois_LatAR1_FHS.stan")
#    )

  # model fit
  mfit <- latent_ar1_FHS(
    n = as.numeric(args[1]),
    sim_params = sim_params,
    mod_file = here("Stan/Pois_LatAR1_FHS.stan")
  )

  # save the fit
  fp <- paste(
    here("Data/model_checks/"),
    "pois_latAr1_FHS_n",
    args[1],
    ".rds",
    sep = ""
  )
  saveRDS(mfit, file = fp)



