######
# This script simulates conditional Poisson time series with
# a latent AR(1) process to induce autocorrelation and 3 non-zero
# regression coefficients, then fits them using Finnish-Horseshoe
# priors for the regression coefficients to induce sparsity
######

# libraries
  library(rstan)
  # library(parallel)
  # library(snow)
  # library(Rmpi)
  devtools::load_all()

# number of datasets to simulate
  sims <- 100

# define parameters
  sim_params <- list(
    n = 100,           # length of the time series
    b = c(1, -1, 0.5), # non-zero effects
    P = 50,            # total number of candidate variables (minus the intercept)
    sigma_e = 0.5      # sd of random innovations
  )


# function to simulate the data, fit the model, and report RMSE stats
  latent_ar1_FHS <- function(sims, sim_params, mod_file, q_cutoff = 0.9, n_test = 30){

    # draws for randomized parameters
    phi <- runif(1, min = 0.2, max = 1)

    # useful variables
    n_train <- sim_params$n - n_test

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

    # simulate response vector
    y <- with(sim_params, rpois(n, mu))

    # split data into training and testing data
    y_train <- y[1:n_test]
    y_test <- y[(n_test + 1):sim_params$n]
    X_train <- X[1:n_test, ]
    X_test <- X[(n_test + 1):sim_params$n, ]

    # estimate for tau0
    tau0 <- tau0(
      m0 = 5,
      M = sim_params$P,
      N = n_test,
      sigma = mean(y_train)^(-1)
    )

    # data list for model fit
    datlist <- with(sim_params,{
      list(
        N = n_test,
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
    Pois_latar1_FHS <- stan_model("Stan/Pois_LatAR1_FHS.stan")

    # fit the model
    mfit <- sampling(
      Pois_latar1_FHS,
      data = datlist,
      chains = 3,
      iter = 2000,
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

    # for each draw from the predictive dist of mu_hat, draw a value for y
    y_hat <- apply(
      mu_hat, MARGIN = 1,
      FUN = function(mu, n){
        rpois(n = n, lambda = exp(mu))
      },
      n = ncol(mu_hat)
    )

    # RMSE
    RMSE <- RMSE_bayes(y_test, y_hat)

    # How many non-zero params
    beta_post <- rstan::extract(mfit, pars = "beta")$beta
    qdf <- data.frame(
      low = apply(beta_post, 2, quantile, probs = (1 - q_cutoff)/2),
      high = apply(beta_post, 2, quantile, probs = 1 - (1 - q_cutoff)/2)
    )
    nzs <- sum(qdf$low > 0 | qdf$high < 0)



  }














