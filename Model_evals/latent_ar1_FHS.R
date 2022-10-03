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
    n = 30,            # length of the time series
    b = c(1, -1, 0.5), # non-zero effects
    P = 50,            # total number of candidate variables
    sigma_e = 0.5      # sd of random innovations
  )


# function to simulate the data, fit the model, and report RMSE stats
  latent_ar1_FHS <- function(sims, sim_params, mod_file, q_cutoff = 0.8){

    # draws for randomized parameters
    phi <- runif(1)

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

    # estimate for tau0
    tau0 <- tau0(
      m0 = 5,
      M = sim_params$P,
      N = sim_params$n,
      sigma = mean(y)^(-1)
    )

    # data list for model fit
    datlist <- with(sim_params,{
      list(
        N = n,
        P = 1,
        S = ncol(X) - 1,
        y = y,
        X_alpha = as.matrix(X[,1]),
        X_beta = X[,-1],
        tau0 = tau0,
        slab_scl = 1,
        slab_df = 6
      )
    })

    # compile stan model
    latar1_sparse <- stan_model("Stan/Pois_LatAR1_sparse_reg.stan")

    # fit the model
    mfit <- sampling(
      latar1_sparse,
      data = datlist,
      chains = 3,
      iter = 2000,
      control = list(adapt_delta = 0.99)
    )

    # compute RMSE for parameters of interest
    # first extract linear predictor and exponentiate
    mu_post <- exp(rstan::extract(mfit, pars = "eta")$eta)

    # create posterior predictive draws
    y_hat <- matrix(nrow = nrow(mu_post), ncol = ncol(mu_post))
    for(i in 1:nrow(y_hat)){
      y_hat[i, ] <- rpois(ncol(y_hat), lambda = mu_post[i, ])
    }

    # RMSE
    RMSE <- RMSE_bayes(y, y_hat)

    # How many non-zero params
    beta_post <- rstan::extract(mfit, pars = "beta")$beta
    qdf <- data.frame(
      low = apply(beta_post, 2, quantile, probs = (1 - q_cutoff)/2),
      high = apply(beta_post, 2, quantile, probs = 1 - (1 - q_cutoff)/2)
    )
    nzs <- sum(qdf$low > 0 | qdf$high < 0)



  }














