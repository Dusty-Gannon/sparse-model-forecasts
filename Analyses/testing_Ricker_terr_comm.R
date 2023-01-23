###################################################
# Fitting Ricker models to simulated community data
###################################################


# libraries
  library(rstan)
  library(here)
  library(parallel)
  devtools::load_all()


##### Generate some parameter values that will be consistent through the following tests #####

  steps <- 100                                       # number of time steps
  nsp <- 10                                          # number of species in the community
  lambdas <- runif(nsp, min = 1.2, max = 1.8)          # intrinsic growth rates
  alphas <- runif(nsp, min = 0.005, max = 0.01)        # intraspecific competition
  sigma <- 0.1                                      # initial abundances

  A_mat <- comp_matrix2(
    n_sp = nsp,
    rho = 0,
    num_ngs = 2,
    alpha = alphas,
    ng_range = c(0.1, 0.5)
  )


##### Simulate with randomly selected heterospecific abundances #####

  # initialize response
  N_foc <- vector(mode = "double", length = steps)
  N_foc[1] <- N0[1]

  # randomly generate heterospecific abundances for each time step
  N_het <- t(sapply(
    1:steps,
    FUN = function(X, lambdas, alphas){
      n <- length(lambdas)
      rnorm(n, log(lambdas)/alphas/nsp)
    },
    lambdas = lambdas[-1],
    alphas = alphas[-1]
  ))

  # create response vector
  for(t in 2:steps){
    N_foc[t] <- N_foc[t - 1] * lambda[1] * exp(-alphas[1] * N_foc[t - 1] - N_het[t - 1, ] %*% A_mat[1, -1]) +
      rnorm(1) * sigma / sqrt(N_foc[t - 1])
  }

  # compile data list
  tau_0 <- tau0(
    y = log(N_foc[2:steps]/N_foc[1:(steps - 1)]),
    m0 = 5,
    M = ncol(N_het),
    N = steps - 1,
    fam = "gaussian"
  )

  datlist_randhet <- list(
    N = steps,
    P = ncol(N_het),
    y = N_foc,
    X_beta = N_het,
    error_scl = 0.5,
    tau0 = tau_0,
    slab_scl = 1,
    slab_df = 10
  )

  # fit the model
  growth_mod <- stan_model(here("Stan/pop_growth_rate_FHS.stan"))
  mfit <- sampling(
    growth_mod,
    data = datlist_randhet,
    cores = 3, chains = 3,
    control = list(adapt_delta = 0.99, max_treedepth = 15)
  )


##### Simulate full community #####

  N <- matrix(nrow = nsp, ncol = steps)
  N[, 1] <- rnorm(nsp, mean = 10)
  Sigma <- runif(nsp, min = 0.05, max = 0.2)

  for(t in 2:steps){

    N[, t] <- N[, t - 1] * lambdas * exp(- A_mat %*% N[, t - 1] + rnorm(nsp) * Sigma / N[, t - 1])

  }

  N_het_full <- t(N[-1, ])
  N_foc_full <- as.double(N[1, ])
  tau_0_full <- tau0(
    y = log(N_foc_full[2:steps] / N_foc_full[1:(steps - 1)]),
    m0 = 2,
    M = nsp - 1,
    N = steps - 1,
    fam = "gaussian"
  )

  datlist_full <- list(
    N = steps,
    P = nsp - 1,
    y = N_foc_full,
    X_beta = N_het_full,
    error_scl = 0.5,
    tau0 = tau_0_full,
    slab_scl = 0.5,
    slab_df = 4
  )

  mfit_full <- sampling(
    growth_mod,
    data = datlist_full,
    cores = 3, chains = 3,
    control = list(adapt_delta = 0.99, max_treedepth = 15)
  )

