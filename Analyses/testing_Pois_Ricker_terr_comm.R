###################################################
# Fitting Ricker models to simulated community data
###################################################


# libraries
  library(rstan)
  library(here)
  library(parallel)
  devtools::load_all()


##### Generate some parameter values that will be consistent through the following tests #####

  set.seed(9999)

  steps <- 100                                      # number of time steps
  nsp <- 50                                         # number of species in the community
  lambdas <- runif(nsp, min = 1, max = 2)           # intrinsic growth rates
  alphas <- runif(nsp, min = 0.005, max = 0.01)        # intraspecific competition
  N0 <- rpois(nsp, 20)                              # initial abundances

  # competition matrix
  A_mat <- comp_matrix2(
    n_sp = nsp,
    alpha = alpha,
    rho = 0,
    num_ngs = 3,
    ng_range = c(0.6, 0.6)
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
      rpois(n, log(lambdas)/alphas)
    },
    lambdas = lambdas[-1],
    alphas = alphas[-1]
  ))

  # create response vector
  for(t in 2:steps){
    mu_t <- N_foc[t - 1] * lambdas[1] * exp(-(alphas[1] * N_foc[t - 1] + N_het[t - 1, ] %*% A_mat[1, -1]))
    N_foc[t] <- rpois(1, lambda = mu_t)
  }

  # compile data list
  tau_0 <- tau0(
    N_foc,
    m0 = 5,
    M = ncol(N_het),
    N = steps,
    fam = "poisson"
  )

  datlist_randhet <- list(
    N = steps,
    P = ncol(N_het),
    y = N_foc,
    X_beta = N_het,
    lambda = as.numeric(lambdas[1]),
    alpha = as.numeric(alphas[1]),
    tau0 = tau_0,
    slab_scl = 1,
    slab_df = 10
  )

  # fit the model
  ricker_mod <- stan_model(here("Stan/Pois_ricker_FHS.stan"))
  mfit_randhet <- sampling(
    ricker_mod,
    data = datlist_randhet,
    cores = 3, chains = 3,
    control = list(adapt_delta = 0.99, max_treedepth = 12)
  )




