### ARMA(p,q) process ###

# libraries
library(rstan)

n <- 500

phi <- 0.5
theta <- 0.1

X <- cbind(
  rep(1, n),
  seq(-2, 2, length.out = n),
  matrix(
    rnorm(n * 49),
    nrow = n,
    ncol = 49
  )
)

beta <- c(0.5, 1, 2, rep(0, 48))

mu <- as.double(X %*% beta)
sigma_e <- 1

y <- stats::arima.sim(
  n = n,
  model = list(ar = c(0.5), ma = 0.1),
  sd = sigma_e,
  mean = mu
)

# testing stan model


arma_pq <- stan_model("Stan/ARMA-p-q_FHS.stan")
datlist <- list(
  N = n,
  P0 = 1,
  P = ncol(X) - 1,
  p = 1,
  q = 1,
  y = y,
  X_alpha = matrix(1, nrow = n, ncol = 1),
  X_beta = X[,-1],
  tau0 = 0.002,
  slab_scl = 2,
  slab_df = 10
)

mfit <- sampling(
  arma_pq,
  data = datlist,
  chains = 3, cores = 3
)




