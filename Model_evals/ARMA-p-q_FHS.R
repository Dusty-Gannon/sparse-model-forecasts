### AR(p) process ###
# libraries
library(rstan)
library(here)

# length of time series
n <- 200

# time series parameters
phi <- c(0.5, -0.3, 0.3)

# model matrix with two non-zero coefficients
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

# simulate AR(1) process

y <- arima.sim(
  n = n,
  model = list(ar = phi),
  mean = mu,
  sd = sigma_e
)

# compile stan model
arp <- stan_model(here("Stan/AR-p_FHS.stan"))

# compile data
datlist <- list(
  N = n,
  P0 = 1,
  P = ncol(X) - 1,
  p = 3,
  y = y,
  X_alpha = matrix(1, nrow = n, ncol = 1),
  X_beta = X[,-1],
  tau0 = 0.002,
  slab_scl = 2,
  slab_df = 10
)

mfit_arp <- sampling(
  arp,
  data = datlist,
  chains = 3, cores = 3
)



### ARMA(p,q) process ###


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




