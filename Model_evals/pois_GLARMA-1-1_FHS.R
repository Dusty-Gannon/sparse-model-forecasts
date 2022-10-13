# libraries
library(rstan)

# simulating from a GLARMA(1,1) process

# set parameters
  # time series length
  n <- 50

  # intercept
  alpha <- 1

  # nz slopes
  beta_nz <- c(1, -1, 0.5)

  # combine
  beta <- c(
    alpha,
    beta_nz,
    rep(0, 46)
  )

  # AR and MA parameters
  phi <- 0.6
  theta <- 0.5

# model matrix
  X <- matrix(
    rnorm(n = n * (length(beta) - 1)),
    nrow = n
  )
  X <- cbind(
    rep(1, n),
    X
  )

# initiate the process
  y <- vector(mode = "double", length = n)
  mu <- vector(mode = "double", length = n)
  z <- vector(mode = "double", length = n)
  y[1] <- rpois(1, lambda = exp(X[, 1] %*% beta))
  mu[1] <- y[1]
  z[1] <- 0

# continue the process

  for(t in 2:n){
    # Pearson residual
    e_tm1 <- (y[t-1] - mu[t-1])/sqrt(mu[t-1])

    # arma term
    z_t <- phi * (z[t-1] + e_tm1) + theta * e_tm1
    lmu_t <- X[, t] %*% beta + z_t

    # fill in series
    mu[t] <- exp(lmu_t)

    y[t] <- rpois(1, lambda = mu[t])

  }




