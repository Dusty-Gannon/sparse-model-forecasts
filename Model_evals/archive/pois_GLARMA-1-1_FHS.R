# libraries
library(here)
library(rstan)
devtools::load_all()

post_mode <- function(x){
  id <- which.max(density(x)$y)
  return(density(x)$x[id])
}

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
  y_star <- vector(mode = "double", length = n)
  mu <- vector(mode = "double", length = n)
  z <- vector(mode = "double", length = n)
  y[1] <- rpois(1, lambda = exp(X[1, ] %*% beta))
  y_star[1] <- max(y[1], 0.1)
  mu[1] <- y_star[1]
  z[1] <- 0

# continue the process

  for(t in 2:n){

    # arma term
    z_t <- phi * (log(y_star[t-1]) - X[t-1, ] %*% beta) + theta * log(y_star[t-1]/mu[t-1])
    mu[t] <- exp(X[t, ] %*% beta + z_t)
    y[t] <- rpois(1, lambda = mu[t])
    y_star[t] <- max(y[t], 0.1)

  }


# compute tau0
  tau_0 <- tau0(
    m0 = 5,
    M = length(beta) - length(alpha),
    N = n,
    sigma = mean(y)^(-1)
  )

# compile data to feed into stan
  datlist <- list(
    N = n,
    P0 = length(alpha),
    P = length(beta) - length(alpha),
    y = y,
    y_star = y_star,
    X_alpha = as.matrix(X[, 1]),
    X_beta = X[, -1],
    tau0 = tau_0,
    slab_scl = 1,
    slab_df = 10
  )



# fitting with Stan
  pois_garma11 <- stan_model(here("Stan/Pois_GLARMA-1-1_FHS.stan"))
  mfit <- sampling(
    pois_garma11,
    data = datlist,
    cores = 3,
    chains = 3,
    control = list(adapt_delta = 0.95, max_treedepth = 12)
  )








