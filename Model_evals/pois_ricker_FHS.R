# libraries
library(here)
library(rstan)
devtools::load_all()

post_mode <- function(x){
  id <- which.max(density(x)$y)
  return(density(x)$x[id])
}

# simulating from a Ricker model

# set parameters
  # time series length
  n <- 100

  # intercept
  alpha <- 1

  # nz slopes
  beta_nz <- -0.05


# initiate the process
  y <- vector(mode = "double", length = n)
  y[1] <- 10

# continue the process

  for(t in 2:n){

    # arma term
    y[t] = rpois(
      1, lambda = y[t - 1] * exp(alpha + beta_nz * y[t - 1])
    )

  }

# model matrix
  X <- matrix(
    rnorm(n = (n - 1) * 20),
    nrow = n - 1
  )
  X <- cbind(
    y[1:(n - 1)],
    X
  )


# compute tau0
  tau_0 <- tau0(
    m0 = 5,
    M = ncol(X),
    N = n - 1,
    sigma = mean(y[2:n])^(-1)
  )

# compile data to feed into stan
  datlist <- list(
    N = n - 1,
    P0 = length(alpha),
    P = ncol(X),
    y = y[2:n],
    X_alpha = matrix(rep(1, n - 1), ncol = 1),
    X_beta = X,
    tau0 = tau_0,
    slab_scl = 1,
    slab_df = 10
  )



# fitting with Stan
  pois_ricker <- stan_model(here("Stan/Pois_ricker_FHS.stan"))
  mfit <- sampling(
    pois_ricker,
    data = datlist,
    cores = 3,
    chains = 3,
    control = list(adapt_delta = 0.95, max_treedepth = 12)
  )








