# libraries
library(here)
library(rstan)
devtools::load_all()

post_mode <- function(x){
  id <- which.max(density(x)$y)
  return(density(x)$x[id])
}

# compile the stan model
  pois_garma11 <- stan_model(here("Stan/Pois_GLARMA-1-1_simple.stan"))

# function to simulate from a GLARMA(1,1) process and fit the model with
# 10 different fixed values of c and one random for augmenting the raw data

threshold_eval <- function(stan_mod){

  # simulating from a GLARMA(1,1) process

  # set parameters
    # time series length
    n <- 100

    # reg. coefficients
    beta <- c(-0.5, 1, 0.5)

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

    # threshold parameters
    thresh <- seq(0.05, 0.95, length.out = 10)

  # initiate the process
    y <- vector(mode = "double", length = n)
    y_star_true <- vector(mode = "double", length = n)
    mu <- vector(mode = "double", length = n)
    z <- vector(mode = "double", length = n)
    y[1] <- rpois(1, lambda = exp(X[1, ] %*% beta))
    y_star_true[1] <- max(y[1], runif(1))
    mu[1] <- y_star_true[1]
    z[1] <- 0

  # continue the process
    for(t in 2:n){

      # arma term
      z_t <- phi * (log(y_star_true[t-1]) - X[t-1, ] %*% beta) + theta * log(y_star_true[t-1]/mu[t-1])
      mu[t] <- exp(X[t, ] %*% beta + z_t)
      y[t] <- rpois(1, lambda = mu[t])
      y_star_true[t] <- max(y[t], runif(1))

    }

  # now define y_star according to different thresholds
    y_star <- vector(mode = "list", length = length(thresh) + 1)
    for(j in 1:length(thresh)){
      y_star[[j]] <- sapply(
        y,
        function(x, x2){
          max(x, x2)
        },
        x2 = thresh[j]
      )
    }
    # randomized threshold
    y_star[[length(y_star)]] <- sapply(
      1:n,
      function(i, x, x2){
        max(x[i], x2[i])
      },
      x = y,
      x2 = runif(n)
    )

    # compile data EXCEPT for y_star
    dat_remdat <- list(
      N = n,
      P = length(beta),
      y = y,
      X = X
    )

    # function to fit the model
    fit <- function(x, remain_dat, mod){

      # compile data to feed into stan
      datlist <- list(
        N = remain_dat$N,
        P = remain_dat$P,
        y = remain_dat$y,
        X = remain_dat$X,
        y_star = x
      )

      # fit mod
      mfit <- sampling(
        mod,
        data = datlist,
        cores = 3,
        chains = 3
      )

      return(mfit)

    }

    # make y_star a list

    # apply the fit to each y_star
    fits_all <- lapply(
      y_star,
      FUN = fit,
      remain_dat = dat_remdat,
      mod = stan_mod
    )

}
