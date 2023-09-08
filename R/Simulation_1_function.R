


#' Simulate a time series using temporally autocorrelated covariates
#'
#' This function simulates a time series by creating autoregressive and
#' seasonal drivers and constructing a response based on a linear combination
#' of those drivers.
#'
#' @param K Number of drivers.
#' @param n Length of time series
#' @param freq Frequency of sampling per time step (e.g., \code{n = 5} and
#' \code{freq = 365}) would create 5 years of daily data.
#' @param trend_fraction Fraction of drivers that have a trend.
#' @param NumLargeEffect Number of strong drivers
#' @param ProbCycle Probability that a driver experiences a seasonal cycle.
#' @param sigma Standard deviation of the noise component.
#'
#' @return List with the response, \code{y}, the model matrix, \code{X},
#' and the regression coefficients used to construct the response, \code{beta}.
#'
basic_timeseries <- function(K, n, freq, trend_fraction, NumLargeEffect, ProbCycle, sigma = 0.5){

  # get total number of samples
  N <- n * freq

  # add trends to some fraction
  trends <- c(
    rnorm(round(K * trend_fraction), mean = 0, sd = 0.2),
    rep(0, K - round(K * trend_fraction))
  )

  # give all of them random means
  means <- rnorm(K)

  # construct the model matrix for all of them
  X_covs <- cbind(
    rep(1, N),
    (1:N) / freq
  )

  # ar coefficients
  ars <- runif(K, max = 0.8)

  # variances
  sds <- runif(K, min = 0.3, max = 2)

  # give some covariates yearly cycles
  S <- sample(
    1:round(freq / 2),
    size = K,
    replace = T
  ) * rbinom(K, size = 1, prob = ProbCycle)

  # construct the variables
  xvars <- purrr::map(
    1:K,
    ~ {
      if(S[.x] > 0){
        runif(1, max = 0.5) * cos(pi * 1:N / S[.x]) +
        runif(1, max = 0.5) * sin(pi * 1:N / S[.x]) +
          means[.x] + trends[.x] * 1:N / freq +
          rnorm(N, sd = sds[.x])
      } else {
        rnorm(
          N,
          mean = means[.x] + trends[.x] * 1/N * freq,
          sd = sds[.x]
        )
      }
    }
  )

  # convert to a matrix
  X <- cbind(
    rep(1,n),
    Reduce(cbind, xvars)
  )

  # construct beta
  k_prime <- K - NumLargeEffect # number of small effect variables

  beta <- c(
    runif(NumLargeEffect, min = 0.5, max = 1) * sample(c(-1,1), NumLargeEffect, replace = T),
    rnorm(k_prime, sd = 0.05)
  )

  # distribute these randomly
  beta <- c(
    0,
    beta[sample(1:K, size = K)]
  )

  # construct the response
  y <- as.double(X %*% beta) + rnorm(n, sd = sigma)

  # Store everything together in a list
  TimeSeries <- list(y = y, X = X, beta = beta)
  return(TimeSeries)
}
