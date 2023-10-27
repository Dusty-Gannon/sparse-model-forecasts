


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
#' @param num_strong Number of strong drivers
#' @param prob_cycle Probability that a driver experiences a seasonal cycle.
#' @param sigma Standard deviation of the noise component.
#' @param correlated Whether or not to generate correlated covariates
#' @param rateCorr magnitude of correlation between covariates
#'
#' @return List with the response, \code{y}, the model matrix, \code{X},
#' and the regression coefficients used to construct the response, \code{beta}.
#'
#' @examples
#' # simulating a 3-year time series with data from each day,
#' # constructed with 10 drivers, 3 strong and 7 weak.
#'
#' basic_timeseries(
#'   K = 10,
#'   num_strong = 3,
#'   n = 3,
#'   freq = 365
#' )
#'
#'
basic_timeseries <- function(
    K, num_strong, n, freq,
    trend_fraction = 0.5, prob_cycle = 0.5, sigma = 0.5, correlated=F, rateCorr=2
  ){

  # get total number of samples
  N <- n * freq

  # add trends to some fraction
  trends <- c(
    rnorm(round(K * trend_fraction), mean = 0, sd = 2 / n),
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
  ) * rbinom(K, size = 1, prob = prob_cycle)

  # if correlated covariates are desired, build them here:
  if(correlated){
    # To create a covariance matrix, step one is to create an orthogonal matrix,
    # which can be done using QR decomposition of an arbitrary matrix
    Q <- qr.Q(qr(matrix(rnorm(K^2), nrow = K, ncol = K)))

    # step two is to generate a diagonal matrix with the standard deviations
    D <- diag(x = rgamma(K, shape = 2, rate = rateCorr)) # making the rate parameter smaller increases the possible cov values

    # now use matrix multiplication to generate Sigma (covariance matrix)
    Sigma <- t(Q) %*% D %*% Q

    # to test that this is positive definite, check that all eigenvalues are positive
    # sum(eigen(Sigma)$values < 0)

    # create un-correlated matrix by using random draws
    u=matrix(rnorm(N * K),
             nrow = N,
             ncol = K)
    # then Cholesky decompose Sigma and multiply by the uncorrelated matrix
    X_corr <- u %*% chol(Sigma)

    # construct the variables
    xvars <- purrr::map(
      1:K,
      ~ {
        if(S[.x] > 0){
          runif(1, max = 0.5) * cos(pi * 1:N / S[.x]) +
            runif(1, max = 0.5) * sin(pi * 1:N / S[.x]) +
            means[.x] + trends[.x] * 1:N / freq +
            X_corr[,.x]
        } else {
          mean = means[.x] + trends[.x] * 1/N * freq + X_corr[,.x]

        }
      }
    )

  } else {
    # construct the variables without correlation
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
  }

  # convert to a matrix
  X <- cbind(
    rep(1, N),
    Reduce(cbind, xvars)
  )

  # construct beta
  k_prime <- K - num_strong # number of small effect variables

  beta <- c(
    runif(num_strong, min = 0.5, max = 1) * sample(c(-1,1), num_strong, replace = T),
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
  ts <- list(y = y, X = X, beta = beta)
  return(ts)
}


