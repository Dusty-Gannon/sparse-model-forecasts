#' Run a basic time series with no autoregressive components and all seasonal
#' fluctuations contained within the covariates.
#' Aguments:
#' K = number of covariates
#' n = number of monthly data points
#' trend_fraction = fraction of covariates with a trend added
#' NumLargeEffect = Number of large effect covariates
#' ProbCycle = Probability of covariates experiencing seasonal cycles
#' @return
#' A list containing a y vector for the time series and the model matrix of
#' covariates used in the simulation with the beta vector of 'true' regression
#' coefficients
#'
basic_timeseries <- function(K, n, trend_fraction, NumLargeEffect, ProbCycle){
  # add trends to some fraction
  trends <- c(
    rnorm(round(K * trend_fraction), mean = 0, sd = 0.2),
    rep(0, K - round(K * trend_fraction))
  )

  # give all of them random means
  means <- rnorm(K)

  # construct the model matrix for all of them
  X_covs <- cbind(
    rep(1, n),
    (1:n)/12
  )

  # ar coefficients
  ars <- runif(K, max = 0.8)

  # variances
  sds <- runif(K, min = 0.3, max = 2)

  # give some covariates yearly cycles
  S <- sample(
    c(1, 12),
    size = K,
    replace = T,
    prob = c(1-ProbCycle, ProbCycle)
  )

  # construct the variables
  xvars <- map(
    1:K,
    ~ {
      sim_sARp(
        n = n,
        ar = ars[.x],
        S = S[.x],
        sd = sds[.x],
        X = X_covs,
        beta = c(means[.x], trends[.x]),
        burnin = 100
      )$y
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
  y <- as.double(X %*% beta) + rnorm(n, sd = 0.5)

  # Store everything together in a list
  TimeSeries <- list(y = y, X = X, beta = beta)
  return(TimeSeries)
}
