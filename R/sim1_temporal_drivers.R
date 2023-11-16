


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
#' @param probWeakCorr probability of a weak driver being correlated to a strong driver
#' @param probStrongCorr probability of a strong driver being correlated to a weak driver
#' @param corrLevel Correlation level desired between strong and weak drivers, if they are correlated
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
    trend_fraction = 0.5, prob_cycle = 0.5, sigma = 0.5, probWeakCorr=0.2, probStrongCorr=0.5, corrLevel=0.7
  ){

  # get total number of samples
  N <- n * freq

  # total number of drivers




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

  # ar coefficients (excluded from these time series)
  # ars <- runif(K, max = 0.8)

  # variances
  sds <- runif(K, min = 0.3, max = 2)

  # give some covariates yearly cycles
  S <- sample(
    1:round(freq / 2),
    size = K,
    replace = T
  ) * rbinom(K, size = 1, prob = prob_cycle)


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

  # substitute in the correlated covariates
  # goal: we want to make some weakly explanatory drivers correlate with some strong explanatory drivers
  # for now, weak drivers do not correlate with each other, and strong drivers do not correlate with each other
  correlatedVariable = function(x, r){
    r2 = r^2
    ve = 1-r2
    SD = sqrt(ve)
    e  = rnorm(length(x), mean=0, sd=SD)
    y  = r*x + e
    return(y)
  }

  if(probWeakCorr>0 & probStrongCorr>0){
    # remove the intercept column
    beta0=beta[-1]

    # select the strong drivers to correlate to
    b1=which(abs(beta0)>0.5)
    # randomly, choose if a strong driver has a correlated weak driver
    strong1=runif(length(b1),0,1)<probStrongCorr
    # gives the strong drivers selected to have correlated weak drivers
    strong2=b1[strong1]


    # select the weak drivers to correlate to
    b2=which(abs(beta0)<=0.5)
    # remove the first beta (zero)

    # randomly, choose if a weak driver has a correlated strong driver
    weak1=runif(length(b2),0,1)<probWeakCorr
    # gives the weak drivers selected to have correlated strong drivers
    weak2=b2[weak1]

    # Randomly assign weak and strong drivers together
    if(length(strong2)>0&length(weak2)>0){
      # assign each strong driver at least one weak driver
      strongAssign=sample(weak2, length(strong2), replace=F)
      # calculate the new correlated variables
      newCorrVars1=lapply(xvars[strong2],correlatedVariable,corrLevel)
      # plug the new correlated variables in for the appropriate weak drivers
      xvars[strongAssign]<-newCorrVars1

      if(length(weak2)>length(strong2)){
        # assign the rest of the weak ones to chosen strong ones
        weakAssign=sample(strong2,length(weak2)-length(strong2),replace=T)
        # calculate the new correlated variables
        newCorrVars2=lapply(xvars[weakAssign],correlatedVariable,corrLevel)
        # plug the new correlated variables in for the appropriate weak drivers
        weak3=weak2[-which(weak2%in%strongAssign)]
        xvars[weak3]=newCorrVars2
      }

    }

  }

  # convert to a matrix
  X <- cbind(
    rep(1, N),
    Reduce(cbind, xvars)
  )


  # construct the response
  y <- as.double(X %*% beta) + rnorm(n, sd = sigma)

  # Store everything together in a list
  ts <- list(y = y, X = X, beta = beta)
  return(ts)
}


