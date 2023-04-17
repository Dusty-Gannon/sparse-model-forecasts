

#' Simulate seasonal AR time series
#'
#' @param n Length of the simulated time series
#' @param ar Vector of non-seasonal autoregressive parameters (phi).
#' @param sar Vector or scalar of seasonal autoregressive parameters.
#' @param S Length of the seasons (e.g., 12 time steps for monthly data)
#' @param sd Standard deviation of the random innovations
#' @param X n x K covariate matrix
#' @param beta K-vector of regression coefficients
#' @param burnin Length of the burnin period for the process
#'
#' @return List with y as the response and the parameters used to generate the series
#'
sim_sAR_p <- function(n, ar, sar, S, sd = 1, X = NULL, beta = NULL, burnin = NULL){

  # some useful variables
  p <- length(ar); P <- length(sar)
  ar_poly <- c(1, -ar)
  sar_poly <- vector(mode = "double", length = 1 + S * P)
  sar_poly[1] <- 1
  sar_poly[S * 1:P + 1] <- -sar

  # do the polynomial multiplication
  m <- matrix(nrow = length(ar_poly), ncol = length(sar_poly))
  for(i in 1:length(ar_poly)){
    for(j in 1:length(sar_poly)){
      m[i,j] <- ar_poly[i] * sar_poly[j]
    }
  }

  # pad the end
  m <- cbind(m, matrix(data = 0, nrow = nrow(m), ncol = nrow(m) - 1))

  # combine like terms
  char_poly <- vector(mode = "double")
  for(j in 1:ncol(m)){
    i <- 1
    j2 <- j
    temp_sum <- 0
    while(j2 > 0 & i <= nrow(m)){
      temp_sum <- temp_sum + m[i, j2]
      i <- i + 1
      j2 <- j2 - 1
    }
    char_poly <- c(char_poly, temp_sum)
  }

  # define phi
  phi <- -char_poly[-1]
  p2 <- length(phi)

  if(is.null(burnin)){
    burnin <- length(phi) * 5
  }
  # complete centered process
  if(is.null(X)){
    y <- rep(0, burnin + n)
    for(t in (p2 + 1):length(y)){
      y[t] <- y[(t - 1):(t - p2)] %*% phi + rnorm(1, sd = sd)
    }
  }

  # if covariates / mean are included
  if(!is.null(X)){

    # create a burnin matrix by sampling from the provided matrix
    X_burn <- apply(X, 2, FUN = sample, size = burnin, replace = T)
    X2 <- rbind(X_burn, X)
    y <- vector(mode = "double", length = burnin + n)
    y[1:p2] <- X_burn[1:p2, ] %*% beta
    for(t in (p2 + 1):length(y)){
      y[t] <- X2 %*% beta + y[(t - 1):(t - p2)] %*% phi + rnorm(1, sd = sd)
    }
  }

  return(list(
    y = y[(burnin + 1):length(y)],
    phi = phi,
    X = X,
    beta = beta,
    sd = sd
  ))

}
