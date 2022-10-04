#' Construct AR1 correlation matrix
#'
#' @param n Number of dimensions for the correlation matrix
#' @param rho Correlation between observations y(t) and y(t-1)
#'
#' @return n x n matrix
#' @export
#'
#' @examples
#'
ar1_cormat <- function(n, rho){
  # vector of ones
  ones <- matrix(rep(1, n), ncol = 1)

  # vector of values 0, 1, ..., n-1
  nvec <- matrix(1:n - 1, ncol = 1)

  # construct matrix of powers for rho
  pow <- abs(nvec %*% t(ones) - ones %*% t(nvec))

  return(rho^pow)
}







#' Draw from a Dirichlet distribution
#'
#' @param n Number of random vectors to draw
#' @param theta Parameter vector
#'
#' @return Either a matrix with n rows and as many columns as elements in \code{theta}
#' @export
#'
#' @examples
rdirch <- function(n, theta){
  K <- length(theta)
  g <- matrix(nrow = n, ncol = K)
  for(i in 1:n){
    for(k in 1:K){
      g[i, k] <- rgamma(1, shape = theta[k])
    }
  }

  samps <- matrix(nrow = n, ncol = K)
  for(i in 1:n){
    samps[i, ] <- g[i, ]/sum(g[i, ])
  }
  if(n == 1){
    return(as.double(samps))
  } else{
    return(samps)
  }
}








#' Calculate tau0 based on Piironen and Vehtari (2017)
#'
#' @param m0 Prior guess for number of non-zero effects
#' @param M Number of coefficients getting shrinkage priors
#' @param N Number of observations
#' @param sigma Approximate variance (see Piironen and Vehtari, 2017)
#'
#' @return real
#' @export
#'
#' @examples
#' # for a Poisson regression with 200 observations and 50 candidate covariates
#' tau0(m0 = 2, M = 50, N = 200, sigma = mean(y)^(-1))
#'
tau0 <- function(m0, M, N, sigma){
  m0/(M-m0) * sigma/sqrt(N)
}



