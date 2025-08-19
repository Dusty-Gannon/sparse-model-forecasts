
#' Calculate tau0 based on Piironen and Vehtari (2017)
#'
#' @param y Vector of responses
#' @param m0 Prior guess for number of non-zero effects
#' @param M Number of coefficients getting shrinkage priors
#' @param N Number of observations
#' @param fam Exponential family distribution for computing pseudo-variance (see Piironen and Vehtari, 2017).
#' The options (for now) are \code{fam = c("poisson", "gamma", "gaussian")}
#'
#' @return real
#' @export
#'
#' @examples
#' # for a Poisson regression with 200 observations and 50 candidate covariates
#' tau0(y = y, m0 = 2, M = 50, N = 200, fam = "poisson")
#'
tau0 <- function(y, m0, M, N, fam){
  if(sum(fam %in% c("poisson", "gamma", "gaussian")) == 0){
    stop("This family is not supported (yet). To add it, see Table 1 in Piironen and Vehtari (2017)")
  }
  if(fam == "gaussian"){
    sigma <- sd(y)
  }
  if(fam == "poisson"){
    sigma <- mean(y)^{-0.5}
  }
  if(fam == "gamma"){
    sigma = sd(y)
  }

  m0/(M-m0) * sigma/sqrt(N)
}



