



#' Bayesian Root-mean squared error
#'
#' @param obs Vector of observed values
#' @param ppreds Matrix of posterior predictions,
#' one column per observation and as many rows as draws from the posterior predictive distribution
#'
#' @return Vector of draws from the posterior distribution of the RMSE
#'
#' @examples
#'
RMSE_bayes <- function(obs, ppreds){
  apply(
    ppreds, 1,
    function(x, obs){
      sqrt(mean((x - obs)^2))
    },
    obs = obs
  )
}


#' Residuals plots
#'
#' Plot a vector of residuals against each of the explanatory variables in X
#'
#' @param r Residuals
#' @param X \eqn{n \times P} matrix of \eqn{P} explanatory variables to plot against
#' @param cnames If \code{TRUE}, \code{X} has column names.
#'
#' @return Faceted plot of residuals against each explanatory variable
#' @export
#'
#' @examples
resid_v_x <- function(r, X, cnames = T){

  # are there column names
  if(isTRUE(cnames)){
    # construct long-dataframe
    df <- data.frame(
      r = rep(r, ncol(X)),
      x = as.vector(X),
      variable = rep(colnames(X), each = nrow(X))
    )
  } else{
    # construct long-dataframe
    df <- data.frame(
      r = rep(r, ncol(X)),
      x = as.vector(X)
    )

    # makeshift column names
    index <- 1:ncol(X)
    vname <- paste("V", as.character(index), sep = "")
    df$variable <- rep(vname, each = nrow(X))
  }

  # plot
  ggplot(data = df, aes(x = x, y = r))+
    geom_point()+
    geom_hline(yintercept = 0)+
    facet_wrap(~ variable, ncol = 3, scales = "free")

}


















