
#' Continue a set of Fourier basis functions into the future
#'
#' This function takes a matrix of Fourier basis functions and continues them
#' into the future by \code{h} time steps.
#'
#' @param X Matrix of Fourier basis functions with \eqn{\sin(\cdot)} in the odd
#' columns and \eqn{\cos(\cdot)} in the even columns. This is the standard formatting
#' coming from \code{forecast::fourier()}. The matrix \emphasis{should not} include
#' the column of ones for the intercept.
#' @param h Number of steps ahead to forecast.
#' @param n Number of time steps in the original series that was used to create the
#' Fourier basis.
#' @param append Logical. If \code{TRUE}, the function will append the extend the
#' Fourier basis functions, adding them to \code{X}. If \code{FALSE}, the function
#' will return the Fourier basis functions from the last time step in \code{X} to
#' \code{h} time steps into the future.
#'
#' @returns A matrix of Fourier basis functions.
#' @export
#'
continue_fourier <- function(X, h, n, append = FALSE){

  # number of Fourier terms
  p <- ncol(X) / 2
  X_new <- matrix(0, nrow = h, ncol = ncol(X))
  tsteps <- (nrow(X) + 1):(nrow(X) + h)

  for(j in 1:p){
    X_new[, 2 * j - 1] <- sin(2 * pi * j * tsteps / n)
    X_new[, 2 * j] <- cos(2 * pi * j * tsteps / n)
  }

  if(append){
    return(rbind(X, X_new))
  } else{
    colnames(X_new) <- colnames(X)
    return(X_new)
  }

}
