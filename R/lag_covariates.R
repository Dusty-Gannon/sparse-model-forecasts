
#' Create dataset with lagged variables
#'
#' @param data A dataframe containing the variables
#' to be lagged. Other variables may be included in the dataframe.
#' @param names Names of the columns to lag.
#' @param lags A vector or scalar of the number of lags for each variable.
#' @param time_col An optional character string giving the column name for
#' the time variable. Default is \code{NULL}.
#'
#' @returns A data frame with \eqn{n - p^*} rows and \eqn{K + \sum_k p_k}
#' columns, where \eqn{n} is the number of rows in the original data,
#' \eqn{p^*} is the maximum lag number, \eqn{K} is the number of columns in
#' the original data set, and \eqn{p_k} is the number of lags to include for
#' column \eqn{k}.
#' @export
#'
lag_covariates <- function(data, names, lags, time_col = NULL){

  # get lagged column indexes
  lag_cols <- which(names(data) %in% names)

  # get column indexes to leave alone
  nolag_cols <- which(!(names(data) %in% names))

  # get dimensions for new dataset
  maxlag <- max(lags)
  n <- nrow(data) - maxlag    # num rows
  K <- ncol(data) + sum(lags) # num columns
  if(!is.null(time_col)){
    newdat <- data.frame(data[, time_col])
    names(newdat) <- time_col
    newdat <- newdat[(maxlag + 1):nrow(data), ]
    nolag_cols <- nolag_cols[-which(nolag_cols == which(names(data) == time_col))]
  } else{
    newdat <- data.frame(row.names = (maxlag + 1):nrow(data))
  }

  # build new data set inside loop
  for(j in 1:length(lag_cols)){

    m_j <- matrix(nrow = n, ncol = lags[j] + 1)
    m_j[, 1] <- data[(maxlag + 1):(nrow(data)), lag_cols[j]]
    for(k in 1:lags[j]){
      m_j[, k + 1] <- data[(maxlag + 1 - k):(nrow(data) - k), lag_cols[j]]
    }
    colnames(m_j) <- c(
      names[j],
      paste0(names[j], "_", "l", 1:lags[j])
    )

    # now bind with newdata
    newdat <- cbind(newdat, as.data.frame(m_j))
  }

  # now add the remaining columns
  newdat <- cbind(
    newdat,
    data[(maxlag + 1):nrow(data), nolag_cols]
  )

  return(newdat)

}
