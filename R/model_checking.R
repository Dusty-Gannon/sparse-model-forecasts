



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
      sqrt(mean(x - obs)^2)
    },
    obs = obs
  )
}









#' Quantile Residuals
#'
#' Compute quantile residuals from a glm or glmnet object
#'
#' @param obs Vector of observed values
#' @param pred Vector of predicted values
#' @param fam Error distribution family. Options are \code{"poisson", "gaussian", "neg_binomial"}.
#' @param z If FALSE, the returned residuals should be uniformly distributed
#' @param ... Extra parameters passed to distribution function to compute lower-tail
#' probabilities (e.g., \code{size} for negative binomial distribution)
#'
#' @return Vector of quantile residuals
#' @export
#'
#' @examples
#' pred <- exp(seq(1, 2, length.out = 100))
#' obs <- rpois(100, lambda = pred)
#' quantile_resids(obs, pred, fam = "poisson")
quantile_resids <- function(obs, pred, fam = "poisson", z = TRUE, ...){

  xtra_params <- list(...)

  if(fam == "poisson"){
    q <- purrr::map2_dbl(
      obs,
      pred,
      ~ runif(
        n = 1,
        min = ppois(max(c(.x - 1, 0)), .y),
        max = ppois(.x, .y)
      )
    )
    if(isTRUE(z)){
      z <- qnorm(q)
      return(z)
    } else { return(q) }
  }

  if(fam == "gaussian"){
    q <- purrr::map2_dbl(
      obs,
      pred,
      ~ pnorm(.x, .y, sd = xtra_params$sd)
    )
    if(isTRUE(z)){
      z <- qnorm(q)
      return(z)
    } else { return(q) }
  }

  if(fam == "neg_binomial"){
    q <- purrr::map2_dbl(
      obs,
      pred,
      ~ runif(
        n = 1,
        min = pnbinom(max(c(.x - 1, 0)), mu = .y, size = xtra_params$size),
        max = pnbinom(.x, mu = .y, size = xtra_params$size)
      )
    )
    if(isTRUE(z)){
      z <- qnorm(q)
      return(z)
    } else { return(q) }
  }

}




#' Posterior Predictive Samples
#'
#' Draw samples from the posterior predictive distribution (PPD) of a fitted model
#'
#' @param mu_post Matrix of posterior draws of the mean with posterior samples
#' in rows and mean parameters (multiple for a Bayesian GLM-like model) in columns
#' @param ndraws Number of samples to take from the PPD
#' @param fam Family of the PPD (currently supports \code{"poisson"} and \code{"neg_binomial"})
#' @param ... Draws from posterior distributions of extra parameters
#' (e.g., phi = neg_binomial dispersion parameter) necessary for predictive draws
#'
#' @return Matrix of \code{ndraws} (in rows) from the PPD of each of the means supplied (in columns)
#' @export
#'
#' @examples
post_predict <- function(mu_post, ndraws = NULL, fam = "poisson", ...){

  # turn extra args into a list
  xtra_params <- list(...)
  # if ndraws is left blank, set to number of posterior draws
  if(is.null(ndraws)){
    ndraws <- dim(mu_post)[1]
  }

  # sample means from posteriors
  samps <- sample(1:dim(mu_post)[1], size = ndraws)

  # create a list of posterior predictions
  if(fam == "poisson"){
    pp_draws <- purrr::map(
      samps,
      ~ rpois(
        n = dim(mu_post)[2],
        lambda = mu_post[.x, ]
      )
    )
  }

  if(fam == "neg_binomial"){
    pp_draws <- purrr::map(
      samps,
      ~ rnbinom(
        n = dim(mu_post)[2],
        mu = mu_post[.x, ],
        size = xtra_params$phi[.x]
      )
    )
  }

  # convert to an array of posterior predictions
  pp_mat <- t(simplify2array(pp_draws))

  # make some column names
  index <- 1:dim(mu_post)[2]
  nms <- paste("y_rep", as.character(index), sep = "_")
  colnames(pp_mat) <- nms
  return(pp_mat)

}






#' Posterior predictive checks
#'
#' @param y_rep Matrix of S predicted values sampled from the
#'  posterior predictive distribution for each of N observed values
#' @param obs Vector of N observed values
#' @param plot Should a PPC plot be returned
#' @param pred_level Level of posterior predictive interval
#'
#' @return Either a list or the single value of the proportion of observed values captured by
#' the posterior predictive intervals (should be around pred_level given a good fit)
#' @export
#'
#' @examples
#'
ppcs <- function(y_rep, obs, plot = T, pred_level = 0.95){
  # build df
  df_ppc <- data.frame(
    y = obs,
    mean = colMeans(y_rep),
    pred_low = apply(y_rep, 2, quantile, probs = (1 - pred_level)/2),
    pred_high = apply(y_rep, 2, quantile, probs = 1 - (1 - pred_level)/2)
  )

  # order the df by ascending mean pred
  df_ppc_ord <- df_ppc[order(df_ppc$mean), ]
  df_ppc_ord$x <- 1:nrow(df_ppc_ord)
  df_ppc_ord$original_obs_id <- order(df_ppc$y)

  # create list with plot and summary stat
  if(isTRUE(plot)){
    ppc_plot <- ggplot2::ggplot(data = df_ppc_ord, aes(x = x))+
      geom_errorbar(
        aes(ymin = pred_low, ymax = pred_high),
        color = "grey",
        width = 0,
        size = 1.5
      )+
      geom_point(aes(y = y), color = "red", size = 0.5)+
      theme_bw()+
      xlab("")+
      ylab("Observed (red) and predicted fecundity")

    # calculate proportion of observations captured
    prop_captured <- mean(df_ppc$y >= df_ppc$pred_low & df_ppc$y <= df_ppc$pred_high)

    # find which observations fall outside the predictive intervals
    outside_obs <- which(df_ppc$y <= df_ppc$pred_low | df_ppc$y >= df_ppc$pred_high)

    # message about results
    print(paste(
      "Proportion of observations captured in ",
      pred_level*100,
      "% posterior predictive intervals: ",
      round(prop_captured, 2),
      sep = ""
    ))

    return(list(
      plot = ppc_plot,
      prop_captured = prop_captured,
      df = df_ppc_ord
    ))
  }

  if(isFALSE(plot)){
    return(
      mean(df_ppc$y >= df_ppc$pred_low & df_ppc$y <= df_ppc$pred_high)
    )
  }
}







#' Bayesian Quantile Residuals
#'
#' Compute Bayesian quantile or probability residuals using observed data and the
#' posterior predictive distribution. A probability residual is defined as
#' \eqn{p_i = P(\tilde y_i | y_i | \bf{y})} where \eqn{\tilde y_i} is a realization from the
#' posterior predictive distribution. The corresponding quantile residual, \eqn{q_i}, is defined
#' as \eqn{q_i = \Phi^{-1}(p_i)}, where \eqn{\Phi} is the standard normal cdf.
#'
#'
#' @param obs Observed data in vector of length \eqn{n}
#' @param y_rep Draws from the posterior predictive distribution in a matrix of \eqn{n} columns
#' and \code{ndraws} rows.
#' @param z If true, quantile residuals are returned. Otherwise, probability residuals are returned.
#'
#' @return Vector of residuals
#' @export
#'
#' @examples
#'
bayes_qresids <- function(obs, y_rep, z = TRUE){

  # index for observations
  index <- 1:length(obs)

    # compute P(y < y_rep | y)
    p <- purrr::map_dbl(
          index,
          ~ mean(y_rep[, .x] < obs[.x])
        )

    # convert to z - quantiles
    if(isTRUE(z)){
      z <- qnorm(p)
      return(z)
    } else { return(p) }

}






#' Empirical Kullback-Leibler Divergence
#' This function computes the empirical KL divergence using the estimator of Perez-Cruz (2008)
#'
#' @param prior_samps Vector of samples from the prior (or reference) distribution
#' @param post_samps Vector of samples from the posterior distribution or distribution of interest
#'
#' @return scalar estimate of the KL-divergence
#' @export
#'
#' @examples
#' # Theoretical KL-divergence for p(x) = exp(2) and q(x) = exp(1) is 0.193 (ish)
#'
#' prior_samps <- rexp(500)
#' post_samps <- rexp(500, rate = 2)
#'
#' kl_divergence(prior_samps, post_samps)
#'
kl_divergence <- function(prior_samps, post_samps){

  # order the samples smallest to largest
  prior_samps <- prior_samps[order(prior_samps)]
  post_samps <- post_samps[order(post_samps)]

  # get a suitable epsilon
  thresh <- min(c(
    prior_samps[-1] - prior_samps[-length(prior_samps)],
    post_samps[-1] - post_samps[-length(post_samps)]
  ))
  epsilon <- thresh/2

  # deal with boundary conditions
  x_0 <- min(prior_samps, post_samps) - 2 * abs(min(prior_samps) - min(post_samps))
  x_np1 <- max(prior_samps, post_samps) + 2 * abs(max(prior_samps) - max(post_samps))

  # create empirical cdf function
  ecdf <- function(x, obs){
    n <- length(obs)
    purrr::map_dbl(
      x,
      ~ mean(obs < .x) + (1/n) * 0.5 * sum(.x == obs)
    )
  }

  # create piecewise linear functions based on Perez-Cruz (2008)
  Pc <- function(x, samps, x_0, x_np1){

    # tack on x_np1
    samps2 <- c(samps, x_np1)

    # create empirical probabilities
    P_vec <- ecdf(x = samps2, obs = samps)

    # create a list of lines
    pwlines <- vector(mode = "list", length = length(samps2))
    a_1 <- (P_vec[1] - 0)/(samps2[1] - x_0)
    b_1 <- (-1) * (P_vec[1] * x_0)/(samps2[1] - x_0)
    pwlines[[1]] <- c(a = a_1, b = b_1)

    for(i in 2:length(samps2)){
      # slope and intercept of line passing through two points
      a_i <- (P_vec[i] - P_vec[i - 1]) / (samps2[i] - samps2[i - 1])
      b_i <- (P_vec[i - 1] * samps2[i] - P_vec[i] * samps2[i - 1]) / (samps2[i] - samps2[i - 1])

      pwlines[[i]] <- c(a = a_i, b = b_i)
    }

    # write the piecewise linear function
    pwl_func <- function(x, lookup = pwlines){

      if(x < x_0){
        return(0)
      } else if(x > x_np1){
        return(1)
      } else{
        index <- min(which(x < samps2))
        return(
          lookup[[index]][1] * x + lookup[[index]][2]
        )
      }

    }

    # now return Pc function
    purrr::map_dbl(
      x,
      ~ pwl_func(.x)
    )

  }

  # now use Pc function to compute KLD estimator (eq. 4 in Perez-Cruz, 2008)
  numers <- Pc(post_samps, samps = post_samps, x_0, x_np1) -
    Pc(post_samps - epsilon, samps = post_samps, x_0, x_np1)
  denoms <- Pc(post_samps, samps = prior_samps, x_0, x_np1) -
    Pc(post_samps - epsilon, samps = prior_samps, x_0, x_np1)
  KLD_hat <- mean(
    log(numers / denoms),
    na.rm = T
  )

  return(KLD_hat - 1)

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


















