
#' Unpack AR-p model fit and predict holdout set
#'
#' This function takes the output of two model fits, one regularized and one
#' not regularized, run on simulated AR-p data and summarizes the posterior
#' parameter estimates, predicts the holdout set, and calculates rmse. It is a
#' wrapper for extract_predict_fits
#'
#'
#' @param fits a list of model_pars and stanfit objects, or if model pars is
#' provided separately, this can just be a stanfit object.
#' @param model_pars a list of parameters that was used to run models
#'
#' @return A list including model input parameters and a model fit list for each
#' model that contains: posterior parameter estimates, forcasts of held out
#' observations, and rmses of model forecasts.
#'
#' @export
#'

unpack_ARp_fit <- function(fits, model_pars = NULL){


  model_pars <- fits$model_pars

  if('mfit_nr' %in% names(fits)){
    mod_fit_r <- extract_predict_fit(fit = fits$mfit_r, model_pars, TRUE)
    mod_fit_nr <- extract_predict_fit(fit = fits$mfit_nr, model_pars, FALSE)


    return(list(
      model_pars = model_pars,
      mod_fit_r = mod_fit_r,
      mod_fit_nr = mod_fit_nr
    ))
  }

  mod_fit_r <- extract_predict_fit(fits$mfit_r, model_pars, TRUE)

  return(list(
    model_pars = fits$model_pars ,
    mod_fit_r = mod_fit_r
  ))

}



#' Extract parameter estimates from a stanfit object and predict holdout set
#'
#' This function calculate parameter RMSEs by comparing to know values used in
#' simulation. Forecast a held out set of observations and calculate a forecast
#' RMSE.
#'
#' @param fit a stanfit object from a simulated AR-p dataset.
#'
#' @param reg a boolean indicating if the stanfit is from a regularized (TRUE)
#' or un-regularized (FALSE) model.
#'
#' @return a list including parameter estimates, parameter RMSE and forecast
#' RMSEs, names of parameters with bad fits according to stanpsum function,
#' and a forecast for the holdout observations.
#'
#' @export
#'

extract_predict_fit <- function(fit, model_pars, reg = TRUE){
  # Extract posterior estimates
  if(reg){
    alpha_post <- rstan::extract(fit, pars = "alpha")$alpha
    beta_post <- cbind(alpha_post,
                       rstan::extract(fit, pars = "beta")$beta)
  }else{
    beta_post <- rstan::extract(fit, pars = "beta")$beta

  }

  phi_post <- rstan::extract(fit, pars = "phi")$phi
  sigma_post <- rstan::extract(fit, pars = "sigma")$sigma

  beta_hat <- data.frame(
    mean = apply(beta_post, 2, mean),
    median = apply(beta_post, 2, median),
    min = apply(beta_post, 2, min),
    max = apply(beta_post, 2, max),
    low = apply(beta_post, 2, quantile, probs = 0.025),
    high = apply(beta_post, 2, quantile, probs = 0.975),
    q0.01 = apply(beta_post, 2, quantile, probs = 0.01),
    q0.05 = apply(beta_post, 2, quantile, probs = 0.05),
    q0.1 = apply(beta_post, 2, quantile, probs = 0.1),
    q0.9 = apply(beta_post, 2, quantile, probs = 0.9),
    q0.95 = apply(beta_post, 2, quantile, probs = 0.95),
    q0.99 = apply(beta_post, 2, quantile, probs = 0.99)

  )

  phi_hat <- data.frame(
    mean = apply(phi_post, 2, mean),
    median = apply(phi_post, 2, median),
    min = apply(phi_post, 2, min),
    max = apply(phi_post, 2, max),
    low = apply(phi_post, 2, quantile, probs = 0.025),
    high = apply(phi_post, 2, quantile, probs = 0.975),
    q0.01 = apply(phi_post, 2, quantile, probs = 0.01),
    q0.05 = apply(phi_post, 2, quantile, probs = 0.05),
    q0.1 = apply(phi_post, 2, quantile, probs = 0.1),
    q0.9 = apply(phi_post, 2, quantile, probs = 0.9),
    q0.95 = apply(phi_post, 2, quantile, probs = 0.95),
    q0.99 = apply(phi_post, 2, quantile, probs = 0.99)
  )
  # reverse the estimates to match input vector
  phi_hat <- phi_hat[nrow(phi_hat):1,]

  sigma_hat <- data.frame(
    mean = mean(sigma_post),
    median = median(sigma_post),
    min = min(sigma_post),
    max = max(sigma_post),
    low = quantile(sigma_post, probs = 0.025),
    high = quantile(sigma_post, probs = 0.975),
    q0.01 = quantile(sigma_post, probs = 0.01),
    q0.05 = quantile(sigma_post, probs = 0.05),
    q0.1 = quantile(sigma_post, probs = 0.1),
    q0.9 = quantile(sigma_post, probs = 0.9),
    q0.95 = quantile(sigma_post, probs = 0.95),
    q0.99 = quantile(sigma_post, probs = 0.99)

  )

  par_ests <- list(beta_hat = beta_hat,
                   phi_hat = phi_hat,
                   sigma_hat = sigma_hat)

  # forecast the held-out observations

  y_rep <- rstan::extract(fit, pars = "y_rep")$y_rep
  draws <- nrow(beta_post)

  # matrix of draws from the posterior-predictive distribution
  post_preds <- matrix(nrow = draws, ncol = model_pars$n + model_pars$holdout)

  # fill in first p observations that are considered fixed
  post_preds[, 1:model_pars$p] <- matrix(
    rep(model_pars$y[1:model_pars$p], each = draws),
    nrow = draws, ncol = model_pars$p
  )

  # fill in post. pred. draws from stan
  post_preds[, (model_pars$p + 1):model_pars$n] <- y_rep

  for(i in 1:draws){
    for(t in (model_pars$n + 1):(model_pars$n + model_pars$holdout)){
      y_past <- as.double(post_preds[i, (t - model_pars$p):(t - 1)])
      post_preds[i, t] <- model_pars$X[t, ] %*% beta_post[i, ] +
        phi_post[i, ] %*% y_past + rnorm(1, sd = sigma_post[i])
    }
  }

  forecast <- data.frame(
    time = 1:(model_pars$n + model_pars$holdout),
    y = as.double(model_pars$y),
    estim = apply(post_preds, 2, mean),
    median = apply(post_preds, 2, median),
    low = apply(post_preds, 2, quantile, probs = 0.025),
    high = apply(post_preds, 2, quantile, probs = 0.975)
  )


  # compute prediction root mean squared error for each model
  # RMSE_bayes() is a user-defined function in R/model_checking.R
  rmse_forecast = RMSE_bayes(model_pars$y[(model_pars$n + 1):(model_pars$n + model_pars$holdout)],
                             ppreds = post_preds[, (model_pars$n + 1):(model_pars$n + model_pars$holdout)])

  rmse_beta = RMSE_bayes(model_pars$beta, beta_post)
  rmse_phi = RMSE_bayes(model_pars$phi, phi_post)
  rmse_sigma = RMSE_bayes(model_pars$sigma_e, sigma_post)

  rmse <- list(rmse_forecast = rmse_forecast,
               rmse_beta = rmse_beta,
               rmse_phi = rmse_phi,
               rmse_sigma = rmse_sigma)

  bad_fits <- stan_psum(fit = fit)
  model = 'not_reg'
  if(reg) model = 'reg'

  mod_fit <- list(
    par_ests = par_ests,
    forecast = forecast,
    rmse = rmse,
    bad_fits = bad_fits,
    model = model)

  return(mod_fit)

}

#' Summarize model fit metrics from posterior
#'
#' This function takes a Stan model fit object and calculates the ratio of the
#' effective sample size for each model parameter to the number of iterations
#' and reports which parameters are less than 10%. It also reports any parameters
#' with an Rhat > 1.1 and gets the number of divergent transitions.
#'
#' @param fit a stan fit object
#'
#' @return A list containing the names of any parameters with an Rhat > 1.1 or
#' an n_eff < 10% of the number of iterations and the number of divergent transitions.
#'
#' @export
#'
#'

stan_psum <- function(fit){

  iter <- fit@stan_args[[1]]$iter
  s_init <- as.data.frame(rstan::summary(fit)$summary)
  s_init$pars <- row.names(s_init)
  # effective samples are the number of independent samples with the same estimation power as the N autocorrelated samples
  s_init$n_eff_pct <- s_init$n_eff/iter
  # 10% is often used as a threshold, below which the chains for a parameter did not properly converge
  s_init$n_eff_less10pct <- ifelse(s_init$n_eff_pct < 0.10,
                                   yes = "true", no = "false")
  s <- s_init[,c("pars","Rhat","n_eff_less10pct")]
  s <- s[!grepl('^mu', row.names(s_init)),]
  s <- s[s$Rhat > 1.1 | s$n_eff_less10pct == 'true',]
  spars <- get_sampler_params(fit, inc_warmup = FALSE)
  divtrans <- sum(sapply(spars, function(x) sum(x[,'divergent__'])))

  return(list(par_conv = s,
              divergent = divtrans)
  )

}

#' Calculate true positives rates for models
#'
#' For a model fit to simulated AR-p data, a true positive is a correctly
#' identified beta or phi parameter as non-zero. We define a true positive as
#' any parameter estimate for which 90% of the posterior mass lies above or
#' below zero, given that the true parameter is non-zero.
#'
#' @param par_post posterior estimates of the parameter values
#'
#' @param value true values of the parameters
#'
#' @param threshold the fraction of the posterior mass that must fall above or
#' below zero for a detection to count as positive. Defaults to 90%
#'
#' @param par either 'beta' or 'phi' to determine what the lower bound should
#' be for a positive detection. Defaults to 'beta'.
#'
#' @return A data frame of the true positive rate, true negative rate,
#' false positive rate, and false negative rate.
#'
#' @export
#'

calculate_true_pos_rate <- function(par_post, value, threshold = 0.9,
                                    par = 'beta'){
  min_pos = 0
  if(par == 'beta') min_pos = 1.96 * 0.05

  par_post <- rename(par_post,
         q0.025 = low, q0.975 = high)

  # calculate the boundaries of the posterior mass based on the threshold
  bottom <- paste0('q', (1-threshold))
  top <- paste0('q', (threshold))

  ests <- par_post[,c('mean', 'median', bottom, top)]
  colnames(ests) <- c('mean', 'median', 'bottom', 'top')

  ests <- ests %>%
    mutate(value = value,
           effect = case_when(abs(value) <= min_pos ~ 'zero',
                              value > min_pos ~ 'positive',
                              value < -min_pos ~ 'negative'),
           post_mass = case_when(bottom > 0 ~ 'positive',
                                 top < 0 ~ 'negative',
                                 TRUE ~ 'zero'))
  true_pos = sum((ests$effect == 'positive' & ests$post_mass == 'positive')|
                   (ests$effect == 'negative' & ests$post_mass == 'negative'))/
    sum(ests$effect != 'zero')
  false_pos = sum((ests$effect != 'negative' & ests$post_mass == 'negative')|
                    (ests$effect != 'positive' & ests$post_mass == 'positive'))/
    sum(ests$post_mass != 'zero')
  false_neg = sum((ests$effect == 'negative' & ests$post_mass != 'negative')|
                    (ests$effect == 'positive' & ests$post_mass != 'positive'))/
    sum(ests$post_mass == 'zero')

  TFP <- data.frame(true_pos = true_pos,
                    false_pos = false_pos,
                    false_neg = false_pos) %>%
      rename_with(function(x) paste0(par, '_', x))

  return( TFP )

}

#' Summarize the true positive rates of different types of parameter estimates
#' for a paired set of regularized and not regularized models run on the same
#' simulated dataset.
#'
#' @param fits A list of stanfit objects
#'
#' @param model_pars a list of the parameters used to simulate
#' them.
#'
#' @param threshold the fraction of the posterior mass that must fall above or
#' below zero for a detection to count as positive. Defaults to 90%
#'
#' @return a data frame with the true positive rates for phi and beta parameters
#' calculated for both the regularized and unregularized models.
#'
#' @export
#'

summarize_pos_rate <- function(fits, model_pars, threshold = 0.9){

  pr <- data.frame()
  for(i in 1:length(fits)){
      beta_post <- fits[[i]]$par_ests$beta_hat
      phi_post <- fits[[i]]$par_ests$phi_hat

      pos_rate_beta <- calculate_true_pos_rate(beta_post, model_pars$beta,
                                               threshold = threshold,
                                               par = 'beta')
      pos_rate_phi <- calculate_true_pos_rate(phi_post, model_pars$phi,
                                              threshold = threshold,
                                              par = 'phi')
      pr <- cbind(pos_rate_phi, pos_rate_beta) %>%
        mutate(model = fits[[i]]$model) %>%
        bind_rows(pr)
  }

  pr$TPR_threshold = threshold

  return( pr )

}
