

#' Fit AR-p_beta model
#'
#' This function fits an AR-p_beta model on a simulated dataset with the option
#' to fit both a regularized and non regularized model version
#'
#' @param model_pars a list of parameters describing the AR-p time series as well
#' as a matrix of covariates, vectors of beta and phi parameters, and a simulated
#' time series. Must include:
#'    + n: length of time series
#'    + p: number of AR lags to consider
#'    + tau_0: the prior guess for tau0
#'    + X: matrix of covariates arranged with the first column as the intercept,
#'    columns 2:(nlags covariate + 1) are lagged versions of covariate 1, and the
#'    remaining columns are other covariates.
#'    + y: response timeseries
#'    + holdout: number of observation to hold out for model evaluation
#' @param fit_nr (default = TRUE) - should a non-regularized model also be fit to the data?
#'
#' @return A list including model input parameters and model fit objects for any models run
#'
#' @export
#'


fit_ARp_beta_model <- function(model_pars, iter = 2000, warmup = 1000,
                               fit_nr = TRUE){
  # compile stan model
  # arp_r <- rstan::stan_model("Stan/AR-p_FHS-p-beta.stan")
  arp_r <- rstan::stan_model("Stan/AR-p_err_FHS-p-beta2.stan")

  # compile data (see Stan file for descriptions of each input)
  datlist <- list(
    N = model_pars$n,
    P0 = 1,
    P = ncol(model_pars$X)-1,
    p = model_pars$p,
    y = model_pars$y[1:model_pars$n],
    X_alpha = matrix(model_pars$X[1:model_pars$n, 1],
                     ncol = 1),
    X_beta = model_pars$X[1:model_pars$n, -1],
    tau0 = model_pars$tau_0,
    slab_scl = 1,
    slab_df = 10
  )

  mtd = 18
  ad = 0.99

  # sample the posterior
  mfit_arp_r <- rstan::sampling(
    arp_r,
    data = datlist,
    chains = 4, cores = 4,
    iter = iter, warmup = warmup,
    control = list(adapt_delta = 0.99,
                   max_treedepth = mtd)
  )


  # fitting the non-regularized AR model
  if(fit_nr){

    datlist_nr <- list(
      N = model_pars$n,
      P = ncol(model_pars$X),
      p = model_pars$p,
      y = model_pars$y[1:model_pars$n],
      X = model_pars$X[1:model_pars$n, ]
    )

    # arp_nr <- rstan::stan_model("Stan/AR-p.stan")
    arp_nr <- rstan::stan_model("Stan/AR-p_err2.stan")

    mfit_arp_nr <- rstan::sampling(
      arp_nr,
      data = datlist_nr,
      chains = 4, cores = 4,
      iter = iter, warmup = warmup,
      control = list(adapt_delta = 0.99,
                     max_treedepth = mtd)
    )

    return(list(
      model_pars = model_pars,
      mfit_r = mfit_arp_r,
      mfit_nr = mfit_arp_nr
    ))
  }

  return(list(
    model_pars = model_pars,
    mfit_r = mfit_arp_r
  ))
}


#' Fit seasonal AR-p_beta model
#'
#' This function fits an AR-p_beta model on a simulated dataset with the option
#' to fit both a regularized and non regularized model version
#'
#' @param model_pars a list of parameters describing the AR-p time series as well
#' as a matrix of covariates, vectors of beta and phi parameters, and a simulated
#' time series. Must include:
#'    + n: length of time series
#'    + p: number of AR lags to consider
#'    + tau_0: the prior guess for tau0
#'    + X: matrix of covariates arranged with the first column as the intercept,
#'    columns 2:(nlags covariate + 1) are lagged versions of covariate 1, and the
#'    remaining columns are other covariates.
#'    + y: response timeseries
#'    + holdout: number of observation to hold out for model evaluation
#' @param fit_nr (default = TRUE) - should a non-regularized model also be fit to the data?
#'
#' @return A list including model input parameters and model fit objects for any models run
#'
#' @export
#'


fit_seasonal_ARp_models <- function(model_pars,
                                    fit_gauss = TRUE, fit_flat = TRUE,
                                    iter = 2000, warmup = 1000,
                                    mtd = 18, ad = 0.99){
  # compile stan model
  # arp_r <- rstan::stan_model("Stan/AR-p_FHS-p-beta.stan")
  # compile data (see Stan file for descriptions of each input)

  # list to use in the stan model
  dat_flat <-
    list(
      N = model_pars$n,
      P = ncol(model_pars$X),
      p = 15,
      y = as.double(model_pars$y)[1:model_pars$n],
      X = model_pars$X[1:model_pars$n, ],
      N_new = model_pars$holdout,
      X_new = model_pars$X[(model_pars$n + 1):nrow(model_pars$X), ]
    )

  # don't regularize the trend for the
  # gaussian fit
  dat_gauss <- c(
    dat_flat,
    list(P_0 = length(model_pars$beta) + 1)
  )

  # add in HS-specific inputs
  # tau_0 <- tau0(
  #   y = dat_flat$y,
  #   m0 = 20,
  #   M = ncol(model_pars$X) + dat_flat$p*2,
  #   N = model_pars$n,
  #   fam = "gaussian"
  # )

  dat_hs <- c(
    dat_gauss,
    list(
      tau0_beta = 0.01,
      slab_scl_beta = 0.5,
      slab_df_beta = 4,
      tau0_phi = 0.001,
      slab_scl_phi = 0.5,
      slab_df_phi = 4
    )
  )

  # fit the models
  if(!('ar_err_hs' %in% ls(envir = .GlobalEnv))){
    ar_err_hs <- rstan::stan_model("Stan/AR-p_err3_FHS_DG.stan")
  }
  hs_fit <- sampling(
    ar_err_hs,
    dat = dat_hs,
    chains = 4, cores = 4,
    iter = iter, warmup = warmup,
    control = list(adapt_delta = ad,
                   max_treedepth = mtd)
  )
  fits <- list(hs_fit = hs_fit)

  if(fit_gauss){
    if(!('ar_err_gauss' %in% ls(envir = .GlobalEnv))){
      ar_err_gauss <- rstan::stan_model("Stan/AR-p_err3_Gauss_DG.stan")
    }
    gauss_fit <- sampling(
      ar_err_gauss,
      data = dat_gauss,
      chains = 4, cores = 4,
      iter = iter, warmup = warmup,
      control = list(adapt_delta = ad,
                     max_treedepth = mtd)
    )
    fits[['gauss_fit']] <- gauss_fit
  }

  if(fit_flat){
    if(!('ar_err_flat' %in% ls(envir = .GlobalEnv))){
      ar_err_flat <- rstan::stan_model("Stan/AR-p_err3_Flat_DG.stan")
    }
    flat_fit <- sampling(
      ar_err_flat,
      data = dat_flat,
      chains = 4, cores = 4
    )
    fits[['flat_fit']] <- flat_fit
  }

  return(list(
      model_pars = model_pars,
      fits = fits
    ))

}




#' Fit arima model
#'
#' This function fits an AR-p_beta model on a simulated dataset with the option
#' to fit both a regularized and non regularized model version
#'
#' @param model_pars a list of parameters describing the AR-p time series as well
#' as a matrix of covariates, vectors of beta and phi parameters, and a simulated
#' time series. Must include:
#'    + n: length of time series
#'    + p: number of AR lags to consider
#'    + tau_0: the prior guess for tau0
#'    + X: matrix of covariates arranged with the first column as the intercept,
#'    columns 2:(nlags covariate + 1) are lagged versions of covariate 1, and the
#'    remaining columns are other covariates.
#'    + y: response timeseries
#'    + holdout: number of observation to hold out for model evaluation
#' @return the forecast rmse
#'
#' @export
#'


fit_arima_model <- function(model_pars){

  # determine how many covariates to include:
  n_beta <- length(model_pars$beta)-1
  n <- (length(model_pars$y)-model_pars$holdout)
  if(n_beta == 0){X = NA}
  if(n_beta > 0){X = model_pars$X[, 3:(3+n_beta-1)]}

 #fit the models

  if(is.na(X)){
    # fit_ar <- forecast::tbats(model_pars$y[1:n])
    fit_ar <- forecast::auto.arima(model_pars$y[1:n])
    ar_for <- forecast::forecast(fit_ar, h = model_pars$holdout)
    ar_beta <- NA
  }

  if(!is.na(X)){
    fit_ar <- forecast::auto.arima(model_pars$y[1:n],
                                   xreg = X[1:n,])
    ar_for <- forecast::forecast(fit_ar, xreg = X[(n+1):nrow(X),])
    ar_beta <- fit_ar$coef[(length(fit_ar$coef)-n_beta+1):length(fit_ar$coef)]
  }

  # calculate forecast RMSE
  rmse <- sqrt(mean((ar_for$mean - model_pars$y[(n+1):length(model_pars$y)])^2))

  return(list(rmse = rmse,
              ar_forecast = ar_for))

}


