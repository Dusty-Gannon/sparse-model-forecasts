


#' Simulate AR-p data with lagged covariates
#'
#' This function simulates an AR process in which some lags and explanatory
#' variables, as well as their lags are important, but not others. Generates
#' a matrix of correlated covariates with n-lags of the first covariate.
#'
#' @param input_pars a list of parameters describing the AR-p model to be generated.
#' Default values for all parameters are embedded in the function, so only
#' supply parameters you wish to change.
#'   + n (300): length of time series to be simulated
#'   + p (16): number of AR lags to consider
#'   + beta_p (5): number of lags to consider in lagged covariate
#'   + beta_n (45): number of additional covariates to include
#'   + b0 (0.5): intercept
#'   + n_phi (3): number of non-zero autoregressive terms
#'   + n_lags (2): number of non-zero lags in lagged covariate
#'   + n_beta (2): number of non-zero covariate parameters in addition to lagged covariate
#'   + non_zero_coef_guess (0): guess for the number of non-zero coefficients
#'   + holdout (50): number of observation to hold out for model evaluation
#' @param draw_beta can be 'near_zero' or 'zero'. Indicates if the not-significant
#' beta parameters should be drawn from a normal distribution that results in numbers
#' that are near zero, or if they should just be set to zero
#' @param draw_phi can be 'near_zero' or 'zero'. Indicates if the not-significant
#' phi parameters should be drawn from a normal distribution that results in numbers
#' that are near zero, or if they should just be set to zero.
#' @param burnin number of steps to run the AR timeseries simulation before using numbers
#'
#' @return A list including input parameters or any default values used as well
#' as a matrix of covariates, vectors of beta and phi parameters, and a simulated
#' time series.
#'
#' @export
#'
#' @examples
#'
#' input_pars <- list(
#'   n = 300,      # length of time series
#'   p  = 16,      # number of AR lags to consider
#'   beta_p = 5,   # number of beta lags to consider in lagged covariate (beta_1)
#'   beta_n = 45,  # number of additional covariates to include
#'   b0 = 0.5,     # intercept
#'   n_phi = 3,    # number of non-zero autoregressive terms
#'   n_lags = 2,   # number of non-zero lags in covariate
#'   n_beta = 2,   # number of non-zero covariate parameters
#'   non_zero_coef_guess = 5, # guess for the number of non-zero coefficients
#'   sigma_e = 1,  #standard deviation of the innovations
#'   holdout = 50
#' )
#'
#' simulate_AR_p_beta_p_timeseries(input_pars)
#'

simulate_AR_p_beta_p_timeseries <- function(input_pars = NULL,
                                            draw_beta = 'near_zero',
                                            draw_phi = 'zero',
                                            burnin = NULL){

  model_pars <- list( # default parameters if any are not supplied
    n = 300,      # length of time series
    p  = 10,      # number of AR lags to consider
    beta_p = 5,   # number of beta lags to consider in lagged covariate (beta_1)
    beta_n = 45,  # number of additional covariates to include
    b0 = 0,       # intercept
    sar = 0.3,    # seasonal AR term
    S = 10,       # length of AR seasonal cycle
    n_phi = 3,    # number of non-zero autoregressive terms
    n_lags = 2,   # number of non-zero lags in covariate
    n_beta = 2,   # number of non-zero covariate parameters
    non_zero_coef_guess = 5, # guess for the number of non-zero coefficients
    sigma_e = 1,  #standard deviation of the innovations
    holdout = 50
  )

  for(nm in names(input_pars)){
    model_pars[[nm]] <- input_pars[[nm]]
  }

  # draw parameters from distributions:
  if(draw_beta == 'near_zero'){
   model_pars$beta = c(model_pars$b0,
                      sample(c(rnorm(model_pars$n_lags, 0, 3),
                               rep(0, model_pars$beta_p - model_pars$n_lags)),
                             replace = FALSE ),
                      sample(c(rnorm(model_pars$n_beta, 0, 3),
                               rnorm(model_pars$beta_n - model_pars$n_beta, 0, 0.03)),
                             replace = FALSE))
  }

  if(draw_beta == 'zero'){
   model_pars$beta = c(model_pars$b0,
                      sample(c(rnorm(model_pars$n_lags, 0, 3),
                               rep(0, model_pars$beta_p - model_pars$n_lags)),
                             replace = FALSE ),
                      sample(c(rnorm(model_pars$n_beta, 0, 3),
                             rep(0, model_pars$beta_n - model_pars$n_beta)),
                             replace = FALSE))
  }

  # do not allow a holdout size of greater than 30% of the data
  # if(model_pars$holdout > 0.3 * model_pars$n){
  #   model_pars$holdout = floor(0.3) * model_pars$n
  # }

  ### generate the model matrix with some correlated variables ###
  # To create a covariance matrix, step one is to create an orthogonal matrix,
  # which can be done using QR decomposition of an arbitrary matrix
  P <- model_pars$beta_n + 1
  Q <- qr.Q(qr(matrix(rnorm(P^2), nrow = P, ncol = P)))

  # step two is to generate a diagonal matrix with the standard deviations
  D <- diag(x = rgamma(P, shape = 2, rate = 2))

  # now use matrix multiplication to generate Sigma
  Sigma <- t(Q) %*% D %*% Q

  # to test that this is positive definite, check that all eigenvalues are positive
  # sum(eigen(Sigma)$values < 0)

  # then Cholesky decompose Sigma to multiply by independent standard normal draws
  # and create the matrix
  X <- matrix(rnorm(n = (model_pars$n + model_pars$beta_p + model_pars$holdout) * P),
              nrow = (model_pars$n + model_pars$beta_p + model_pars$holdout),
              ncol = P) %*% chol(Sigma)

  # generate lagged columns of beta 1
  beta_1 <- matrix(rep(NA, (model_pars$n + model_pars$beta_p + model_pars$holdout) *
                         model_pars$beta_p),
                   nrow = (model_pars$n + model_pars$beta_p + model_pars$holdout),
                   ncol = model_pars$beta_p)

  for(i in 1:model_pars$beta_p){
    beta_1[,i] <- c(rep(NA, i), X[1:(nrow(X)-i),1])
  }

  X <- cbind(
    rep(1, (model_pars$n + model_pars$beta_p + model_pars$holdout)),
    beta_1,
    X[,2:ncol(X)]
  )

  X <- X[(model_pars$beta_p+1):nrow(X),]


  # Draw phi parameters
  if(draw_phi == 'near_zero'){
     phi <- c(rnorm(model_pars$p, 0, 0.03))
  }
  if(draw_phi == 'zero'){
     phi <- c(rep(0, model_pars$p))
  }

  phi[1:model_pars$n_phi] <- ar_param_sim(model_pars$n_phi)

  # simulate the process
  sar_sim <- sim_sARp(n = model_pars$n + model_pars$holdout, ar = phi,
                      sar = model_pars$sar, S = model_pars$S,
                      sd = model_pars$sigma_e,
                      X = X, beta = model_pars$beta,
                      lagged_beta = model_pars$beta_p,
                      burnin = burnin)

  model_pars$phi <- sar_sim$phi
  model_pars$p <- length(model_pars$phi)
  # compute prior guess for tau0 based on a guess of number of
  #  non-zero coefficients
  #  see ?tau0() for documentation
  tau_0 <- tau0(
    y = sar_sim$y[1:model_pars$n],
    m0 = model_pars$non_zero_coef_guess,
    M = ncol(X) + model_pars$p*2,
    N = model_pars$n,
    fam = "gaussian"
  )

  model_pars <- c(model_pars,
                  list(X = X,
                       y = sar_sim$y,
                       tau_0 = tau_0
                  ))

  return(model_pars)

}




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
  ad = 0.95
  if(model_pars$n <=100) ad = 0.99
#  if(model_pars$n < 150) ad = 0.95

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







