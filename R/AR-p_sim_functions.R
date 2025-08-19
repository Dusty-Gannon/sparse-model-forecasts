

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

  X <- generate_AR_covariates(model_pars$n + model_pars$holdout,
                              beta_n = model_pars$beta_n,
                              beta_p = model_pars$beta_p)

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


#' Simulate seasonal AR-p timeseries
#'
#' This function simulates an AR process from seasonal components and
#' generates a model matrix with sparse Fourier components
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
#'   n = 365,           # length of time series
#'   phi = c(0.5, 0.1), # vector of AR terms
#'   sd = 1,            # standard deviation of the innovations
#'   beta_n = 5,        # number of covariates to include
#'   beta_p = 0,        # number of lags of covariate 1 to include
#'   K = 100,           # number of fourier components to include
#'   holdout = 50
#' )
#'
#' simulate_seasonal_AR_p_timeseries(input_pars)
#'

simulate_seasonal_AR_p_timeseries <- function(input_pars = NULL){

  model_pars <- list(  # default parameters if any are not supplied
    n = 365,           # length of time series
    # phi = c(0.5, 0.1), # default vector of AR terms
    n_phi = 5,         # the number of AR terms to randomly select
    sd = 1,            # standard deviation of the innovations
    beta_n = 5,        # number of covariates to include
    beta_p = 0,        # number of lags of covariate 1 to include
    beta_sig = NULL,   # the number of significant betas. NULL defaults to all significant.
    beta_select = NULL,# the number of beta terms to keep in the model
    K = 100,           # number of fourier components to include
    holdout = 50
  )

  for(nm in names(input_pars)){
    model_pars[[nm]] <- input_pars[[nm]]
  }

  # calculate the full timeseries length to generate
  n_tot <- model_pars$n + model_pars$holdout

  # randomly draw beta's
  N_beta <- model_pars$beta_n + model_pars$beta_p
  if(is.null(model_pars$beta_sig)){
    model_pars$beta_sig = N_beta
  }
  if(is.null(model_pars$beta_select)){
    model_pars$beta_select = N_beta
  }

  beta = c(rnorm(1,0,1),
           sample(c(rnorm(model_pars$beta_sig, 0, 3),
                    rnorm(model_pars$beta_n + model_pars$beta_p -
                          model_pars$beta_sig, 0, 0.03)),
                  replace = FALSE))

  # create matrix of covariates:
  X_full <- generate_AR_covariates(n = n_tot,
                              model_pars$beta_n, model_pars$beta_p,
                              method = 'seasonal')

  X <- X_full$X
  mu <- X %*% beta
  # how many/which of the covariates did you measure?
  beta_keep <- c(1, # keep the intercept
                 1 + sort(sample(1:N_beta, model_pars$beta_select,
                                  replace = FALSE)))

  model_pars$beta <- beta[beta_keep]
  model_pars$seasonality_full <- X_full$seasonality
  model_pars$seasonality_missing <- X_full$seasonality %>%
    filter(!x %in% beta_keep) %>%
    group_by(freq) %>%
    summarize(beta = sum(beta))

  # create seasonally fluctuating mean
  # mu2 <-  1.5 * 1:n_tot/365 + 2 * cos(pi * 1:n_tot / 365) + sin(pi * 1:n_tot / 30) +
  #   0.5 * sin(pi * 1:n_tot / 10)

  model_pars$phi <- ar_param_sim(model_pars$n_phi)


  # add the AR-2 errors
  y <- ts(mu + arima.sim(list(ar = model_pars$phi), n = n_tot,
                          sd = model_pars$sd),
          frequency = 365)

  # create the model matrix using K fourier components
  X_tot <- cbind(
    X[,beta_keep],
    1:n_tot / 365,                      # trend?
    forecast::fourier(y, model_pars$K)  # sparse seasonality terms
  )

  model_pars$X <- X_tot
  model_pars$y <- y

  return(model_pars)

}

#' Generate a model matrix for AR-p simulations
#'
#' This function generates a model matrix that may contain any/all of the following:
#' - an intercept column of 1's
#' - Correlated covariate variables
#' - lagged covariate effects
#'
#' @param n length of time series
#' @param beta_n the number of covariates to be generated
#' @param beta_p the number of lags to consider in the first covariate
#'
#' @return A covariate matrix including an intercept, a trend, and any desired
#' covariates or Fourier components.
#'
#' @export
#'

generate_AR_covariates <- function(n,
                                   beta_n = 0,
                                   beta_p = 0,
                                   method = 'correlated'){

  ### Generate correlated and lagged covariate variables ###
  # To create a covariance matrix, step one is to create an orthogonal matrix,
  # which can be done using QR decomposition of an arbitrary matrix
  if(beta_n == 0){
    return(X <- matrix(rep(1, n), nrow = n))
  }

  if(method == 'correlated'){
    P <- beta_n
    Q <- qr.Q(qr(matrix(rnorm(P^2), nrow = P, ncol = P)))

    # step two is to generate a diagonal matrix with the standard deviations
    D <- diag(x = rgamma(P, shape = 2, rate = 2))

    # now use matrix multiplication to generate Sigma
    Sigma <- t(Q) %*% D %*% Q

    # to test that this is positive definite, check that all eigenvalues are positive
    # sum(eigen(Sigma)$values < 0)

    # then Cholesky decompose Sigma to multiply by independent standard normal draws
    # and create the matrix
    X <- matrix(rnorm(n = (n + beta_p) * P),
                nrow = (n + beta_p),
                ncol = P) %*% chol(Sigma)
  }

  generate_seasonal_covariate <- function(n_tot, max_season = 365, sd = 0.1){
    a <- rnorm(4)
    freq <- round(runif(2, 2, 100))
    mu <-  a[1] * 1:n_tot/max_season + a[2] * cos(2*pi * 1:n_tot / max_season) +
      a[3] * sin(2*pi * 1:n_tot / (max_season / freq[1])) +
      a[4] * sin(2*pi * 1:n_tot / (max_season / freq[2])) + rnorm(n_tot, 0, sd)
    return(list(mu = mu,
                seasonality = data.frame(freq = c(0, 1, freq),
                                         beta = a)))
  }

  if(method == 'seasonal'){
    X <- matrix(nrow = n, ncol = 0)
    seasonality = data.frame()
    for(i in 1:beta_n){
      x_i = generate_seasonal_covariate(n)
      # plot(x_i)
      X <- cbind(X, matrix(x_i$mu, nrow = n))
      seasonality <- bind_rows(seasonality,
                               mutate(x_i$seasonality, x = i+1))
    }
  }

  # generate lagged columns of beta 1
  if(beta_p >= 1){
    beta_1 <- matrix(rep(NA, (n + beta_p) * beta_p),
                     nrow = (n + beta_p),
                     ncol = beta_p)

    for(i in 1:beta_p){
      beta_1[,i] <- c(rep(NA, i), X[1:(nrow(X)-i),1])
    }

    X <- cbind(
      beta_1,
      X[,2:ncol(X)]
    )

  }

  X <- cbind(rep(1, nrow(X) - beta_p),
             X[(beta_p+1):nrow(X),])

  if(method == 'seasonal') return(list(X = X, seasonality = seasonality))
  return(list(X = X))
}



