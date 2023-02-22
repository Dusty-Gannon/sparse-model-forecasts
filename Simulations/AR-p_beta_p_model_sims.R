#############################################################
# This function simulates an AR process in which some lags and
# explanatory variables are important, but not others. We use
# regularizing priors to shrink the unimportant parameters
# towards zero.
#############################################################

# libraries
library(rstan)
library(here)
library(ggplot2)
devtools::load_all()

#### Simulating the data ####
model_params <- list(
  n = 300,      # length of time series
  p  = 16,      # number of AR lags to consider
  beta_p = 5,   # number of beta lags to consider in lagged covariate (beta_1)
  beta_n = 45,  # number of additional covariates to include
  b0 = 0.5,     # intercept
  # coefficient vector for b1 lags:
  b1 = c(-0.5, rep(0, 4), -0.5),
  beta_important = c(0.8, -0.2), # number of non-zero elements to include
  phi = c(0.6, rep(0, 4), -0.2, rep(0, 9), 0.3),
  holdout = 50
)

simulate_AR_p_beta_p_timeseries <- function(model_params){

  ### generate the model matrix with some correlated variables ###
  # To create a covariance matrix, step one is to create an orthogonal matrix,
  # which can be done using QR decomposition of an arbitrary matrix
  P <- model_params$beta_n + 1
  Q <- qr.Q(qr(matrix(rnorm(P^2), nrow = P, ncol = P)))

  # step two is to generate a diagonal matrix with the standard deviations
  D <- diag(x = rgamma(P, shape = 2, rate = 2))

  # now use matrix multiplication to generate Sigma
  Sigma <- t(Q) %*% D %*% Q

  # to test that this is positive definite, check that all eigenvalues are positive
  # sum(eigen(Sigma)$values < 0)

  # then Cholesky decompose Sigma to multiply by independent standard normal draws
  # and create the matrix
  X <- matrix(rnorm(n = model_params$n * P),
              nrow = model_params$n, ncol = P) %*% chol(Sigma)


  # generate lagged columns of beta 1
  beta_1 <- matrix(rep(NA, model_params$n * model_params$beta_p),
                   nrow = model_params$n, ncol = model_params$beta_p)

  for(i in 1:model_params$beta_p){
    beta_1[,i] <- c(rep(NA, i), X[1:(nrow(X)-i),1])
  }

  X <- cbind(
    rep(1, model_params$n),
    beta_1,
    X[,2:ncol(X)]
  )

  X <- X[(model_params$beta_p+1):model_params$n,]

  # reassign n because of days lost to beta_1 lag:
  n <- nrow(X)

  beta <- c(model_params$b0, model_params$b1, model_params$beta_important)
  beta <- c(beta, rep(0, ncol(X)-length(beta)))

  # mean of the process
  mu <- as.double(X %*% beta)
  sigma_e <- 1

  # simulate the AR process
  y <- arima.sim(
    n = n,
    model = list(ar = model_params$phi),
    mean = mu,
    sd = sigma_e
  )

#### Fit the model, holding out last 10 observations for forecasting ####

  # compile stan model
  arp_r <- stan_model(here("Stan/AR-p_FHS-p-beta.stan"))

  # compute prior guess for tau0 based on a guess of 5
  #  non-zero coefficients
  #  see ?tau0() for documentation
  tau_0 <- tau0(
    y = y[1:(n - model_params$holdout)],
    m0 = 5,
    M = ncol(X) + length(model_params$phi),
    N = n - model_params$holdout,
    fam = "gaussian"
  )

  # compile data (see Stan file for descriptions of each input)
  datlist <- list(
    N = n - model_params$holdout,
    P0 = 1,
    P = ncol(X)-1,
    p = length(model_params$phi),
    y = y[1:(n - model_params$holdout)],
    X_alpha = matrix(X[1:(n - model_params$holdout), 1], ncol = 1),
    X_beta = X[1:(n - model_params$holdout), -1],
    tau0 = tau_0,
    slab_scl = 1,
    slab_df = 10
  )

  # sample the posterior
  mfit_arp_r <- sampling(
    arp_r,
    data = datlist,
    chains = 3, cores = 3
  )

  # visualize the results
  rawdata_df <- data.frame(
    time = 1:n,
    response = as.double(y)
  )

  # plot for raw time series
  plot_ts <- ggplot(data = rawdata_df, aes(x = time, y = response)) +
    geom_line() +
    theme_classic() +
    ggtitle("Raw data")

  # plots for posterior estimates and CIs compared to simulated values
  plot_beta <- ggplot(data = df_beta, aes(x = param)) +
    geom_errorbar(aes(ymin = low, ymax = high), width = 0) +
    geom_point(aes(y = beta_hat)) +
    geom_point(aes(y = beta_true), shape = 2, color = "steelblue", size = 3) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    annotate("point", x = 40, y = 2, shape = 2, color = "steelblue", size = 3) +
    annotate("text", label = "True value", x = 41, y = 2, hjust = 0) +
    annotate("point", x = 40, y = 1.5) +
    annotate("text", x = 41, y = 1.5, label = "Estimate", hjust = 0) +
    ylab("value") +
    ggtitle("Covariate coefficients")

  plot_phi <- ggplot(data = df_phi, aes(x = param)) +
    geom_errorbar(aes(ymin = low, ymax = high), width = 0) +
    geom_point(aes(y = phi_hat)) +
    geom_point(aes(y = phi_true), shape = 2, color = "steelblue", size = 3) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    annotate("point", x = 2, y = 0.5, shape = 2, color = "steelblue", size = 3) +
    annotate("text", label = "True value", x = 3, y = 0.5, hjust = 0) +
    annotate("point", x = 2, y = 0.375) +
    annotate("text", x = 3, y = 0.375, label = "Estimate", hjust = 0) +
    ylab("value") +
    ggtitle("AR coefficients")

  ## save figure ##
  lo <- rbind(
    c(1, 1, 2, 2),
    c(3, 3, 3, 3)
  )

  gridExtra::grid.arrange(plot_phi, plot_ts, plot_beta, layout_matrix = lo)
  p <- gridExtra::grid.arrange(plot_phi, plot_ts, plot_beta, layout_matrix = lo)
  print(p)

  return(list(fit = mfit_arp_r, p = p))
}


extract_AR_p_beta_p_mod_ests <- function(model_params, fit_arp_r){

  n <- model_params$n - model_params$beta_p
  beta <- c(model_params$b0, model_params$b1, model_params$beta_important)
  beta <- c(beta, rep(0, model_params$beta_n + model_params$beta_p + 1 - length(beta)))

    #### bind parameter posteriors to truth ####
  post <- Reduce(cbind, rstan::extract(fit_arp_r, pars = c("alpha", "beta", "phi")))

  # create data frames for plotting
  df_pars <- data.frame(
    param = c(paste0("beta", 0:(length(beta) - 1)),
              paste0("phi", 1:model_params$p)),
    beta_true = c(beta, model_params$phi),
    beta_hat = apply(post, 2, mean),
    low = apply(post, 2, quantile, probs = 0.025),
    high = apply(post, 2, quantile, probs = 0.975),
    n = n
  )

  return(df_pars)

}

plot_AR_p_beta_p_sims <- function(model_params, mfit_arp_r){

  n <- model_params$n - model_params$beta_p

  return(p)
}






