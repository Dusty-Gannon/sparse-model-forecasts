#############################################################
# This script simulates an AR process in which some lags and
# explanatory variables are important, but not others. We use
# regularizing priors to shrink the unimportant parameters
# towards zero.
#############################################################

# libraries
library(rstan)
library(here)
library(ggplot2)
library(devtools)
devtools::load_all()


#### Simulating the data ####

  # length of time series
  n <- 150

  # time series parameters, constrained for stationarity
  phi <- c(0.6, rep(0, 4), -0.2, rep(0, 9), 0.3)

  # coefficient vector with two non-zero elements and the intercepts
  beta <- c(0.5, 1, 2, runif(48, min=0, max=0.1))

  ### generate the model matrix with some correlated variables ###

  # To create a covariance matrix, step one is to create an orthogonal matrix,
  # which can be done using QR decomposition of an arbitrary matrix
  P <- length(beta) - 1
  Q <- qr.Q(qr(matrix(rnorm(P^2), nrow = P, ncol = P)))

  # step two is to generate a diagonal matrix with the standard deviations
  D <- diag(x = rgamma(P, shape = 2, rate = 2))

  # now use matrix multiplication to generate Sigma
  Sigma <- t(Q) %*% D %*% Q

  # to test that this is positive definite, check that all eigenvalues are positive
  # sum(eigen(Sigma)$values < 0)

  # then Cholesky decompose Sigma to multiply by independent standard normal draws
  # and create the matrix
  X <- cbind(
    rep(1, n),
    matrix(rnorm(n = n * P), nrow = n, ncol = P) %*% chol(Sigma)
  )

  # mean of the process
  mu <- as.double(X %*% beta)
  sigma_e <- 2

  # simulate the AR process
  y <- arima.sim(
    n = n,
    model = list(ar = phi),
    mean = mu,
    sd = sigma_e
  )

#### Fit the model, holding out last 10 observations for forecasting ####

  # compile stan model
  arp_r <- stan_model(here("Stan/AR-p_FHS-p-beta.stan"))

  # number of observation to hold out for forecast testing
  holdout <- 50

  # compute prior guess for tau0 based on a guess of 5
  #  non-zero coefficients
  #  see ?tau0() for documentation
  tau_0 <- tau0(
    y = y[1:(n - holdout)],
    m0 = 5,
    M = P + 20,
    N = n - holdout,
    fam = "gaussian"
  )

  # compile data (see Stan file for descriptions of each input)
  datlist <- list(
    N = n - holdout,
    P0 = 1,
    P = P,
    p = 20,
    y = y[1:(n - holdout)],
    X_alpha = matrix(X[1:(n - holdout), 1], ncol = 1),
    X_beta = X[1:(n - holdout), -1],
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


#### Compare parameter posteriors to truth ####
  alpha_post_r <- rstan::extract(mfit_arp_r, pars = "alpha")$alpha
  beta_post_r <- rstan::extract(mfit_arp_r, pars = "beta")$beta
  phi_post_r <- rstan::extract(mfit_arp_r, pars = "phi")$phi

  # create data frames for plotting
  df_beta <- data.frame(
    param = paste0("beta", 1:(length(beta) - 1)),
    beta_true = beta[-1],
    beta_hat = apply(beta_post_r, 2, mean),
    low = apply(beta_post_r, 2, quantile, probs = 0.025),
    high = apply(beta_post_r, 2, quantile, probs = 0.975)
  )

  df_phi <- data.frame(
    param = paste0("phi", 1:ncol(phi_post_r)),
    # note that the ordering of phi is reversed in the stan model
    phi_true = c(phi, rep(0, 4))[20:1],
    phi_hat = apply(phi_post_r, 2, mean),
    low = apply(phi_post_r, 2, quantile, probs = 0.025),
    high = apply(phi_post_r, 2, quantile, probs = 0.975)
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
  png(
    filename = here("Figures/param_checks_AR-p.png"),
    width = 2700,
    height = 2100,
    res = 300, units = "px"
  )
  gridExtra::grid.arrange(plot_phi, plot_ts, plot_beta, layout_matrix = lo)
  dev.off()




#### Compare forecasting to a non-regularized fitted model #####

  # fitting the non-regularized AR model Gaussian priors
  datlist_nr <- list(
    N = n - holdout,
    P = P + 1,
    p = 20,
    y = y[1:(n - holdout)],
    X = X[1:(n - holdout), ]
  )

  arp_nr <- stan_model(here("Stan/AR-p.stan"))

  mfit_arp_nr <- sampling(
    arp_nr,
    data = datlist_nr,
    chains = 3,
    cores = 3
  )


#### Compare forecasting to a model with flat priors (no priors) #####

  # fitting the non-regularized AR model with flat priors
  datlist_fl <- list(
    N = n - holdout,
    P = P + 1,
    p = 20,
    y = y[1:(n - holdout)],
    X = X[1:(n - holdout), ]
  )

  arp_fl <- stan_model(here("Stan/AR-p_flat.stan"))

  mfit_arp_fl <- sampling(
    arp_fl,
    data = datlist_fl,
    chains = 3,
    cores = 3
  )





  # COMPARISON STARTS

  # forecast the held-out observations
  beta_post_nr <- rstan::extract(mfit_arp_nr, pars = "beta")$beta
  beta_post_fl <- rstan::extract(mfit_arp_fl, pars = "beta")$beta
  phi_post_nr <- rstan::extract(mfit_arp_nr, pars = "phi")$phi
  phi_post_fl <- rstan::extract(mfit_arp_fl, pars = "phi")$phi
  sigma_post_nr <- rstan::extract(mfit_arp_nr, pars = "sigma")$sigma
  sigma_post_r <- rstan::extract(mfit_arp_r, pars = "sigma")$sigma
  sigma_post_fl <- rstan::extract(mfit_arp_r, pars = "sigma")$sigma
  y_rep_nr <- rstan::extract(mfit_arp_nr, pars = "y_rep")$y_rep
  y_rep_r <- rstan::extract(mfit_arp_r, pars = "y_rep")$y_rep
  y_rep_fl <- rstan::extract(mfit_arp_fl, pars = "y_rep")$y_rep

  draws <- nrow(beta_post_nr)

  # matrix of draws from the posterior-predictive distribution
  # non-regularized model
  post_preds_nr <- matrix(nrow = draws, ncol = n)

  # regularized model
  post_preds_r <- matrix(nrow = draws, ncol = n)

  # flat model
  post_preds_fl <- matrix(nrow = draws, ncol = n)

  # fill in first p observations that are considered fixed
  post_preds_r[, 1:datlist$p] <- matrix(
    rep(y[1:datlist$p], each = draws), nrow = draws, ncol = datlist$p
  )
  post_preds_nr[, 1:datlist_nr$p] <- matrix(
    rep(y[1:datlist_nr$p], each = draws), nrow = draws, ncol = datlist_nr$p
  )
  post_preds_fl[, 1:datlist_fl$p] <- matrix(
    rep(y[1:datlist_nr$p], each = draws), nrow = draws, ncol = datlist_nr$p
  )

  # fill in post. pred. draws from stan
  post_preds_r[, (datlist$p + 1):(n - holdout)] <- y_rep_r
  post_preds_nr[, (datlist_nr$p + 1):(n - holdout)] <- y_rep_nr
  post_preds_fl[, (datlist_fl$p + 1):(n - holdout)] <- y_rep_fl

  for(i in 1:draws){
    for(t in (n - holdout + 1):n){
      # regularized model
      y_past_r <- as.double(post_preds_r[i, (t - datlist$p):(t - 1)])
      post_preds_r[i, t] <- alpha_post_r[i, 1] + X[t, -1] %*% beta_post_r[i, ] +
        phi_post_r[i, ] %*% y_past_r +
        rnorm(1, sd = sigma_post_r[i])

      # non-regularized model
      y_past_nr <- as.double(post_preds_nr[i, (t - datlist_nr$p):(t - 1)])
      post_preds_nr[i, t] <- X[t, ] %*% beta_post_nr[i, ] +
        phi_post_nr[i, ] %*% y_past_nr +
        rnorm(1, sd = sigma_post_nr[i])

      # flat model
      y_past_fl <- as.double(post_preds_fl[i, (t - datlist_fl$p):(t - 1)])
      post_preds_fl[i, t] <- X[t, ] %*% beta_post_fl[i, ] +
        phi_post_nr[i, ] %*% y_past_nr +
        rnorm(1, sd = sigma_post_fl[i])
    }
  }

  forecast_df <- data.frame(
    time = 1:n,
    y = as.double(y),
    estim_r = apply(post_preds_r, 2, mean),
    low_r = apply(post_preds_r, 2, quantile, probs = 0.025),
    high_r = apply(post_preds_r, 2, quantile, probs = 0.975),
    estim_nr = apply(post_preds_nr, 2, mean),
    low_nr = apply(post_preds_nr, 2, quantile, probs = 0.025),
    high_nr = apply(post_preds_nr, 2, quantile, probs = 0.975),
    estim_fl = apply(post_preds_fl, 2, mean),
    low_fl = apply(post_preds_fl, 2, quantile, probs = 0.025),
    high_fl = apply(post_preds_fl, 2, quantile, probs = 0.975)
  )

  # compute prediction root mean squared error for each model
  # RMSE_bayes() is a user-defined function in R/model_checking.R
  rmse_df <- data.frame(
    model = rep(c("Horseshoe", "Gaussian","Flat"), each = draws),
    rmse = c(
      RMSE_bayes(y[(n - holdout + 1):n], ppreds = post_preds_r[, (n - holdout + 1):n]),
      RMSE_bayes(y[(n - holdout + 1):n], ppreds = post_preds_nr[, (n - holdout + 1):n]),
      RMSE_bayes(y[(n - holdout + 1):n], ppreds = post_preds_fl[, (n - holdout + 1):n])
    )
  )

  # Plots comparing simulated values against post. pred intervals
  r_forecast <- ggplot(forecast_df[50:n, ], aes(x = time, y = y))+
    geom_ribbon(aes(ymin = low_r, ymax = high_r), fill = "brown", alpha = 0.5) +
    geom_point() +
    geom_line(linetype = "dashed") +
    geom_vline(xintercept = 100) +
    theme_classic() +
    ggtitle("Horseshoe priors")

  nr_forecast <- ggplot(forecast_df[50:n, ], aes(x = time, y = y)) +
    geom_ribbon(aes(ymin = low_nr, ymax = high_nr), fill = "brown", alpha = 0.5) +
    geom_point() +
    geom_line(linetype = "dashed") +
    geom_vline(xintercept = 100) +
    theme_classic() +
    ggtitle("Gaussian Priors")

  fl_forecast <- ggplot(forecast_df[50:n, ], aes(x = time, y = y)) +
    geom_ribbon(aes(ymin = low_fl, ymax = high_fl), fill = "brown", alpha = 0.5) +
    geom_point() +
    geom_line(linetype = "dashed") +
    geom_vline(xintercept = 100) +
    theme_classic() +
    ggtitle("Flat Priors")



  # posterior predictive distributions of RMSE
  rmse <- ggplot(rmse_df, aes(x = rmse, fill = model, color = model)) +
    geom_density(alpha = 0.5) +
    theme_classic() +
    scale_color_manual(values = c("black", "brown", "blue")) +
    scale_fill_manual(values = c("grey", "brown", "blue")) +
    xlab("Forecasting RMSE")+
    xlim(0,10)

# save the plot
  png(
    here("Figures/forecast_comparson_AR-p.png"),
    height = 3600, width = 1500,
    units = "px", res = 300
  )
    gridExtra::grid.arrange(rmse, nr_forecast, r_forecast, fl_forecast, ncol = 1)
  dev.off()









