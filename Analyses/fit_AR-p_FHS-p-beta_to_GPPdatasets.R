#############################################################
# This script runs an AR-p model on several actual GPP time
# series that searches for important explanatory variables
# and important time lags in both the GPP series and in one
# covariate, intended to represent discharge. The purpose
# is to get a sense of what parameter space we should be
# exploring in our simulated datasets.
#############################################################

# libraries
library(rstan)
library(here)
library(ggplot2)
devtools::load_all()

dd <- read_csv('Data/aquatic_sim_data/NWIS_data.csv')
sites <- unique(dd$site_name)
#### Load in data sets ####
dat <- filter(dd, site_name == sites[4]) %>%
  mutate(logQ = log(Q),
         light = scale(light)[,1],
         Q = scale(logQ)[,1],
         temp = scale(temp)[,1]) %>%
  select(date, GPP, temp, light, Q)


n_lag_Q <- 5 # number of potential lags in discharge
Q <- data.frame(matrix(rep(NA, nrow(dat) * n_lag_Q),
                       nrow = nrow(dat), ncol = n_lag_Q))

for(i in 1:n_lag_Q){
  Q[,i] <- c(rep(NA, i), dat$Q[1:(nrow(Q)-i)])
}

# predictor matrix:
X = mutate(dat, intercept = 1) %>%
  select(intercept, temp, light, Q) %>%
  bind_cols(Q) %>%
  slice(-c(1:n_lag_Q)) %>%
  as.matrix()

y = dat$GPP[(n_lag_Q+1):nrow(dat)]

# length of time series
n <- nrow(X)

#### Fit the model ####

  # compile stan model
  arp_r <- stan_model(here("Stan/AR-p_FHS-p-beta.stan"))

  # compute prior guess for tau0 based on a guess of 5
  #  non-zero coefficients
  #  see ?tau0() for documentation
  tau_0 <- tau0(
    y = y,
    m0 = 5,
    M = ncol(X) -1 + 20,
    N = n,
    fam = "gaussian"
  )

  # compile data (see Stan file for descriptions of each input)
  datlist <- list(
    N = n,
    P0 = 1,
    P = ncol(X)-1,
    p = 20,
    y = y,
    X_alpha = matrix(X[, 1], ncol = 1),
    X_beta = X[, -1],
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
    param = paste0("beta", 1:(ncol(X) - 1)),
    beta_hat = apply(beta_post_r, 2, mean),
    low = apply(beta_post_r, 2, quantile, probs = 0.025),
    high = apply(beta_post_r, 2, quantile, probs = 0.975)
  )

  df_phi <- data.frame(
    param = paste0("phi", 1:ncol(phi_post_r)),
    # note that the ordering of phi is reversed in the stan model
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
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ylab("value") +
    ggtitle("Covariate coefficients")

  plot_phi <- ggplot(data = df_phi, aes(x = param)) +
    geom_errorbar(aes(ymin = low, ymax = high), width = 0) +
    geom_point(aes(y = phi_hat)) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ylab("value") +
    ggtitle("AR coefficients")

  ## save figure ##
  lo <- rbind(
    c(1, 1, 1, 1),
    c(2, 2, 3, 3)
  )
  gridExtra::grid.arrange(plot_ts, plot_phi, plot_beta, layout_matrix = lo)

