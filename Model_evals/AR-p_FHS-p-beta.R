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
library(dplyr)
devtools::load_all("/Users/amypatterson/Documents/Laramie_postdoc/Sparse_time_series_project/sponges")

################################################################
#### Simulating the data ####

#set.seed(198440)# maybe
set.seed(198441) # maybe
# length of time series
n <- 150

# time series parameters, constrained for stationarity
#phi <- c(0.6, rep(0, 4), -0.2, rep(0, 9), 0.3)
phi <- 0

# coefficient vector with two non-zero elements and the intercepts
beta <- c(0.5, 1, 2, runif(9, min=0, max=0.1), 1.5, runif(7, min=0, max=0.1),1.8,runif(13,min=0,max=0.1))

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


################################################################
####  Perform the 3 fits ####
#### Fit the model, holding out last few observations for forecasting ####

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
  M = P + length(phi),
  N = n - holdout,
  fam = "gaussian"
)

# compile data (see Stan file for descriptions of each input)
datlist <- list(
  N = n - holdout,
  P0 = 1,
  P = P,
  p = length(phi),
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
  chains = 3,
  cores = 3,
  iter = 4000
)


#### Compare forecasting to a non-regularized fitted model #####

# fitting the non-regularized AR model Gaussian priors
datlist_nr <- list(
  N = n - holdout,
  P = P + 1,
  p = length(phi),
  y = y[1:(n - holdout)],
  X = X[1:(n - holdout), ]
)

arp_nr <- stan_model(here("Stan/AR-p.stan"))

mfit_arp_nr <- sampling(
  arp_nr,
  data = datlist_nr,
  chains = 3,
  cores = 3,
  iter = 4000,
  control = list(max_treedepth = 15)
)


#### Compare forecasting to a model with flat priors (no priors) #####

# fitting the non-regularized AR model with flat priors
datlist_fl <- list(
  N = n - holdout,
  P = P + 1,
  p = length(phi),
  y = y[1:(n - holdout)],
  X = X[1:(n - holdout), ]
)

arp_fl <- stan_model(here("Stan/AR-p_flat.stan"))

mfit_arp_fl <- sampling(
  arp_fl,
  data = datlist_fl,
  chains = 3,
  cores = 3,
  iter = 4000,
  control = list(max_treedepth = 15)
)



################################################################
#### Get parameter estimates for all models ####
#### Compare parameter posteriors to truth ####
alpha_post_r <- rstan::extract(mfit_arp_r, pars = "alpha")$alpha # this is the intercept
beta_post_r <- rstan::extract(mfit_arp_r, pars = "beta")$beta
phi_post_r <- rstan::extract(mfit_arp_r, pars = "phi")$phi

beta_post_nr <- rstan::extract(mfit_arp_nr, pars = "beta")$beta
phi_post_nr <- rstan::extract(mfit_arp_nr, pars = "phi")$phi

beta_post_fl <- rstan::extract(mfit_arp_fl, pars = "beta")$beta
phi_post_fl <- rstan::extract(mfit_arp_fl, pars = "phi")$phi


# create data frames for plotting
df_beta <- data.frame(
  param = paste0("beta", 1:(length(beta) - 1)),
  beta_true = beta[-1],
  beta_hat_r = apply(beta_post_r, 2, mean),
  low_r = apply(beta_post_r, 2, quantile, probs = 0.025),
  high_r = apply(beta_post_r, 2, quantile, probs = 0.975),
  beta_hat_nr = apply(beta_post_nr, 2, mean)[2:length(beta)],
  low_nr = apply(beta_post_nr, 2, quantile, probs = 0.025)[2:length(beta)],
  high_nr = apply(beta_post_nr, 2, quantile, probs = 0.975)[2:length(beta)],
  beta_hat_fl = apply(beta_post_fl, 2, mean)[2:length(beta)],
  low_fl = apply(beta_post_fl, 2, quantile, probs = 0.025)[2:length(beta)],
  high_fl = apply(beta_post_fl, 2, quantile, probs = 0.975)[2:length(beta)]
)

df_phi <- data.frame(
  param = paste0("phi", 1:ncol(phi_post_r)),
  # note that the ordering of phi is reversed in the stan model
  phi_true = c(phi, rep(0, 4))[ncol(phi_post_r):1],
  phi_hat_r = apply(phi_post_r, 2, mean),
  low_r = apply(phi_post_r, 2, quantile, probs = 0.025),
  high_r = apply(phi_post_r, 2, quantile, probs = 0.975),
  phi_hat_nr = apply(phi_post_nr, 2, mean),
  low_nr = apply(phi_post_nr, 2, quantile, probs = 0.025),
  high_nr = apply(phi_post_nr, 2, quantile, probs = 0.975),
  phi_hat_fl = apply(phi_post_fl, 2, mean),
  low_fl = apply(phi_post_fl, 2, quantile, probs = 0.025),
  high_fl = apply(phi_post_fl, 2, quantile, probs = 0.975)
)

################################################################
# plots for posterior estimates and CIs compared to simulated values
plot_beta <- ggplot(data = df_beta, aes(x = 1:nrow(df_beta))) +
  geom_rect(aes(xmin=20.5,xmax=21.5,ymin=-1,ymax=3),fill="gray") +
  geom_rect(aes(xmin=11.5,xmax=12.5,ymin=-1,ymax=3),fill="white",linewidth=1,color = "black") +
  #geom_rect(aes(xmin=0.5,xmax=1.5,ymin=-1,ymax=3),fill="gray") +
  geom_errorbar(aes(ymin = low_fl, ymax = high_fl, x=(1:nrow(df_beta))-0.2), width = 0) +
  geom_point(aes(y = beta_hat_fl, x=(1:nrow(df_beta))-0.2)) +
  geom_errorbar(aes(ymin = low_nr, ymax = high_nr), width = 0,color="#a52a2aff") +
  geom_point(aes(y = beta_hat_nr),color = "#a52a2aff") +
  geom_errorbar(aes(ymin = low_r, ymax = high_r, x=(1:nrow(df_beta))+0.2), width = 0,color="#33406fff") +
  geom_point(aes(y = beta_hat_r, x=(1:nrow(df_beta))+0.2), color = "#33406fff") +
  geom_point(aes(y = beta_true), shape = "l", color = "green3", size = 5) +
  theme_classic() +
  theme(axis.title.y = element_blank(),axis.text.y = element_blank(),axis.ticks.y = element_blank()) +
  annotate("point", x = 9, y = 1.8, shape = "l", color = "green3", size = 5) +
  annotate("text", label = "Actual", x = 9, y = 2, hjust = 0, size = 3) +
  annotate("point", x = 8, y = 1.8, color = "#33406fff") +
  annotate("text", x = 8, y = 2, label = "Horseshoe", hjust = 0, size = 3) +
  annotate("point", x = 7, y = 1.8, color = "#a52a2aff") +
  annotate("text", x = 7, y = 2, label = "Gaussian", hjust = 0, size = 3) +
  annotate("point", x = 6, y = 1.8) +
  annotate("text", x = 6, y = 2, label = "Flat", hjust = 0, size = 3) +
  ylab("value") +
  ggtitle("a) Covariate coefficients")

plot_phi <- ggplot(data = df_phi, aes(x = 1:nrow(df_phi))) +
  geom_errorbar(aes(ymin = low_fl, ymax = high_fl, x=(1:nrow(df_phi))-0.2), width = 0,color="#33406fff") +
  geom_point(aes(y = phi_hat_fl, x=(1:nrow(df_phi))-0.2),color="#33406fff") +
  geom_errorbar(aes(ymin = low_nr, ymax = high_nr), width = 0,color="#a52a2aff") +
  geom_point(aes(y = phi_hat_nr),color="#a52a2aff") +
  geom_errorbar(aes(ymin = low_r, ymax = high_r, x=(1:nrow(df_phi))+0.2), width = 0) +
  geom_point(aes(y = phi_hat_r), x=(1:nrow(df_phi))+0.2) +
  geom_point(aes(y = phi_true), shape = "|", color = "goldenrod2", size = 4) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 0, hjust = 1)) +
  annotate("point", x = 6, y = 0.35, shape = "|", color = "goldenrod2", size = 4) +
  annotate("text", label = "Actual", x = 6, y = 0.4, hjust = 0, size = 3) +
  annotate("point", x = 5, y = 0.35) +
  annotate("text", x = 5, y = 0.4, label = "Horseshoe", hjust = 0, size = 3) +
  annotate("point", x = 4, y = 0.35, color = "#a52a2aff") +
  annotate("text", x = 4, y = 0.4, label = "Gaussian", hjust = 0, size = 3) +
  annotate("point", x = 3, y = 0.35, color = "#33406fff") +
  annotate("text", x = 3, y = 0.4, label = "Flat", hjust = 0, size = 3) +
  ylab("value") +
  xlab("parameters") +
  ggtitle("AR coefficients")

#names(df_beta)=names(df_phi)
# df_betaphi=rbind(df_phi,df_beta)
#
# plot_betaphi <- ggplot(data = df_betaphi, aes(x = 1:nrow(df_betaphi))) +
#   geom_rect(aes(xmin=11.5,xmax=12.5,ymin=-1,ymax=3),fill="gray") +
#   geom_rect(aes(xmin=18.5,xmax=19.5,ymin=-1,ymax=3),fill="gray") +
#   geom_errorbar(aes(ymin = low_fl, ymax = high_fl, x=(1:nrow(df_betaphi))-0.2), width = 0) +
#   geom_point(aes(y = phi_hat_fl, x=(1:nrow(df_betaphi))-0.2)) +
#   geom_errorbar(aes(ymin = low_nr, ymax = high_nr), width = 0,color="#a52a2aff") +
#   geom_point(aes(y = phi_hat_nr),color="#a52a2aff") +
#   geom_errorbar(aes(ymin = low_r, ymax = high_r, x=(1:nrow(df_betaphi))+0.2), width = 0, color="#33406fff") +
#   geom_point(aes(y = phi_hat_r), x=(1:nrow(df_betaphi))+0.2, color="#33406fff") +
#   geom_point(aes(y = phi_true), shape = "|", color = "goldenrod3", size = 4) +
#   theme_classic() +
#   theme(axis.text.y = element_blank(),axis.ticks.y = element_blank(),axis.title.y=element_blank()) +
#   annotate("point", x = 6, y = 1.65, shape = "|", color = "goldenrod3", size = 4) +
#   annotate("text", label = "Actual", x = 6, y = 1.75, hjust = 0, size = 3) +
#   annotate("point", x = 5, y = 1.65, color = "#33406fff") +
#   annotate("text", x = 5, y = 1.75, label = "Horseshoe", hjust = 0, size = 3) +
#   annotate("point", x = 4, y = 1.65, color = "#a52a2aff") +
#   annotate("text", x = 4, y = 1.75, label = "Gaussian", hjust = 0, size = 3) +
#   annotate("point", x = 3, y = 1.65) +
#   annotate("text", x = 3, y = 1.75, label = "Flat", hjust = 0, size = 3) +
#   ylab("value") +
#   xlab("parameters") +
#   ggtitle("Parameter Estimates")
#


################################################################
#### Get predictions for all models ####
# COMPARISON STARTS

# forecast the held-out observations
beta_post_nr <- rstan::extract(mfit_arp_nr, pars = "beta")$beta
beta_post_fl <- rstan::extract(mfit_arp_fl, pars = "beta")$beta
beta_post_r <- rstan::extract(mfit_arp_r, pars = "beta")$beta
phi_post_nr <- rstan::extract(mfit_arp_nr, pars = "phi")$phi
phi_post_fl <- rstan::extract(mfit_arp_fl, pars = "phi")$phi
phi_post_r <- rstan::extract(mfit_arp_r, pars = "phi")$phi
sigma_post_nr <- rstan::extract(mfit_arp_nr, pars = "sigma")$sigma
sigma_post_r <- rstan::extract(mfit_arp_r, pars = "sigma")$sigma
sigma_post_fl <- rstan::extract(mfit_arp_fl, pars = "sigma")$sigma
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
  rep(y[1:datlist_fl$p], each = draws), nrow = draws, ncol = datlist_fl$p
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
      phi_post_fl[i, ] %*% y_past_fl +
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

# n_xAxis is the maximum value of the x-axis (i.e. how many predicted observations are included in the figure)
n_xAxis <- 115

################################################################
#### Plot predictions all models ####

# Show x after time...
xmin=80

# Plots comparing simulated values against post. pred intervals
# variables to determine spread of y-axis (min and max observed +/- rmse for flat prior model)
y_minVal <- min(forecast_df[xmin:n_xAxis, "low_fl"])-1
y_maxVal <- max(forecast_df[xmin:n_xAxis, "high_fl"])+1

r_forecast <- ggplot(forecast_df[xmin:n_xAxis, ], aes(x = time, y = y))+
  geom_ribbon(aes(ymin = low_r, ymax = high_r), fill = "#33406fff", alpha = 0.5) +
  #geom_point() +
  geom_line(linetype = "solid") +
  geom_vline(xintercept = 100) +
  annotate("rect",xmin=80,xmax=99.9,ymin=y_minVal,ymax=y_maxVal,alpha=0.6,fill="white") +
  ylim(c(y_minVal, y_maxVal)) +
  xlim(c(xmin,n_xAxis)) +
  theme_classic() +
  ggtitle("b) Horseshoe priors") +
  geom_segment(x=n_xAxis+0.5,y=(forecast_df$low_r[n_xAxis]+forecast_df$high_r[n_xAxis])/2+1.5,xend=n_xAxis+0.5,yend=forecast_df$high_r[n_xAxis],arrow=arrow(length=unit(0.03,"npc"),ends="last"))+
  geom_segment(x=n_xAxis+0.5,y=(forecast_df$low_r[n_xAxis]+forecast_df$high_r[n_xAxis])/2-1.5,xend=n_xAxis+0.5,yend=forecast_df$low_r[n_xAxis],arrow=arrow(length=unit(0.03,"npc"),ends="last"))+
  geom_text(x=n_xAxis+1.3,y=(forecast_df$low_r[n_xAxis]+forecast_df$high_r[n_xAxis])/2,label=round(forecast_df$high_r[n_xAxis]-forecast_df$low_r[n_xAxis],digits=1),size=4)+
  coord_cartesian(clip = "off")+
  scale_x_continuous(expand = expansion(mult = c(0, 0.09)))



nr_forecast <- ggplot(forecast_df[xmin:n_xAxis, ], aes(x = time, y = y)) +
  geom_ribbon(aes(ymin = low_nr, ymax = high_nr), fill = "#a52a2aff", alpha = 0.5) +
  #geom_point() +
  geom_line(linetype = "solid") +
  geom_vline(xintercept = 100) +
  annotate("rect",xmin=80,xmax=99.9,ymin=y_minVal,ymax=y_maxVal,alpha=0.6,fill="white") +
  annotate("point",x=111.5,y=0.3, pch=21,size=10)+
  ylim(c(y_minVal, y_maxVal)) +
  xlim(c(xmin,n_xAxis)) +
  theme_classic() +
  ggtitle("c) Gaussian Priors")+
  geom_segment(x=n_xAxis+0.5,y=(forecast_df$low_nr[n_xAxis]+forecast_df$high_nr[n_xAxis])/2+1.5,xend=n_xAxis+0.5,yend=forecast_df$high_nr[n_xAxis],arrow=arrow(length=unit(0.03,"npc"),ends="last"))+
  geom_segment(x=n_xAxis+0.5,y=(forecast_df$low_nr[n_xAxis]+forecast_df$high_nr[n_xAxis])/2-1.5,xend=n_xAxis+0.5,yend=forecast_df$low_nr[n_xAxis],arrow=arrow(length=unit(0.03,"npc"),ends="last"))+
  geom_text(x=n_xAxis+1.4,y=(forecast_df$low_nr[n_xAxis]+forecast_df$high_nr[n_xAxis])/2,label=round(forecast_df$high_nr[n_xAxis]-forecast_df$low_nr[n_xAxis],digits=1),size=4)+
  coord_cartesian(clip = "off")+
  scale_x_continuous(expand = expansion(mult = c(0, 0.09)))


fl_forecast <- ggplot(forecast_df[xmin:n_xAxis, ], aes(x = time, y = y)) +
  geom_ribbon(aes(ymin = low_fl, ymax = high_fl), fill = "#bebebeff", alpha = 0.5) +
  #geom_point() +
  geom_line(linetype = "solid") +
  geom_vline(xintercept = 100) +
  annotate("rect",xmin=80,xmax=99.9,ymin=y_minVal,ymax=y_maxVal,alpha=0.6,fill="white") +
  annotate("point",x=111.5,y=0.3, pch=21,size=10)+
  ylim(c(y_minVal, y_maxVal)) +
  xlim(c(xmin,n_xAxis)) +
  theme_classic() +
  ggtitle("d) Flat Priors")+
  geom_segment(x=n_xAxis+0.5,y=(forecast_df$low_fl[n_xAxis]+forecast_df$high_fl[n_xAxis])/2+1.5,xend=n_xAxis+0.5,yend=forecast_df$high_fl[n_xAxis],arrow=arrow(length=unit(0.03,"npc"),ends="last"))+
  geom_segment(x=n_xAxis+0.5,y=(forecast_df$low_fl[n_xAxis]+forecast_df$high_fl[n_xAxis])/2-1.5,xend=n_xAxis+0.5,yend=forecast_df$low_fl[n_xAxis],arrow=arrow(length=unit(0.03,"npc"),ends="last"))+
  geom_text(x=n_xAxis+1.4,y=(forecast_df$low_fl[n_xAxis]+forecast_df$high_fl[n_xAxis])/2,label=round(forecast_df$high_fl[n_xAxis]-forecast_df$low_fl[n_xAxis],digits=1),size=4)+
  coord_cartesian(clip = "off")+
  scale_x_continuous(expand = expansion(mult = c(0, 0.09)))


################################################################
#### Plot full composite figure ####

lo <- rbind(
  c(4,1),
  c(4,2),
  c(4,3)
)

png(
  filename = here("Figures/newFig1.png"),
  width = 2700,
  height = 2100,
  res = 300, units = "px"
)

gridExtra::grid.arrange(r_forecast, nr_forecast, fl_forecast,plot_beta+coord_flip(),layout_matrix = lo)

dev.off()




