
# Load libraries
library(tidyverse)
library(forecast)
library(patchwork)
devtools::load_all()

tree_dat <- readRDS(here::here("SparseTS_prismdata/ca719_BM_Seq_1800_2012.rds"))

# average all trees for a given year
tree_dat <- tree_dat %>%
  group_by(year) %>%
  summarise(
    mean_rwi = mean(rwi, na.rm = TRUE)
  )

# Create testing and training explanatory variables
train_yrs <- 1:which(tree_dat$year == 1970)
test_yrs <- (which(tree_dat$year == 1970) + 1):nrow(tree_dat)

# now create fourier terms
X_train <- forecast::fourier(
  ts(tree_dat$mean_rwi[train_yrs], frequency = length(train_yrs)),
  K = length(train_yrs) / 2
)

X_test <- continue_fourier(
  X_train,
  h = length(test_yrs),
  n = length(train_yrs)
)

# ---- Sparse model fitting ----

# compile stan model
sparse_mod <- rstan::stan_model(
  here::here("Stan/AR-p_err3_FHS_DG.stan")
)

# base_mod <- rstan::stan_model(
#   here::here("Stan/sparse_reg_FHS.stan")
# )

# create list of objects that are shared across all models
# get estimate of tau0
tau_0_phi <- tau0(
  y = tree_dat$mean_rwi[train_yrs],
  m0 = 1,
  M = 15,
  N = nrow(X_train),
  fam = "gaussian"
)

tau_0_beta <- tau0(
  y = tree_dat$mean_rwi[train_yrs],
  m0 = 5,
  M = ncol(X_train),
  N = nrow(X_train),
  fam = "gaussian"
)

dat_stan <- list(
  N = nrow(X_train),
  P_0 = 1,
  P = ncol(X_train) + 1,
  X = cbind(1, X_train),
  p = 10,
  tau0_phi = tau_0_phi,
  slab_scl_phi = 1,
  slab_df_phi = 50,
  tau0_beta = tau_0_beta,
  slab_scl_beta = 0.1,
  slab_df_beta = 6,
  y = tree_dat$mean_rwi[train_yrs],
  N_new = nrow(X_test),
  X_new = cbind(1, X_test)
)

# fit models
sparse_fit <- rstan::sampling(
  sparse_mod,
  data = dat_stan,
  chains = 4,
  iter = 4000,
  warmup = 2000,
  cores = 4,
  control = list(adapt_delta = 0.99, max_treedepth = 15)
)

# extract posterior predictive draws
y_hat <- rstan::extract(sparse_fit, pars = "y_rep")$y_rep

plot_df_allyrs <- data.frame(
  y = c(tree_dat$mean_rwi, colMeans(y_hat)),
  low = c(
    rep(NA, nrow(tree_dat)),
    apply(y_hat, 2, quantile, probs = 0.025)
  ),
  high = c(
    rep(NA, nrow(tree_dat)),
    apply(y_hat, 2, quantile, probs = 0.975)
  ),
  year = rep(tree_dat$year, 2),
  source = rep(c("Observed", "RHS"), each = nrow(tree_dat))
)

plot_df_allyrs$y[plot_df_allyrs$source == "RHS" & plot_df_allyrs$year <= 1970] <- NA
plot_df_allyrs$low[1:length(train_yrs)] <- NA
plot_df_allyrs$high[1:length(train_yrs)] <- NA

# Create forecasting plot function
forecast_plot <- function(df, horizon, col){
  ggplot(df, aes(x = year, y = y, colour = source, linewidth = source, fill = source)) +
    geom_line() +
    geom_ribbon(aes(ymin = low, ymax = high), alpha = 0.2)
}

all_dat_RHS <- forecast_plot(plot_df_allyrs, horizon = 1970, col = "#33406fff")

# now fit auto.arima model

auto_arima_fit <- fit_seasonal_arima_model(
  model_pars = list(
    n = nrow(X_train),
    p = 10,
    X = matrix(1, nrow = nrow(X_train), ncol = 1),
    y = tree_dat$mean_rwi,
    holdout = length(test_yrs)
  ),
  freq = length(train_yrs)
)

# create dataframe for plotting
fc_aarima <- forecast(
  auto_arima_fit$fit_ar,
  h = length(test_yrs),
  xreg = X_test[, colnames(X_test) %in% auto_arima_fit$terms]
) |> as.data.frame() |> janitor::clean_names()

plot_df_arima <- data.frame(
  y = tree_dat$mean_rwi,
  y_hat = c(
    rep(NA, nrow(X_train)),
    fc_aarima$point_forecast
  ),
  low = c(rep(NA, nrow(X_train)), fc_aarima$lo_95),
  high = c(rep(NA, nrow(X_train)), fc_aarima$hi_95),
  year = tree_dat$year
)

all_dat_aarima <- forecast_plot(plot_df_arima, horizon = 1970, p = 10, col = "#818181")




# ---- Analysis for 1800 - 1925 ----

# Create testing and training explanatory variables
yrs2 <- which(tree_dat$year %in% 1800:1925)
train_yrs2 <- which(tree_dat$year %in% 1800:1900)
test_yrs2 <- which(tree_dat$year %in% 1901:1925)

X_train2 <- forecast::fourier(
  ts(tree_dat$mean_rwi[train_yrs2], frequency = length(train_yrs2)),
  K = length(train_yrs2) / 2
)

X_test2 <- continue_fourier(
  X_train2,
  h = length(test_yrs2),
  n = length(train_yrs2)
)

# create data for stan model
dat_stan2 <- modifyList(
  dat_stan,
  list(
    N = nrow(X_train2),
    P = ncol(X_train2) + 1,
    X = cbind(1, X_train2),
    y = tree_dat$mean_rwi[train_yrs2],
    N_new = nrow(X_test2),
    X_new = cbind(1, X_test2),
    tau0_phi = tau0(
      y = tree_dat$mean_rwi[train_yrs2],
      m0 = 1,
      M = 100,
      N = nrow(X_train2),
      fam = "gaussian"
    ),
    tau0_beta = tau0(
      y = tree_dat$mean_rwi[train_yrs2],
      m0 = 5,
      M = ncol(X_train2),
      N = nrow(X_train2),
      fam = "gaussian"
    )
  )
)

# fit the RHS model
sparse_fit2 <- rstan::sampling(
  sparse_mod,
  data = dat_stan2,
  chains = 4,
  iter = 4000,
  warmup = 3000,
  cores = 4,
  control = list(adapt_delta = 0.99, max_treedepth = 15)
)

# extract posterior predictive draws
y_hat2 <- rstan::extract(sparse_fit2, pars = "y_rep")$y_rep

plot_df_RHS2 <- data.frame(
  y = tree_dat$mean_rwi,
  y_hat = c(
    colMeans(y_hat2),
    rep(NA, nrow(tree_dat) - length(yrs2))
  ),
  low = c(
    apply(y_hat2, 2, quantile, probs = 0.025),
    rep(NA, nrow(tree_dat) - length(yrs2))
  ),
  high = c(
    apply(y_hat2, 2, quantile, probs = 0.975),
    rep(NA, nrow(tree_dat) - length(yrs2))
  ),
  year = tree_dat$year
)

plot_df_RHS2$y_hat[train_yrs2] <- NA
plot_df_RHS2$low[train_yrs2] <- NA
plot_df_RHS2$high[train_yrs2] <- NA
plot_df_RHS2$y[(yrs2[length(yrs2)] + 1):nrow(tree_dat)] <- NA

dat2_plot_RHS <- forecast_plot(plot_df_RHS2, horizon = 1900, p = 10, col = "#33406fff")

# repeat the process with auto.arima

auto_arima_fit2 <- fit_seasonal_arima_model(
  model_pars = list(
    n = nrow(X_train2),
    p = 10,
    X = matrix(1, nrow = nrow(X_train2), ncol = 1),
    y = tree_dat$mean_rwi[train_yrs2],
    holdout = length(test_yrs2)
  ),
  freq = length(train_yrs2)
)

# create dataframe for plotting
fc_aarima2 <- forecast(
  auto_arima_fit2$fit_ar,
  h = length(test_yrs2),
  xreg = X_test2[, colnames(X_test2) %in% auto_arima_fit2$terms]
) |> as.data.frame() |> janitor::clean_names()

plot_df_arima2 <- data.frame(
  y = tree_dat$mean_rwi,
  y_hat = c(
    rep(NA, nrow(X_train2)),
    fc_aarima2$point_forecast,
    rep(NA, nrow(tree_dat) - length(yrs2))
  ),
  low = c(
    rep(NA, nrow(X_train2)),
    fc_aarima2$lo_95,
    rep(NA, nrow(tree_dat) - length(yrs2))
  ),
  high = c(
    rep(NA, nrow(X_train2)),
    fc_aarima2$hi_95,
    rep(NA, nrow(tree_dat) - length(yrs2))
  ),
  year = tree_dat$year
)

# omit the later time points
plot_df_arima2$y[(yrs2[length(yrs2)] + 1):nrow(tree_dat)] <- NA

dat2_plot_aarima <- forecast_plot(plot_df_arima2, horizon = 1900, p = 10, col = "#818181")


# ---- Analysis for 1926 - 2012 ----

# Create testing and training explanatory variables
yrs3 <- which(tree_dat$year %in% 1926:2012)
train_yrs3 <- which(tree_dat$year %in% 1926:1995)
test_yrs3 <- which(tree_dat$year %in% 1996:2012)

X_train3 <- forecast::fourier(
  ts(tree_dat$mean_rwi[train_yrs3], frequency = length(train_yrs3)),
  K = length(train_yrs3) / 2
)

X_test3 <- continue_fourier(
  X_train3,
  h = length(test_yrs3),
  n = length(train_yrs3)
)

# create data for stan model

dat_stan3 <- modifyList(
  dat_stan,
  list(
    N = nrow(X_train3),
    P = ncol(X_train3) + 1,
    X = cbind(1, X_train3),
    y = tree_dat$mean_rwi[train_yrs3],
    N_new = nrow(X_test3),
    X_new = cbind(1, X_test3),
    tau0_phi = tau0(
      y = tree_dat$mean_rwi[train_yrs3],
      m0 = 1,
      M = 100,
      N = nrow(X_train3),
      fam = "gaussian"
    ),
    tau0_beta = tau0(
      y = tree_dat$mean_rwi[train_yrs3],
      m0 = 5,
      M = ncol(X_train3),
      N = nrow(X_train3),
      fam = "gaussian"
    )
  )
)

# fit the RHS model
sparse_fit3 <- rstan::sampling(
  sparse_mod,
  data = dat_stan3,
  iter = 4000,
  warmup = 2000,
  cores = 4,
  control = list(adapt_delta = 0.99, max_treedepth = 15)
)

# extract posterior predictive draws
y_hat3 <- rstan::extract(sparse_fit3, pars = "y_rep")$y_rep

plot_df_RHS3 <- data.frame(
  y = tree_dat$mean_rwi,
  y_hat = c(
    rep(NA, nrow(tree_dat) - length(yrs3)),
    colMeans(y_hat3)
  ),
  low = c(
    rep(NA, nrow(tree_dat) - length(yrs3)),
    apply(y_hat3, 2, quantile, probs = 0.025)
  ),
  high = c(
    rep(NA, nrow(tree_dat) - length(yrs3)),
    apply(y_hat3, 2, quantile, probs = 0.975)
  ),
  year = tree_dat$year
)

plot_df_RHS3$y[1:(yrs3[1] - 1)] <- NA
plot_df_RHS3$y_hat[train_yrs3] <- NA
plot_df_RHS3$low[train_yrs3] <- NA
plot_df_RHS3$high[train_yrs3] <- NA

dat3_plot_RHS <- forecast_plot(plot_df_RHS3, horizon = 1996, p = 10, col = "#33406fff")


# repeat the process with auto.arima

auto_arima_fit3 <- fit_seasonal_arima_model(
  model_pars = list(
    n = nrow(X_train3),
    p = 10,
    X = matrix(1, nrow = nrow(X_train3), ncol = 1),
    y = tree_dat$mean_rwi[train_yrs3],
    holdout = length(test_yrs3)
  ),
  freq = length(train_yrs3)
)

fc_aarima3 <- forecast(
  auto_arima_fit3$fit_ar,
  h = length(test_yrs3),
  xreg = X_test3[, colnames(X_test3) %in% auto_arima_fit3$terms]
) |> as.data.frame() |> janitor::clean_names()


plot_df_arima3 <- data.frame(
  y = tree_dat$mean_rwi,
  y_hat = c(
    rep(NA, nrow(tree_dat) - length(test_yrs3)),
    fc_aarima3$point_forecast
  ),
  low = c(
    rep(NA, nrow(tree_dat) - length(test_yrs3)),
    fc_aarima3$lo_95
  ),
  high = c(
    rep(NA, nrow(tree_dat) - length(test_yrs3)),
    fc_aarima3$hi_95
  ),
  year = tree_dat$year
)

plot_df_arima3$y[1:(yrs3[1] - 1)] <- NA

dat3_plot_aarima <- forecast_plot(plot_df_arima3, horizon = 1996, p = 10, col = "#818181")


# ---- Combine all plots ----

# first, remove y axis labels from arima plots
all_dat_aarima <- all_dat_aarima +
  ylab("")

dat2_plot_aarima <- dat2_plot_aarima +
  ylab("")

dat3_plot_aarima <- dat3_plot_aarima +
  ylab("")

# remove x axis labels from all but dat 3 plots
all_dat_RHS <- all_dat_RHS +
  theme(axis.text.x = element_blank()) +
  xlab("")
all_dat_aarima <- all_dat_aarima +
  theme(axis.text.x = element_blank()) +
  xlab("")
dat2_plot_RHS <- dat2_plot_RHS +
  theme(axis.text.x = element_blank()) +
  xlab("")
dat2_plot_aarima <- dat2_plot_aarima +
  theme(axis.text.x = element_blank()) +
  xlab("")

# angle x axis text for dat 3 plots
dat3_plot_RHS <- dat3_plot_RHS +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("Year")
dat3_plot_aarima <- dat3_plot_aarima +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("Year")

# fix y axis labels
all_dat_RHS <- all_dat_RHS +
  ylab("RWI")
dat2_plot_RHS <- dat2_plot_RHS +
  ylab("RWI")
dat3_plot_RHS <- dat3_plot_RHS +
  ylab("RWI")

# add titles to top plots
all_dat_RHS <- all_dat_RHS +
  ggtitle("RHS model")

all_dat_aarima <- all_dat_aarima +
  ggtitle("Auto-ARIMA model")

(all_dat_RHS + all_dat_aarima) /
  (dat2_plot_RHS + dat2_plot_aarima) /
  (dat3_plot_RHS + dat3_plot_aarima)

ggsave(
  here::here("Figures/tree_forecasts_fourier.png"),
  width = 6,
  height = 5
)
