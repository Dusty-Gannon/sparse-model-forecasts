
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
train_yrs <- 1800:1970
test_yrs <- 1971:tree_dat$year[nrow(tree_dat)]

# now create fourier terms
X_train <- forecast::fourier(
  ts(
    tree_dat$mean_rwi[tree_dat$year %in% train_yrs],
    frequency = length(train_yrs)
  ),
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
  y = tree_dat$mean_rwi[tree_dat$year %in% train_yrs],
  m0 = 1,
  M = 15,
  N = nrow(X_train),
  fam = "gaussian"
)

tau_0_beta <- tau0(
  y = tree_dat$mean_rwi[tree_dat$year %in% train_yrs],
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
  y = tree_dat$mean_rwi[tree_dat$year %in% train_yrs],
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
plot_df_allyrs$low[plot_df_allyrs$source == "RHS" & plot_df_allyrs$year <= 1970] <- NA
plot_df_allyrs$high[plot_df_allyrs$source == "RHS" & plot_df_allyrs$year <= 1970] <- NA

# Create forecasting plot function
forecast_plot <- function(df, horizon, col){
  ggplot(
    df,
    aes(x = year, y = y, colour = source, linewidth = source, fill = source)
  ) +
    geom_line() +
    geom_ribbon(
      aes(ymin = low, ymax = high, alpha = source),
      linewidth = NA
    ) +
    geom_vline(xintercept = horizon, linetype = "dashed", color = "brown") +
    geom_vline(xintercept = 1926, linetype = "dashed", color = "grey") +
    scale_color_manual(
      values = c("Observed" = "black", "RHS" = col[1], "Auto-ARIMA" = col[2])
    ) +
    scale_fill_manual(
      values = c("Observed" = "black", "RHS" = col[1], "Auto-ARIMA" = col[2])
    ) +
    scale_linewidth_manual(
      values = c("Observed" = 0.5, "RHS" = 0.8, "Auto-ARIMA" = 0.8)
    ) +
    scale_alpha_manual(
      values = c("Observed" = 1, "RHS" = 0.2, "Auto-ARIMA" = 0.6)
    ) +
    theme_classic() +
    theme(legend.title = element_blank()) +
    ylab("RWI")
}

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

# combine with data for plotting

plot_df_allyrs <- data.frame(
  y = c(
    rep(NA, nrow(X_train)),
    fc_aarima$point_forecast
  ),
  low = c(rep(NA, nrow(X_train)), fc_aarima$lo_95),
  high = c(rep(NA, nrow(X_train)), fc_aarima$hi_95),
  year = tree_dat$year,
  source = rep("Auto-ARIMA", nrow(tree_dat))
) %>% rbind(plot_df_allyrs, .)

### ---- First ts plot ----
all_dat_plot <- forecast_plot(
  plot_df_allyrs,
  horizon = 1970,
  col = PNWColors::pnw_palette("Sunset", 7)[c(2,5)]
) +
  xlab("") +
  ylim(c(0.5, 3)) +
  theme(
    legend.position = "top",
    axis.text.x = element_blank()
  ) +
  ggtitle("a)")


## ---- Get RMSE of the forecast ----

rmse_all_aa <- sqrt(
  mean(
    (tree_dat$mean_rwi[tree_dat$year %in% test_yrs] -
       fc_aarima$point_forecast)^2,
    na.rm = TRUE
  )
)

rmse_all_RHS <- data.frame(
  rmse = RMSE_bayes(
    obs = tree_dat$mean_rwi[tree_dat$year %in% test_yrs],
    ppreds = y_hat[, which(tree_dat$year %in% test_yrs)]
  )
)
### ---- make figure for RMSE ----

postprob <- mean(rmse_all_RHS$rmse < rmse_all_aa)

rmse_all_plot <- ggplot(rmse_all_RHS, aes(x = rmse)) +
  geom_density(fill = PNWColors::pnw_palette("Sunset", 7)[2]) +
  geom_vline(
    xintercept = rmse_all_aa,
    color = PNWColors::pnw_palette("Sunset", 7)[5]) +
  theme_classic() +
  xlim(c(0, 2)) +
  annotate(
    "text",
    x = 2,
    y = 3,
    label = paste0("P(R[Bayes] <= r[aa]) == ", round(postprob, 2)),
    hjust = 1,
    parse = TRUE,
    size = 2.5
  ) +
  ylab("Density") +
  xlab("")


# ---- Analysis for 1800 - 1925 ----

# Create testing and training explanatory variables
yrs2 <- 1800:1925
train_yrs2 <-1800:1900
test_yrs2 <- 1901:1925

X_train2 <- forecast::fourier(
  ts(
    tree_dat$mean_rwi[tree_dat$year %in% train_yrs2],
    frequency = length(train_yrs2)
  ),
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
    y = tree_dat$mean_rwi[tree_dat$year %in% train_yrs2],
    N_new = nrow(X_test2),
    X_new = cbind(1, X_test2),
    tau0_phi = tau0(
      y = tree_dat$mean_rwi[tree_dat$year %in% train_yrs2],
      m0 = 1,
      M = 100,
      N = nrow(X_train2),
      fam = "gaussian"
    ),
    tau0_beta = tau0(
      y = tree_dat$mean_rwi[tree_dat$year %in% train_yrs2],
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
  warmup = 2000,
  cores = 4,
  control = list(adapt_delta = 0.99, max_treedepth = 15)
)

# extract posterior predictive draws
y_hat2 <- rstan::extract(sparse_fit2, pars = "y_rep")$y_rep

# combine into plotting dataframe
plot_df2 <- data.frame(
  y = c(
    tree_dat$mean_rwi,
    colMeans(y_hat2),
    rep(NA, nrow(tree_dat) - length(yrs2))
  ),
  low = c(
    rep(NA, nrow(tree_dat)),
    apply(y_hat2, 2, quantile, probs = 0.025),
    rep(NA, nrow(tree_dat) - length(yrs2))
  ),
  high = c(
    rep(NA, nrow(tree_dat)),
    apply(y_hat2, 2, quantile, probs = 0.975),
    rep(NA, nrow(tree_dat) - length(yrs2))
  ),
  year = rep(tree_dat$year, 2),
  source = rep(c("Observed", "RHS"), each = nrow(tree_dat))
)

plot_df2$y[plot_df2$source == "RHS" & plot_df2$year %in% train_yrs2] <- NA
plot_df2$low[plot_df2$source == "RHS" & plot_df2$year %in% train_yrs2] <- NA
plot_df2$high[plot_df2$source == "RHS" & plot_df2$year %in% train_yrs2] <- NA
plot_df2$y[!(plot_df2$year %in% yrs2)] <- NA

# repeat the process with auto.arima

auto_arima_fit2 <- fit_seasonal_arima_model(
  model_pars = list(
    n = nrow(X_train2),
    p = 10,
    X = matrix(1, nrow = nrow(X_train2), ncol = 1),
    y = tree_dat$mean_rwi[tree_dat$year %in% train_yrs2],
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

plot_df2 <- data.frame(
  y = c(
    rep(NA, length(train_yrs2)),
    fc_aarima2$point_forecast,
    rep(NA, nrow(tree_dat) - length(yrs2))
  ),
  low = c(
    rep(NA, length(train_yrs2)),
    fc_aarima2$lo_95,
    rep(NA, nrow(tree_dat) - length(yrs2))
  ),
  high = c(
    rep(NA, length(train_yrs2)),
    fc_aarima2$hi_95,
    rep(NA, nrow(tree_dat) - length(yrs2))
  ),
  year = tree_dat$year,
  source = rep("Auto-ARIMA", nrow(tree_dat))
) %>% rbind(plot_df2, .)

# omit the later time points
plot_df2$y[!(plot_df2$year %in% yrs2)] <- NA


### ---- Second plot ----
dat2_plot <- forecast_plot(
  plot_df2,
  horizon = 1900,
  col = PNWColors::pnw_palette("Sunset", 7)[c(2,5)]
) +
  theme(
    legend.position = "none",
    axis.text.x = element_blank()
  ) +
  xlab("") +
  ylim(c(0.5, 3)) +
  ggtitle("b)")

## ---- Get RMSE of the forecast ----

rmse2_aa <- sqrt(
  mean(
    (tree_dat$mean_rwi[tree_dat$year %in% test_yrs2] -
       fc_aarima2$point_forecast)^2,
    na.rm = TRUE
  )
)

rmse2_RHS <- data.frame(
  rmse = RMSE_bayes(
    obs = tree_dat$mean_rwi[tree_dat$year %in% test_yrs2],
    ppreds = y_hat2[, -(1:length(train_yrs2))]
  )
)

postprob2 <- mean(rmse2_RHS$rmse < rmse2_aa)

### ---- Second RMSE plot ----
rmse2_plot <- ggplot(rmse2_RHS, aes(x = rmse)) +
  geom_density(fill = PNWColors::pnw_palette("Sunset", 7)[2]) +
  geom_vline(
    xintercept = rmse2_aa,
    color = PNWColors::pnw_palette("Sunset", 7)[5]) +
  theme_classic() +
  ylab("Density") +
  xlab("") +
  xlim(c(0, 2)) +
  annotate(
    "text",
    x = 2,
    y = 10,
    label = paste0("P(R[Bayes] <= r[aa]) == ", round(postprob2, 2)),
    hjust = 1,
    parse = TRUE,
    size = 2.5
  )


# ---- Analysis for 1926 - 2012 ----

# Create testing and training explanatory variables
yrs3 <- 1926:2012
train_yrs3 <- 1926:1995
test_yrs3 <- 1996:2012

X_train3 <- forecast::fourier(
  ts(
    tree_dat$mean_rwi[tree_dat$year %in% train_yrs3],
    frequency = length(train_yrs3)
  ),
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
    y = tree_dat$mean_rwi[tree_dat$year %in% train_yrs3],
    N_new = nrow(X_test3),
    X_new = cbind(1, X_test3),
    tau0_phi = tau0(
      y = tree_dat$mean_rwi[tree_dat$year %in% train_yrs3],
      m0 = 1,
      M = 100,
      N = nrow(X_train3),
      fam = "gaussian"
    ),
    tau0_beta = tau0(
      y = tree_dat$mean_rwi[tree_dat$year %in% train_yrs3],
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

plot_df3 <- data.frame(
  y = c(
    tree_dat$mean_rwi,
    rep(NA, nrow(tree_dat) - length(yrs3)),
    colMeans(y_hat3)
  ),
  low = c(
    rep(NA, nrow(tree_dat)),
    rep(NA, nrow(tree_dat) - length(yrs3)),
    apply(y_hat3, 2, quantile, probs = 0.025)
  ),
  high = c(
    rep(NA, nrow(tree_dat)),
    rep(NA, nrow(tree_dat) - length(yrs3)),
    apply(y_hat3, 2, quantile, probs = 0.975)
  ),
  year = rep(tree_dat$year, 2),
  source = rep(c("Observed", "RHS"), each = nrow(tree_dat))
)

# replace training years with NA
plot_df3$y[plot_df3$source == "RHS" & plot_df3$year %in% train_yrs3] <- NA
plot_df3$low[plot_df3$source == "RHS" & plot_df3$year %in% train_yrs3] <- NA
plot_df3$high[plot_df3$source == "RHS" & plot_df3$year %in% train_yrs3] <- NA


# repeat the process with auto.arima

auto_arima_fit3 <- fit_seasonal_arima_model(
  model_pars = list(
    n = nrow(X_train3),
    p = 10,
    X = matrix(1, nrow = nrow(X_train3), ncol = 1),
    y = tree_dat$mean_rwi[tree_dat$year %in% train_yrs3],
    holdout = length(test_yrs3)
  ),
  freq = length(train_yrs3)
)

fc_aarima3 <- forecast(
  auto_arima_fit3$fit_ar,
  h = length(test_yrs3),
  xreg = X_test3[, colnames(X_test3) %in% auto_arima_fit3$terms]
) |> as.data.frame() |> janitor::clean_names()


plot_df3 <- data.frame(
  y = c(
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
  year = tree_dat$year,
  source = rep("Auto-ARIMA", nrow(tree_dat))
) %>% rbind(plot_df3, .)

plot_df3$y[!(plot_df3$year %in% yrs3)] <- NA

### ---- Third plot ----

dat3_plot <- forecast_plot(
  plot_df3,
  horizon = 1996,
  col = PNWColors::pnw_palette("Sunset", 7)[c(2,5)]
) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  xlab("Year") +
  ylim(c(0.5, 3)) +
  ggtitle("c)")

## ---- Get RMSE of the forecast ----

rmse3_aa <- sqrt(
  mean(
    (tree_dat$mean_rwi[tree_dat$year %in% test_yrs3] -
       fc_aarima3$point_forecast)^2,
    na.rm = TRUE
  )
)

rmse3_RHS <- data.frame(
  rmse = RMSE_bayes(
    obs = tree_dat$mean_rwi[tree_dat$year %in% test_yrs3],
    ppreds = y_hat3[, -(1:length(train_yrs3))]
  )
)

postprob3 <- mean(rmse3_RHS$rmse < rmse3_aa)

### ---- Third RMSE plot ----
rmse3_plot <- ggplot(rmse3_RHS, aes(x = rmse)) +
  geom_density(fill = PNWColors::pnw_palette("Sunset", 7)[2]) +
  geom_vline(
    xintercept = rmse3_aa,
    color = PNWColors::pnw_palette("Sunset", 7)[5]) +
  theme_classic() +
  ylab("Density") +
  xlab("RMSE") +
  xlim(c(0, 2)) +
  annotate(
    "text",
    x = 2,
    y = 3,
    label = paste0("P(R[Bayes] <= r[aa]) == ", round(postprob3, 2)),
    hjust = 1,
    parse = TRUE,
    size = 2.5
  )


# ---- Combine all plots ----
lo <- "
  111144
  222255
  333366
"

all_dat_plot + dat2_plot + dat3_plot +
  rmse_all_plot + rmse2_plot + rmse3_plot +
  plot_layout(
    design = lo
  )

ggsave(
  here::here("Figures/tree_forecasts_fourier.png"),
  width = 5.75,
  height = 5
)
