
# Description: Tree growth model fits and forecasts using a trend + Fourier
#              basis as predictors. A single train/test window (1895-1990 /
#              1991-2012) is used to compare a regularized horseshoe (RHS)
#              Bayesian AR(p)-errors model against stepAIC term selection
#              with auto-ARIMA errors.

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
train_yrs <- 1900:1990
test_yrs <- 1991:tree_dat$year[nrow(tree_dat)]
yrs <- c(train_yrs, test_yrs)

# now create fourier terms
# intercept column added later in data construction
X_train <- forecast::fourier(
  ts(
    tree_dat$mean_rwi[tree_dat$year %in% train_yrs],
    frequency = length(train_yrs)
  ),
  K = 30
)

X_test <- continue_fourier(
  X_train,
  h = length(test_yrs),
  n = length(train_yrs)
)

# now add a trend and trend^2 component
X_train <- cbind(
  trend = (train_yrs - mean(train_yrs)) / (sd(train_yrs) / sd(X_train[,ncol(X_train)])),
  X_train
)
X_train <- cbind(
  X_train[,1],
  trend2 = X_train[,1]^2,
  X_train[, 2:ncol(X_train)]
)
# add to test data as well
X_test <- cbind(
  trend = (test_yrs - mean(train_yrs)) / (sd(train_yrs) / sd(X_train[,ncol(X_train)])),
  X_test
)
X_test <- cbind(
  X_test[,1],
  trend2 = X_test[,1]^2,
  X_test[,2:ncol(X_test)]
)

# ---- Sparse model fitting ----

# compile stan model
sparse_mod <- rstan::stan_model(
  here::here("Stan/AR-p_err.stan")
)

# create list of objects that are shared across all models
# get estimate of tau0
tau_0_phi <- tau0(
  y = tree_dat$mean_rwi[tree_dat$year %in% train_yrs],
  m0 = 2,
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
  p = 15,
  tau0_phi = tau_0_phi,
  slab_scl_phi = 0.3,
  slab_df_phi = 15,
  tau0_beta = tau_0_beta,
  slab_scl_beta = 0.1,
  slab_df_beta = 10,
  y = tree_dat$mean_rwi[tree_dat$year %in% train_yrs],
  N_new = nrow(X_test),
  X_new = cbind(1, X_test)
)

# fit the RHS model
sparse_fit <- rstan::sampling(
  sparse_mod,
  data = dat_stan,
  chains = 4,
  iter = 8000,
  warmup = 6000,
  cores = 4,
  control = list(adapt_delta = 0.99, max_treedepth = 15)
)

# extract posterior predictive draws
y_hat <- rstan::extract(sparse_fit, pars = "y_rep")$y_rep

plot_df_allyrs <- data.frame(
  y = c(tree_dat$mean_rwi[tree_dat$year %in% yrs], colMeans(y_hat)),
  low = c(
    rep(NA, length(yrs)),
    apply(y_hat, 2, quantile, probs = 0.025)
  ),
  high = c(
    rep(NA, length(yrs)),
    apply(y_hat, 2, quantile, probs = 0.975)
  ),
  year = rep(yrs, 2),
  source = rep(c("Observed", "RHS"), each = length(yrs))
)

plot_df_allyrs$low[plot_df_allyrs$source == "RHS" & plot_df_allyrs$year <= max(train_yrs)] <- NA
plot_df_allyrs$high[plot_df_allyrs$source == "RHS" & plot_df_allyrs$year <= max(train_yrs)] <- NA

# colors: RHS, stepAIC + auto-ARIMA
my_colors <- PNWColors::pnw_palette("Cascades", 5)[c(5, 3)]

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
    scale_color_manual(
      values = c("Observed" = "black", "RHS" = col[1], "stepAIC + auto-ARIMA" = col[2])
    ) +
    scale_fill_manual(
      values = c("Observed" = "black", "RHS" = col[1], "stepAIC + auto-ARIMA" = col[2])
    ) +
    scale_linewidth_manual(
      values = c("Observed" = 0.8, "RHS" = 0.5, "stepAIC + auto-ARIMA" = 0.5)
    ) +
    scale_alpha_manual(
      values = c("Observed" = 1, "RHS" = 0.2, "stepAIC + auto-ARIMA" = 0.4)
    ) +
    theme_classic() +
    theme(legend.title = element_blank()) +
    ylab("RWI")
}

# now fit auto.arima model with stepAIC-selected terms

auto_arima_fit <- fit_seasonal_arima_model(
  model_pars = list(
    p = 15,
    X = NULL,
    y = tree_dat$mean_rwi[tree_dat$year %in% yrs],
    holdout = length(test_yrs)
  ),
  K = 30,
  freq = length(train_yrs),
  xreg = rbind(X_train[,1:2], X_test[,1:2])
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
    auto_arima_fit$fit_ar$fitted,
    fc_aarima$point_forecast
  ),
  low = c(rep(NA, nrow(X_train)), fc_aarima$lo_95),
  high = c(rep(NA, nrow(X_train)), fc_aarima$hi_95),
  year = yrs,
  source = rep("stepAIC + auto-ARIMA", length(yrs))
) %>% rbind(plot_df_allyrs, .)

### ---- Forecast plot ----
all_dat_plot <- forecast_plot(
  plot_df_allyrs,
  horizon = max(train_yrs),
  col = my_colors
) +
  xlab("Year") +
  ylim(c(0.5, 2.7)) +
  theme(
    legend.position = "none"
  ) +
  ggtitle("b) Endogenous only")


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
    ppreds = y_hat[, (ncol(y_hat) - length(test_yrs) + 1):ncol(y_hat)]
  )
)
## ---- make figure for RMSE ----

postprob <- mean(rmse_all_RHS$rmse < rmse_all_aa)

rmse_all_plot <- ggplot(rmse_all_RHS, aes(x = rmse)) +
  geom_density(fill = my_colors[1]) +
  geom_vline(xintercept = rmse_all_aa, color = my_colors[2], linewidth = 1.2) +
  theme_classic() +
  xlim(c(0, 1.5)) +
  ylim(c(0, 5.5)) +
  annotate(
    "text",
    x = 1.5,
    y = 4.2,
    label = paste0("P(R[Bayes] <= r[aa]) == ", round(postprob, 2)),
    hjust = 1,
    parse = TRUE,
    size = 2.5
  ) +
  ylab("Density") +
  xlab("Prediction RMSE")

## ---- In-text figures ----

wide_preds_all <- plot_df_allyrs %>%
  filter(year %in% test_yrs) %>%
  pivot_wider(
    id_cols = year,
    names_from = source,
    values_from = c(y, low, high)
  )

writeLines(c(
  "Prediction interval coverage for RHS:",
  mean(
    wide_preds_all$low_RHS < wide_preds_all$y_Observed &
      wide_preds_all$y_Observed < wide_preds_all$high_RHS
  )
))

writeLines(c(
  "Prediction interval coverage for stepAIC + auto-ARIMA:",
  with(wide_preds_all, {
    mean(
      `low_stepAIC + auto-ARIMA` < y_Observed &
        y_Observed < `high_stepAIC + auto-ARIMA`
    )
  })
))

# ---- Combine forecast + RMSE panels ----

panel_fourier <- all_dat_plot + rmse_all_plot +
  plot_layout(widths = c(2, 1))

# save the combined panel object so it can be cheaply recombined with the
# env-covariates analysis's panel into one 4-panel figure without rerunning
# the Stan sampling in either script
saveRDS(
  panel_fourier,
  here::here("Figures/tree_forecasts_fourier_panel.rds")
)

panel_fourier

ggsave(
  here::here("Figures/tree_forecasts_fourier.pdf"),
  plot = panel_fourier,
  width = 7,
  height = 4,
  device = "pdf",
  dpi = 300
)
