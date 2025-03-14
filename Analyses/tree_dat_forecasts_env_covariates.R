# Author: Dusty Gannon
# Created: 2019-08-19
# Last edited:
# Description: This script documents the tree growth model fits and forecasts
#              using predictors constructed from PRISM data.

# ---- Setup ----

# Import libraries
library(tidyverse)
devtools::load_all()

# Load data
tree_dat <- readRDS(here::here("SparseTS_prismdata/ca719_BM_Seq_1800_2012.rds"))
prism_wy <- readRDS(here::here("SparseTS_prismdata/prism_wateryear.rds"))
prism_winter <- readRDS(here::here("SparseTS_prismdata/prism_winter.rds"))
prism_summer <- readRDS(here::here("SparseTS_prismdata/prism_summer.rds"))

# compile stan models
rhs_reg <- rstan::stan_model(here::here("Stan/sparse_reg_FHS.stan"))



# ---- Data wrangling ----

# join all the seasonal data
prism <- left_join(
  prism_wy, prism_winter, by = join_by("wateryear" == "year")
) %>% left_join(
  prism_summer, by = join_by("wateryear" == "year")
)

prism <- rename(prism, year = wateryear)

# first, we want to standardize the columns
prism_std <- prism %>%
  drop_na() %>%
  mutate(
    across(
     wy_ppt:summer_vpdmax,
      .fns = ~ scale(.x)[,1]
    )
  )

# now lag the covariates of interest
prism_lagged <- lag_covariates(
  prism_std,
  names = names(prism_std)[-1],
  lags = 5,
  time_col = "year"
)

# average all trees for a given year,
# then subset for years in the prism data
tree_dat <- tree_dat %>%
  group_by(year) %>%
  summarise(
    mean_rwi = mean(rwi, na.rm = TRUE)
  ) %>% filter(
    year %in% prism_lagged$year
  )



# ---- Model 1: all years ----

train_yrs <- 1900:1990
test_yrs <- 1991:2012

train_rows <- which(tree_dat$year %in% train_yrs)
test_rows <- which(tree_dat$year %in% test_yrs)

# compile data for stan
dat_stan <- list(
  N = length(train_yrs),
  P0 = 1,
  P = ncol(prism_lagged),
  y = tree_dat$mean_rwi[train_rows],
  X = cbind(
    1,
    as.matrix(
      prism_lagged[train_rows, -which(names(prism_lagged) == "year")]
    )
  ),
  tau0 = tau0(
    y = tree_dat$mean_rwi[train_rows],
    m0 = 5,
    M = ncol(prism_lagged) - 1,
    N = length(train_yrs),
    fam = "gaussian"
  ),
  slab_scl = 0.5,
  slab_df = 6,
  N_new = length(test_yrs),
  X_new = cbind(
    1,
    as.matrix(
      prism_lagged[test_rows, -which(names(prism_lagged) == "year")]
    )
  )
)

# fit the model
rhs_fit_dat_all <- rstan::sampling(
  rhs_reg,
  data = dat_stan,
  cores = 4,
  iter = 5000,
  warmup = 4000,
  control = list(adapt_delta = 0.99, max_treedepth = 12)
)

# now fit the same model using stepAIC
aicdat_all <- left_join(tree_dat, prism_lagged) %>%
  select(!year)
init_all <- lm(mean_rwi ~ 1, data = aicdat_all[train_rows, ])

# define formula for scope
form_full <- formula(
  paste0(
    "mean_rwi ~ ",
    paste(names(aicdat_all)[-1], collapse = " + ")
  )
)

# stepwise model selection
aic_fit_dat_all <- MASS::stepAIC(
  init_all,
  scope = form_full,
  direction = "forward"
)

preds_aic_dat_all <- predict(aic_fit_dat_all, newdata = aicdat_all, se = T)

## ---- Combine observed and predicted into one dataframe ----
# extract draws
y_pred <- rstan::extract(fit_dat_all, pars = "y_rep")$y_rep
beta_post <- rstan::extract(fit_dat_all, pars = "beta")$beta

# compute residual variance
aic_sigma2 <- summary(aic_fit_dat_all)$sigma^2

df_fcplot_all <- data.frame(
  year = rep(tree_dat$year, 3),
  y = c(tree_dat$mean_rwi, colMeans(y_pred), preds_aic_dat_all$fit),
  low = c(
    rep(NA, nrow(tree_dat)),
    apply(y_pred, 2, quantile, probs = 0.025),
    preds_aic_dat_all$fit - 2 * sqrt(preds_aic_dat_all$se.fit^2 + aic_sigma2)
  ),
  high = c(
    rep(NA, nrow(tree_dat)),
    apply(y_pred, 2, quantile, probs = 0.975),
    preds_aic_dat_all$fit + 2 * sqrt(preds_aic_dat_all$se.fit^2 + aic_sigma2)
  ),
  source = rep(c("Observed", "RHS", "Auto-ARIMA"), each = nrow(tree_dat))
)

# now set the training region to NA
df_fcplot_all$y[
  df_fcplot_all$source != "Observed" &
    df_fcplot_all$year %in% train_yrs
] <- NA

df_fcplot_all$low[
  df_fcplot_all$source != "Observed" &
    df_fcplot_all$year %in% train_yrs
] <- NA

df_fcplot_all$high[
  df_fcplot_all$source != "Observed" &
    df_fcplot_all$year %in% train_yrs
] <- NA

### ---- Forecast plot function ----
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

# set colors
my_colors <- PNWColors::pnw_palette("Sunset", 7)[c(2,6)]

### ---- Forecast plot, all data ----
fc_plot_all <- forecast_plot(df_fcplot_all, horizon = 1990, col = my_colors)


### ---- Coefficients plot, all data ----

# extract coefficient estimates from each method
beta_post_all <- rstan::extract(fit_dat_all, pars = "beta")$beta

# create dataframe from aic approach
beta_aic <- coef(aic_fit_dat_all)
ses_aic <- vcov(aic_fit_dat_all) |> diag() |> sqrt()

df_estims_aic_dat_all <- data.frame(
  var = names(beta_aic),
  estim = beta_aic,
  low = beta_aic - 2 * ses_aic,
  high = beta_aic + 2 * ses_aic,
  method = "AIC"
)
df_estims_aic_dat_all$var <- str_remove_all(
  df_estims_aic_dat_all$var,
  pattern = "[()]"
)

# now create plotting dataframe
df_estims_plot_all <- data.frame(
  var = c("Intercept", names(aicdat_all)[-1]),
  estim = colMeans(beta_post_all),
  low = apply(beta_post_all, 2, quantile, probs = 0.025),
  high = apply(beta_post_all, 2, quantile, probs = 0.975),
  method = "RHS"
)

# expand aic df
df_estims_aic_dat_all <- df_estims_plot_all %>%
  select(var) %>%
  left_join(., df_estims_aic_dat_all)

df_estims_aic_dat_all$method <- "AIC"
df_estims_aic_dat_all[is.na(df_estims_aic_dat_all)] <- 0

# now stack the two dataframes
df_estims_plot_all <- rbind(
  df_estims_plot_all,
  df_estims_aic_dat_all
)

coefs_plot_all <- df_estims_plot_all %>% filter(var != "Intercept") %>%
  ggplot(., aes(x = var, y = estim, color = method, fill = method)) +
    geom_col(position = "dodge") +
    coord_polar() +
    theme_void() +
    theme(
      axis.text.x = element_blank()
    ) +
  scale_color_manual(values = rev(my_colors)) +
  scale_fill_manual(values = rev(my_colors))



# ---- Analysis for 1926 - 2012 ----

yrs2 <- 1926:2012
train_yrs2 <- 1926:1995
test_yrs2 <- 1996:2012
rows2 <- which(tree_dat$year %in% yrs2)
train_rows2 <- which(tree_dat$year %in% train_yrs2)
test_rows2 <- which(tree_dat$year %in% test_yrs2)

# compile data for stan
dat_stan2 <- list(
  N = length(train_yrs2),
  P0 = 1,
  P = ncol(prism_lagged),
  y = tree_dat$mean_rwi[train_rows2],
  X = cbind(
    1,
    as.matrix(
      prism_lagged[train_rows2, -which(names(prism_lagged) == "year")]
    )
  ),
  tau0 = tau0(
    y = tree_dat$mean_rwi[train_rows2],
    m0 = 5,
    M = ncol(prism_lagged) - 1,
    N = length(train_yrs2),
    fam = "gaussian"
  ),
  slab_scl = 0.5,
  slab_df = 6,
  N_new = length(test_yrs2),
  X_new = cbind(
    1,
    as.matrix(
      prism_lagged[test_rows2, -which(names(prism_lagged) == "year")]
    )
  )
)

# fit the model
rhs_fit2 <- rstan::sampling(
  rhs_reg,
  data = dat_stan2,
  cores = 4,
  iter = 5000,
  warmup = 4000,
  control = list(adapt_delta = 0.99, max_treedepth = 12)
)

# now fit the same model using stepAIC

# fit the initial model
init2 <- lm(mean_rwi ~ 1, data = aicdat_all[train_rows2, ])

# stepwise model selection
aic_fit_dat2 <- MASS::stepAIC(
  init2,
  scope = form_full,
  direction = "forward"
)

# # compute residual variance
# aic_sigma2 <- summary(aic_fit_dat_all)$sigma^2

y_pred2 <- rstan::extract(rhs_fit2, pars = "y_rep")$y_rep
preds_aic2 <- predict(aic_fit_dat2, newdat = aicdat_all[rows2, ], se = T)

df_fcplot2 <- data.frame(
  year = rep(tree_dat$year[rows2], 3),
  y = c(tree_dat$mean_rwi[rows2], colMeans(y_pred2), preds_aic2$fit),
  low = c(
    rep(NA, length(yrs2)),
    apply(y_pred2, 2, quantile, probs = 0.025),
    rep(NA, length(yrs2))
  ),
  high = c(
    rep(NA, length(yrs2)),
    apply(y_pred2, 2, quantile, probs = 0.975),
    rep(NA, length(yrs2))
  ),
  source = rep(c("Observed", "RHS", "Auto-ARIMA"), each = length(yrs2))
)

## ---- Forecast plot for set 2 ----

fc_plot2 <- forecast_plot(df_fcplot2, horizon = 1995, col = my_colors)


