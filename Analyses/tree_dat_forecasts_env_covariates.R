# Author: Dusty Gannon
# Created: 2 March 2025
# Last edited: 16 November 2025
# Description: This script documents the tree growth model fits and forecasts
#              using predictors constructed from PRISM data. A single
#              train/test window (1895-1990 / 1991-2013) is used to compare
#              a regularized horseshoe (RHS) Bayesian model against stepAIC
#              variable selection with automatic autoregressive error model
#              selection.

# ---- Setup ----

# Import libraries
library(tidyverse)
library(patchwork)
devtools::load_all()

# Load data
tree_dat <- readRDS(here::here("SparseTS_prismdata/ca719_BM_Seq_1800_2012.rds"))
prism_wy <- readRDS(here::here("SparseTS_prismdata/prism_wateryear.rds"))
prism_winter <- readRDS(here::here("SparseTS_prismdata/prism_winter.rds"))
prism_summer <- readRDS(here::here("SparseTS_prismdata/prism_summer.rds"))

## ---- Global variables ----

# compile stan models
rhs_reg <- rstan::stan_model(here::here("Stan/AR-p_err.stan"))

# set colors
my_colors <- PNWColors::pnw_palette("Sunset", 7)[c(2,5)]

# ---- User-defined functions ----

## ---- Forecast plot function ----
forecast_plot <- function(df, horizon, col){
  ggplot(
    df,
    aes(x = year, y = y, colour = source)
  ) +
    geom_ribbon(
      aes(ymin = low, ymax = high, fill = source, alpha = source),
      linetype = 0,
    ) +
    geom_line(aes(linewidth = source)) +
    geom_vline(xintercept = horizon, linetype = "dashed", color = "brown") +
    geom_vline(xintercept = 1926, linetype = "dashed", color = "grey") +
    scale_color_manual(
      values = c("Observed" = "black", "RHS" = col[1], "stepAIC + auto-ARIMA" = col[2])
    ) +
    scale_fill_manual(
      values = c("Observed" = "black", "RHS" = col[1], "stepAIC + auto-ARIMA" = col[2])
    ) +
    scale_alpha_manual(
      values = c("Observed" = 0, "RHS" = 0.3, "stepAIC + auto-ARIMA" = 0.4),
    ) +
    scale_linewidth_manual(
      values = c("Observed" = 1, "RHS" = 0.5, "stepAIC + auto-ARIMA" = 0.5)
    ) +
    theme_classic() +
    theme(legend.title = element_blank()) +
    ylab("RWI")
}


# Coefficients plot function

coef_plot <- function(
    dat,
    xlabs = TRUE,
    ylims = c(-0.5, 0.5),
    ypos_labs = c(-0.05, -0.3)
  ) {

  strips <- strip_df
  strips$ymin <- ylims[1]
  strips$ymax <- ylims[2]

  p <- ggplot(dat) +
    geom_rect(
      data = strips,
      aes(
        xmin = xmin,
        xmax = xmax,
        ymin = ymin,
        ymax = ymax,
        alpha = I(alpha),
        fill = I(color)
      ),
      linetype = 0
    ) +
    geom_point(aes(x = var, y = estim), size = 1) +
    geom_errorbar(aes(
      x = var,
      ymin = low,
      ymax = high
    ), width = 0, color = "gray2", alpha = 0.8) +
    theme_classic() +
    theme(
      axis.text.x = element_blank(),
      axis.title.y = element_text(angle = 0, vjust = 0.5)
    ) +
    xlab("")

  if(xlabs){
    p <- p + theme(
      plot.margin = unit(c(1, 1, 6, 1), "lines")
    ) +
      coord_cartesian(ylim = ylims, clip = "off") +
      geom_text(
        data = lab_df,
        aes(
          x = x,
          y = I(ypos_labs[1]),
          label = var_lab
        ),
        hjust = 1,
        size = 2.5,
        angle = 45) +
      geom_text(
        data = lab2_df,
        aes(
          x = x,
          y = I(ypos_labs[2]),
          label = var_lab
        ),
        hjust = 0.5,
        size = 4
      )
  }

  return(p)

}

# ---- Data wrangling ----

# join all the seasonal data
prism <- left_join(
  prism_wy, prism_winter, by = join_by("wateryear" == "year")
) %>% left_join(
  prism_summer, by = join_by("wateryear" == "year")
)

prism <- rename(prism, year = wateryear)


# ---- Model fits: 1895-1990 training, 1991-2013 testing ----

train_yrs <- 1895:1990
test_yrs <- 1991:2013

train_rows <- which(prism$year %in% train_yrs)
test_rows <- which(prism$year %in% test_yrs)

# first split the data, then use same column
# standardization across training and testing split
prism_train <- prism[train_rows, ]
prism_test <- prism[test_rows, ]
# compute the sds
tr_means <- select(prism_train, !year) %>%
  colMeans()
tr_sds <- select(prism_train, !year) %>%
  apply(., 2, sd)

prism_train_std <- prism_train
prism_test_std <- prism_test
for(var in names(tr_means)){
  prism_train_std[var] = (prism_train_std[var] - tr_means[var]) / tr_sds[var]
  prism_test_std[var] = (prism_test_std[var] - tr_means[var]) / tr_sds[var]
}

# now recombine
prism_std <- rbind(prism_train_std, prism_test_std)

# now lag the covariates of interest
prism_lagged <- lag_covariates(
  prism_std,
  names = names(prism_std)[-1],
  lags = 5,
  time_col = "year"
) %>% drop_na()

# reset training and testing years based on dropping some from lagged variables
train_yrs <- prism_lagged$year[1]:1990
test_yrs <- 1991:prism_lagged$year[nrow(prism_lagged)]
train_rows <- which(prism_lagged$year %in% train_yrs)
test_rows <- which(prism_lagged$year %in% test_yrs)

# average all trees for a given year,
# then subset for years in the prism data
tree_dat <- tree_dat %>%
  group_by(year) %>%
  summarise(
    mean_rwi = mean(rwi, na.rm = TRUE)
  ) %>% filter(
    year %in% prism_lagged$year
  )

X_train <- cbind(
  1,
  as.matrix(
    prism_lagged[train_rows, -which(names(prism_lagged) == "year")]
  )
)
X_test <- cbind(
  1,
  as.matrix(
    prism_lagged[test_rows, -which(names(prism_lagged) == "year")]
  )
)

# compile data for stan
dat_stan <- list(
  N = length(train_yrs),
  P_0 = 1,
  P = ncol(X_train),
  p = 10,
  y = tree_dat$mean_rwi[train_rows],
  X = X_train,
  tau0_phi = 1 / (10 - 1) * length(train_rows)^(-0.5),
  slab_scl_phi = 0.4,
  slab_df_phi = 10,
  tau0_beta = (5 / ((ncol(X_train) - 1) - 5)) * length(train_rows)^(-0.5),
  slab_scl_beta = 0.5,
  slab_df_beta = 6,
  N_new = length(test_yrs),
  X_new = X_test
)

# fit the RHS model
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

selected_all <- colnames(model.matrix(aic_fit_dat_all))[-1]
X_test_aarima_all <- X_test[, selected_all]
# add intercept back
X_test_aarima_all <- cbind(
  "(Intercept)" = 1,
  X_test_aarima_all
)
X_train_aic <- X_train[, selected_all]

# add in the ARMA model using package forecast
aic_aarima_fit_all <- forecast::auto.arima(
  y = aicdat_all$mean_rwi[train_rows],
  max.order = 10,
  stationary = T,
  xreg = model.matrix(aic_fit_dat_all)
)

preds_aic_dat_all <- forecast::forecast(
  aic_aarima_fit_all,
  h = length(test_rows),
  xreg = X_test_aarima_all
)

## ---- Combine observed and predicted into one dataframe ----
# extract draws
y_pred <- rstan::extract(rhs_fit_dat_all, pars = "y_rep")$y_rep

df_fcplot_all <- data.frame(
  year = rep(tree_dat$year, 3),
  y = c(
    tree_dat$mean_rwi,
    colMeans(y_pred),
    c(as.vector(preds_aic_dat_all$fitted), as.vector(preds_aic_dat_all$mean))
  ),
  low = c(
    rep(NA, nrow(tree_dat)),
    c(
      rep(NA, length(train_rows)),
      apply(y_pred[, test_rows], 2, quantile, probs = 0.025)
    ),
    c(
      rep(NA, length(train_rows)),
      as.vector(preds_aic_dat_all$lower[, "95%"])
    )
  ),
  high = c(
    rep(NA, nrow(tree_dat)),
    c(
      rep(NA, length(train_rows)),
      apply(y_pred[, test_rows], 2, quantile, probs = 0.975)
    ),
    c(
      rep(NA, length(train_rows)),
      as.vector(preds_aic_dat_all$upper[, "95%"])
    )
  ),
  source = rep(c("Observed", "RHS", "stepAIC + auto-ARIMA"), each = nrow(tree_dat))
)

# now set the training region ribbons to NA (training fit is shown as a
# line only, uncertainty ribbons are reserved for the test/forecast period)
df_fcplot_all$low[
  df_fcplot_all$source != "Observed" &
    df_fcplot_all$year %in% train_yrs
] <- NA

df_fcplot_all$high[
  df_fcplot_all$source != "Observed" &
    df_fcplot_all$year %in% train_yrs
] <- NA


## ---- Coefficients plot dataframe ----

# extract coefficient estimates from each method
beta_post_all <- rstan::extract(rhs_fit_dat_all, pars = "beta")$beta

# RHS estimates
df_estims_rhs_all <- data.frame(
  var = c("Intercept", names(aicdat_all)[-1]),
  estim = colMeans(beta_post_all),
  low = apply(beta_post_all, 2, quantile, probs = 0.025),
  high = apply(beta_post_all, 2, quantile, probs = 0.975),
  method = "RHS"
)

# AIC estimates
beta_aic <- coef(aic_aarima_fit_all)
ses_aic <- vcov(aic_aarima_fit_all) |> diag() |> sqrt()
nu_aic <- length(train_rows) - length(beta_aic) -
  length(aic_aarima_fit_all$model$phi) -
  length(aic_aarima_fit_all$model$theta)
t_star_aic <- qt(0.975, nu_aic)

df_estims_aic_dat_all <- data.frame(
  var = names(beta_aic),
  estim = beta_aic,
  low = beta_aic - t_star_aic * ses_aic,
  high = beta_aic + t_star_aic * ses_aic,
  method = "AIC"
)
df_estims_aic_dat_all$var <- str_remove_all(
  df_estims_aic_dat_all$var,
  pattern = "[()]"
)

# expand aic df to the full variable list (variables stepAIC dropped get 0s)
df_estims_aic_dat_all <- df_estims_rhs_all %>%
  select(var) %>%
  left_join(., df_estims_aic_dat_all)
df_estims_aic_dat_all$method <- "AIC"
df_estims_aic_dat_all[is.na(df_estims_aic_dat_all)] <- 0

# now stack the two dataframes
df_estims_plot_all <- rbind(
  df_estims_rhs_all,
  df_estims_aic_dat_all
)

# now create labels for the coefficients plot and drop the intercept row
# (these are variable-name-based and don't depend on the analysis window)
df_estims_plot_all <- df_estims_plot_all %>%
  mutate(
    agg_period = case_when(
      str_detect(var, "wy_") ~ "Water-year",
      str_detect(var, "winter_") ~ "Winter",
      str_detect(var, "summer_") ~ "Summer"
    ),
    var_lab = case_when(
      str_detect(var, "ppt") ~ "Precip",
      str_detect(var, "tmin") ~ "min Temp",
      str_detect(var, "tmax") ~ "max Temp",
      str_detect(var, "tmean") ~ "mean Temp",
      str_detect(var, "td_mean") ~ "mean Dew Point",
      str_detect(var, "vpdmax") ~ "max VPD",
      str_detect(var, "vpdmin") ~ "min VPD",
    )
  ) %>%
  filter(var != "Intercept")

# count number of unique labels
n_labs <- length(unique(df_estims_plot_all$var_lab)) *
  length(unique(df_estims_plot_all$agg_period))

lab_df <- tibble(
  var_lab = unique(df_estims_plot_all$var_lab) %>%
    rep(length(unique(df_estims_plot_all$agg_period))),
  x = seq(0, n_labs - 1) * 6 + 1
)

lab2_df <- tibble(
  var_lab = unique(df_estims_plot_all$agg_period),
  x = seq(1, 3) * 42 - 21
)

fill_colors <- PNWColors::pnw_palette("Shuksan2", 5)[c(1:2,5)]

strip_df <- tibble(
  ymin = -0.5,
  ymax = 0.5,
  xmin = seq(
    from = 0,
    to = (n_labs - 1) * 6,
    by = 6
  ),
  xmax = seq(
    from = 6,
    to = n_labs * 6,
    by = 6
  ),
  alpha = rep(c(1, 0.7), length.out = n_labs),
  color = rep(fill_colors, each = n_labs / 3)
)


## ---- Get RMSE of the forecast ----

rmse_all_aa <- sqrt(
  mean(
    (tree_dat$mean_rwi[tree_dat$year %in% test_yrs] -
       preds_aic_dat_all$mean)^2,
    na.rm = TRUE
  )
)

rmse_all_RHS <- data.frame(
  rmse = RMSE_bayes(
    obs = tree_dat$mean_rwi[tree_dat$year %in% test_yrs],
    ppreds = y_pred[, which(tree_dat$year %in% test_yrs)]
  )
)

## ---- make figure for RMSE ----

postprob <- mean(rmse_all_RHS$rmse < rmse_all_aa)

rmse_all_plot <- ggplot(rmse_all_RHS, aes(x = rmse)) +
  geom_density(fill = PNWColors::pnw_palette("Sunset", 7)[2]) +
  geom_vline(
    xintercept = rmse_all_aa,
    color = PNWColors::pnw_palette("Sunset", 7)[5]) +
  theme_classic() +
  xlim(c(0, 1.5)) +
  annotate(
    "text",
    x = 1.5,
    y = 4,
    label = paste0("P(R[Bayes] <= r[aa]) == ", round(postprob, 2)),
    hjust = 1,
    parse = TRUE,
    size = 2.5
  ) +
  ylab("Density") +
  xlab("Prediction RMSE")



## ---- In-text figures ----

wide_preds_all <- df_fcplot_all %>%
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


# ---- Create figures ----

fc_plot_all <- forecast_plot(df_fcplot_all, horizon = 1990, col = my_colors)

fc_plot_all_final <- fc_plot_all +
  theme(legend.position = "top") +
  ylim(c(0, 3)) +
  xlab("Year")

## ---- Final forecast + RMSE figure ----

panel_env_covs <- fc_plot_all_final + rmse_all_plot +
  plot_layout(widths = c(2, 1))

# save the combined panel object so it can be cheaply recombined with the
# fourier-basis analysis's panel into one 4-panel figure without rerunning
# the Stan sampling in either script
saveRDS(
  panel_env_covs,
  here::here("Figures/tree_growth_forecasts_env_covs_panel.rds")
)

panel_env_covs

ggsave(
  filename = here::here("Figures/tree_growth_forecasts_env_covs.pdf"),
  plot = panel_env_covs,
  width = 7,
  height = 4,
  device = "pdf",
  units = "in",
  dpi = 300
)

## ---- Final coefficients plot ----

sigfigs <- function(breaks){
  sprintf("%.2f", breaks)
}

coef_plot_rhs <- df_estims_plot_all %>%
  filter(method == "RHS") %>%
  coef_plot(xlabs = F) +
  ggtitle("a) RHS sparse model") +
  ylab(expression(hat(beta))) +
  scale_y_continuous(labels = sigfigs)

coef_plot_aic <- df_estims_plot_all %>%
  filter(method == "AIC") %>%
  coef_plot(
    ypos_labs = c(-0.03, -0.25),
    ylims = c(
      min(df_estims_plot_all$low[df_estims_plot_all$method == "AIC"]),
      max(df_estims_plot_all$high[df_estims_plot_all$method == "AIC"])
  )) +
  ggtitle("b) Stepwise AIC") +
  ylab(expression(hat(beta))) +
  scale_y_continuous(labels = sigfigs)

# combine top-to-bottom
coef_plot_rhs / coef_plot_aic

ggsave(
  filename = here::here("Figures/tree_growth_coefs.pdf"),
  width = 6.5,
  height = 9,
  device = "pdf",
  units = "in",
  dpi = 300
)
