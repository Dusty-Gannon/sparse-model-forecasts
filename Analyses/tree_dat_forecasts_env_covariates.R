# Author: Dusty Gannon
# Created: 2 March 2025
# Last edited: 10 April 2025
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
      values = c("Observed" = "black", "RHS" = col[1], "stepAIC + auto.arima" = col[2])
    ) +
    scale_fill_manual(
      values = c("Observed" = "black", "RHS" = col[1], "stepAIC + auto.arima" = col[2])
    ) +
    scale_alpha_manual(
      values = c("Observed" = 0, "RHS" = 0.3, "stepAIC + auto.arima" = 0.4),
    ) +
    scale_linewidth_manual(
      values = c("Observed" = 1, "RHS" = 0.5, "stepAIC + auto.arima" = 0.5)
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


# ---- Model 1: all years ----

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
  source = rep(c("Observed", "RHS", "stepAIC + auto.arima"), each = nrow(tree_dat))
)

# # now set the training region to NA
# df_fcplot_all$y[
#   df_fcplot_all$source != "Observed" &
#     df_fcplot_all$year %in% train_yrs
# ] <- NA
#
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

# # don't need to rescale because AIC mod fit on scaled columns, too
# beta_post_all <- cbind(
#   beta_post_all[,1],
#   sweep(beta_post_all[,-1], 2, tr_sds, FUN = "/")
# )

# create dataframe from aic approach
beta_aic <- coef(aic_aarima_fit_all)
ses_aic <- vcov(aic_aarima_fit_all) |> diag() |> sqrt()
nu_aic <- length(train_rows) - length(beta_aic) -
  length(aic_aarima_fit_all$model$phi) -
  length(aic_aarima_fit_all$model$theta)
t_star_all <- qt(0.975, nu_aic)

df_estims_aic_dat_all <- data.frame(
  var = names(beta_aic),
  estim = beta_aic,
  low = beta_aic - t_star_all * ses_aic,
  high = beta_aic + t_star_all * ses_aic,
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


# ---- Model 2: Analysis for 1926-2012 ----

train_yrs2 <- 1926:1995
test_yrs2 <- 1996:2012

train_rows2 <- which(prism$year %in% train_yrs2)
test_rows2 <- which(prism$year %in% test_yrs2)

# split first, then use the *training-period-specific* means/sds to
# standardize training and testing data (same approach as Model 1 -- this
# model's training window (1926-1995) differs from Model 1's (1895-1990),
# so it needs its own standardization constants rather than reusing
# tr_means/tr_sds/prism_lagged from Model 1)
prism_train2 <- prism[train_rows2, ]
prism_test2 <- prism[test_rows2, ]

tr_means2 <- select(prism_train2, !year) %>%
  colMeans()
tr_sds2 <- select(prism_train2, !year) %>%
  apply(., 2, sd)

prism_train2_std <- prism_train2
prism_test2_std <- prism_test2
for(var in names(tr_means2)){
  prism_train2_std[var] = (prism_train2_std[var] - tr_means2[var]) / tr_sds2[var]
  prism_test2_std[var] = (prism_test2_std[var] - tr_means2[var]) / tr_sds2[var]
}

# now recombine
prism_std2 <- rbind(prism_train2_std, prism_test2_std)

# now lag the covariates of interest
prism_lagged2 <- lag_covariates(
  prism_std2,
  names = names(prism_std2)[-1],
  lags = 5,
  time_col = "year"
) %>% drop_na()

# reset training and testing years based on dropping some rows from the
# lagged variables (the first few years of the 1926-1995 training window
# lose their lags, since this split -- unlike Model 1's -- doesn't include
# the pre-1926 rows needed to compute them)
train_yrs2 <- max(train_yrs2[1], prism_lagged2$year[1]):1995
test_yrs2 <- 1996:min(test_yrs2[length(test_yrs2)], prism_lagged2$year[nrow(prism_lagged2)])

# prism_lagged2 and tree_dat do NOT share the same row-to-year
# correspondence: prism_lagged2 starts wherever the lag trimming above
# left it (e.g. 1931), while tree_dat/aicdat_all still start from Model
# 1's much earlier range. So each object needs its own row indices,
# matched by *year value*, not by row position.
train_rows2_X <- which(prism_lagged2$year %in% train_yrs2)
test_rows2_X <- which(prism_lagged2$year %in% test_yrs2)

train_rows2 <- which(tree_dat$year %in% train_yrs2)
test_rows2 <- which(tree_dat$year %in% test_yrs2)
rows2 <- which(tree_dat$year %in% c(train_yrs2, test_yrs2))
yrs2 <- tree_dat$year[rows2]

X_train2 <- cbind(
  1,
  as.matrix(
    prism_lagged2[train_rows2_X, -which(names(prism_lagged2) == "year")]
  )
)
X_test2 <- cbind(
  1,
  as.matrix(
    prism_lagged2[test_rows2_X, -which(names(prism_lagged2) == "year")]
  )
)

# compile data for stan
dat_stan2 <- list(
  N = length(train_rows2),
  P_0 = 1,
  P = ncol(X_train2),
  p = 10,
  y = tree_dat$mean_rwi[train_rows2],
  X = X_train2,
  tau0_phi = 1 / (ncol(X_train2) - 1 - 1) * length(train_rows2)^(-0.5),
  slab_scl_phi = 0.4,
  slab_df_phi = 10,
  tau0_beta = 5 / (ncol(X_train2) - 5) * length(train_rows2)^(-0.5),
  slab_scl_beta = 0.5,
  slab_df_beta = 6,
  N_new = length(test_rows2),
  X_new = X_test2
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

## ---- Combine all results into df for plotting ----
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
  source = rep(c("Observed", "RHS", "stepAIC + auto.arima"), each = length(yrs2))
)

# replace the training periods with NAs
# df_fcplot2$y[
#   df_fcplot2$source != "Observed" &
#     df_fcplot2$year %in% train_yrs2
# ] <- NA

df_fcplot2$low[df_fcplot2$year %in% train_yrs2] <- NA
df_fcplot2$high[df_fcplot2$year %in% train_yrs2] <- NA



## ---- Coefficients plot dataframe for set 2----

# extract coefficient estimates from each method
beta_post2 <- rstan::extract(rhs_fit2, pars = "beta")$beta

# create dataframe from aic approach
beta_aic2 <- coef(aic_fit_dat2)
ses_aic2 <- rep(NA, length(beta_aic2))

df_estims_aic2 <- data.frame(
  var = names(beta_aic2),
  estim = beta_aic2,
  low = beta_aic2 - 2 * ses_aic2,
  high = beta_aic2 + 2 * ses_aic2,
  method = "AIC"
)
df_estims_aic2$var <- str_remove_all(
  df_estims_aic2$var,
  pattern = "[()]"
)

# now create plotting dataframe
df_estims_rhs2 <- data.frame(
  var = c("Intercept", names(aicdat_all)[-1]),
  estim = colMeans(beta_post2),
  low = apply(beta_post2, 2, quantile, probs = 0.025),
  high = apply(beta_post2, 2, quantile, probs = 0.975),
  method = "RHS"
)

# expand aic df
df_estims_aic2 <- df_estims_rhs2 %>%
  select(var) %>%
  left_join(., df_estims_aic2)

df_estims_aic2$method <- "AIC"
df_estims_aic2[is.na(df_estims_aic2)] <- 0

# now stack the two dataframes
df_estims_plot2 <- rbind(
  df_estims_rhs2,
  df_estims_aic2
)

# now create some labels
df_estims_plot2 <- df_estims_plot2 %>%
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
  )

# remove the intercept rows
df_estims_plot2 <- df_estims_plot2 %>%
  filter(var != "Intercept")

# count number of unique labels
n_labs <- length(unique(df_estims_plot2$var_lab)) *
  length(unique(df_estims_plot2$agg_period))

lab_df <- tibble(
  var_lab = unique(df_estims_plot2$var_lab) %>%
    rep(length(unique(df_estims_plot2$agg_period))),
  x = seq(0, n_labs - 1) * 6 + 1
)

lab2_df <- tibble(
  var_lab = unique(df_estims_plot2$agg_period),
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

df_estims_RHS2 <- df_estims_plot2 %>%
  filter(var != "Intercept") %>%
  filter(method == "RHS")

# ---- Create figures ----

fc_plot_all <- forecast_plot(df_fcplot_all, horizon = 1990, col = my_colors)

fc_plot2 <- forecast_plot(df_fcplot2, horizon = 1995, col = my_colors)

library(patchwork)
# add a and b panel labels


## ---- Final forecast figures ----

fc_plot_all_final <- fc_plot_all +
  theme(
    axis.text.x = element_blank(),
    legend.position = "top"
  ) +
  xlab("") +
  ylim(c(0, 4)) +
  ggtitle("a)")

fc_plot2_final <- fc_plot2 +
  xlim(c(min(tree_dat$year), max(tree_dat$year))) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  ylim(c(0, 4)) +
  xlab("") +
  ggtitle("b)")

fc_plot_all_final / fc_plot2_final

ggsave(
  filename = here::here("Figures/tree_growth_forecasts_env_covs.pdf"),
  width = 5,
  height = 5,
  device = "pdf",
  units = "in",
  dpi = 300
)

## ---- Final coefficients plot ----

sigfigs <- function(breaks){
  sprintf("%.2f", breaks)
}

### ---- Coefficients plot for 1900-2012 ----

coef_plot_all_rhs <- df_estims_plot_all %>%
  filter(var != "Intercept" & method == "RHS") %>%
  coef_plot(xlabs = F) +
  ggtitle("1900 - 1990", subtitle = "a) RHS sparse model") +
  ylab(expression(hat(beta))) +
  scale_y_continuous(labels = sigfigs)

coef_plot_all_aic <- df_estims_plot_all %>%
  filter(var != "Intercept" & method == "AIC") %>%
  coef_plot(., xlabs = F, ylims = c(-1, 1)) +
  ggtitle("b) Stepwise AIC") +
  theme(plot.title = element_text(size = 11)) +
  ylab(expression(hat(beta))) +
  scale_y_continuous(labels = sigfigs)

### ---- Coefficients plot for 1926-2012 ----

coef_plot_rhs2 <- coef_plot(
  df_estims_RHS2,
  xlabs = F
) +
  ggtitle("1931 - 1995", subtitle = "c) RHS sparse model") +
  ylab(expression(hat(beta)))

coef_plot_aic2 <- coef_plot(
  df_estims_aic2 %>%
    filter(var != "Intercept"),
  ypos_labs = c(-0.05, -0.68),
  ylims = c(
    min(df_estims_aic2$estim),
    max(df_estims_aic2$estim)
  )
) +
  ggtitle("d) Stepwise AIC") +
  theme(plot.title = element_text(size = 11)) +
  ylab(expression(hat(beta)))


# combine all top-to-bottom

coef_plot_all_rhs / coef_plot_all_aic /
  coef_plot_rhs2 / coef_plot_aic2


ggsave(
  filename = here::here("Figures/tree_growth_coefs.pdf"),
  width = 6.5,
  height = 9,
  device = "pdf",
  units = "in",
  dpi = 300
)
