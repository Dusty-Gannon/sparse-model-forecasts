# Author: Dusty Gannon
# Created: 2 March 2025
# Last edited: 04 April 2025
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
rhs_reg <- rstan::stan_model(here::here("Stan/sparse_reg_FHS.stan"))

# set colors
my_colors <- PNWColors::pnw_palette("Sunset", 7)[c(2,6)]

# ---- User-defined functions ----

## ---- Forecast plot function ----
forecast_plot <- function(df, horizon, col){
  ggplot(
    df,
    aes(x = year, y = y, colour = source)
  ) +
    geom_line(aes(linewidth = source)) +
    geom_ribbon(
      aes(ymin = low, ymax = high, fill = source, alpha = source),
      linetype = 0,
    ) +
    geom_vline(xintercept = horizon, linetype = "dashed", color = "brown") +
    geom_vline(xintercept = 1926, linetype = "dashed", color = "grey") +
    scale_color_manual(
      values = c("Observed" = "black", "RHS" = col[1], "stepAIC" = col[2])
    ) +
    scale_fill_manual(
      values = c("Observed" = "black", "RHS" = col[1], "stepAIC" = col[2])
    ) +
    scale_alpha_manual(
      values = c("Observed" = 0, "RHS" = 0.3, "stepAIC" = 0.4),
    ) +
    scale_linewidth_manual(
      values = c("Observed" = 1, "RHS" = 0.5, "stepAIC" = 0.5)
    ) +
    theme_classic() +
    theme(legend.title = element_blank()) +
    ylab("RWI")
}


# Coefficients plot function

coef_plot <- function(
    dat,
    ylabs = TRUE,
    xlims = c(-0.5, 0.5),
    xpos_labs = c(-0.05, -0.3)
  ) {

  strips <- strip_df
  strips$xmin <- xlims[1]
  strips$xmax <- xlims[2]

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
    geom_point(aes(x = estim, y = var), size = 1) +
    geom_errorbarh(aes(
      y = var,
      xmin = low,
      xmax = high
    ), height = 0, color = "gray2", alpha = 0.8) +
    theme_classic() +
    theme(
      axis.text.y = element_blank()
    ) +
    ylab("") +
    xlab(expression(hat(beta)))

  if(ylabs){
    p <- p + theme(
      plot.margin = unit(c(1, 1, 1, 6), "lines")
    ) +
      coord_cartesian(xlim = xlims, clip = "off") +
      geom_text(
        data = lab_df,
        aes(
          x = I(xpos_labs[1]),
          y = y,
          label = var_lab
        ),
        hjust = 1,
        size = 2.5) +
      geom_text(
        data = lab2_df,
        aes(
          x = I(xpos_labs[2]),
          y = y,
          label = var_lab
        ),
        hjust = 0.5,
        size = 4,
        angle = 90
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
y_pred <- rstan::extract(rhs_fit_dat_all, pars = "y_rep")$y_rep
beta_post <- rstan::extract(rhs_fit_dat_all, pars = "beta")$beta

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
  source = rep(c("Observed", "RHS", "stepAIC"), each = nrow(tree_dat))
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


# ---- Model 2: Analysis for 1926-2012 ----

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
  source = rep(c("Observed", "RHS", "stepAIC"), each = length(yrs2))
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
ses_aic2 <- vcov(aic_fit_dat2) |> diag() |> sqrt()

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
  y = seq(0, n_labs - 1) * 6 + 1
)

lab2_df <- tibble(
  var_lab = unique(df_estims_plot2$agg_period),
  y = seq(1, 3) * 42 - 21
)

fill_colors <- PNWColors::pnw_palette("Shuksan2", 5)[c(1:2,5)]

strip_df <- tibble(
  xmin = -0.5,
  xmax = 0.5,
  ymin = seq(
    from = 0,
    to = (n_labs - 1) * 6,
    by = 6
  ),
  ymax = seq(
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
  filename = here::here("Figures/tree_growth_forecasts_env_covs.png"),
  width = 6,
  height = 5,
  device = "png",
  units = "in",
  dpi = 300
)

## ---- Final coefficients plot ----

### ---- Coefficients plot for 1926-2012 ----
coef_plot_rhs2 <- coef_plot(
  df_estims_RHS2,
  xpos_labs = c(-0.05, -0.5)
) +
  ggtitle("a) RHS sparse model")

coef_plot_aic2 <- coef_plot(
  df_estims_aic2 %>%
    filter(var != "Intercept"),
  ylabs = F,
  xlims = c(
    min(df_estims_aic2$estim),
    max(df_estims_aic2$estim)
  )
) +
  ggtitle("b) Stepwise AIC")

# combine the two side-by-side
coef_plot_rhs2 + coef_plot_aic2

ggsave(
  filename = here::here("Figures/tree_growth_coefs_1926-2012.png"),
  width = 6.5,
  height = 6,
  device = "png",
  units = "in",
  dpi = 300
)

### ---- Coefficients plot for 1900-2012 ----
coef_plot_all_rhs <- df_estims_plot_all %>%
  filter(var != "Intercept" & method == "RHS") %>%
  coef_plot(xpos_labs = c(-0.05, -0.5)) +
  ggtitle("a) RHS sparse model")

coef_plot_all_aic <- df_estims_plot_all %>%
  filter(var != "Intercept" & method == "AIC") %>%
  coef_plot(., ylabs = F, xlims = c(-1, 1)) +
  ggtitle("b) Stepwise AIC")

# combine the two side-by-side
coef_plot_all_rhs + coef_plot_all_aic

ggsave(
  filename = here::here("Figures/tree_growth_coefs_1900-2012.png"),
  width = 6.5,
  height = 6,
  device = "png",
  units = "in",
  dpi = 300
)
