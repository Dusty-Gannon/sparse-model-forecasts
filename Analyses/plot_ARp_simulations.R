# Plot summary data from AR-p beta-p simulation runs
# These simulations iterate through different time series lengths
# and different standard deviations for the random innovations
# where there are sparse covariates, lags in one covariate, and sparse AR terms
# comparing the fits of regularized and not regularized models

library(tidyverse)
library(grDevices)


mod_cols <- c("#a52a2aff", "#33406fff")

dd <- read_csv('Data/aquatic_sim_data/ARp_err_sims_02_12_condensed.csv')
dd <- read_csv('Data/aquatic_sim_data/ARp_err_sims_01_03_condensed.csv') %>%
  bind_rows(dd)
dd <- read_csv('Data/aquatic_sim_data/ARp_err_sims_10_31_condensed.csv') %>%
  bind_rows(dd)
dd <- read_csv('Data/aquatic_sim_data/ARp_err_sims_04_22_condensed.csv') %>%
  bind_rows(dd)

dd <- dd %>%
  mutate(prior = case_when(model == 'not_reg' ~ 'Gaussian',
                           model == 'reg' ~ 'Horseshoe'))

divergent_trans_cap <- 20

# Plot divergent transitions:
total_runs <- dd %>%
  group_by(n, sigma, model, betas) %>%
  summarize(total_runs = n())

con_runs <- filter(dd, divergent_trans <= divergent_trans_cap) %>%
  group_by(n, sigma, model, betas) %>%
  summarize(runs = n())

png('Manuscript/Figures/ARp_err_model_convergence.png',
    width = 8, height = 5, units = 'in', res = 300)
  left_join(con_runs, total_runs) %>%
    mutate(converged_runs = runs/total_runs*100,
           N_betas = factor(betas)) %>%
    ggplot(aes(n, converged_runs, col = model, lty = N_betas)) +
    geom_line(size = 0.9) +
    scale_color_manual('Prior', values = mod_cols) +
    theme_classic() +
    # theme(legend.position = c(0.75, 0.25),
    #       legend.box = 'horizontal')+
    # xlim(60,150)+
    ylim(0,100)+
    ylab('Percent convergence') +
    xlab('Time series length')

dev.off()

# Subsample model fits so that the same number is in each category.
# min_size = max(min(c(con_runs$runs)), 50)
# min_size = min(c(con_runs$runs))
min_size_gauss = min(c(con_runs$runs[con_runs$model == 'gauss']))
good_fits_gauss <-  dd %>%
  group_by(mod_run) %>%
  summarize(divergent_trans = max(divergent_trans)) %>%
  ungroup() %>%
  filter(divergent_trans <= divergent_trans_cap)

min_size = min(c(con_runs$runs[con_runs$model == 'hs']))
good_fits <- dd %>%
  mutate(divergent_trans = case_when(model == 'gauss' ~ 0,
                                     TRUE ~ divergent_trans)) %>%
  group_by(mod_run) %>%
  summarize(divergent_trans = max(divergent_trans)) %>%
  ungroup() %>%
  filter(divergent_trans <= divergent_trans_cap)


dd <- tidyr::fill(dd, mod_run )

dat <- data.frame()

for(i in 1:nrow(total_runs)){
  tmp <- dd %>%
    filter(mod_run %in% good_fits$mod_run,
           n == total_runs$n[i],
           betas == total_runs$betas[i],
           model == total_runs$model[i])
  rows <- sample(1:nrow(tmp), min_size, replace = FALSE)

  dat <- bind_rows(dat, tmp[rows,])
}

dat <- filter(dat, model != 'gauss') %>%
  mutate(model = case_when(model == 'arima' ~ 'Seasonal \nAuto Arima',
                           model == 'hs' ~ 'Horseshoe \nPrior'))

# subsample only runs where the gaussian prior model converged.
gauss_dat <- data.frame()

for(i in 1:nrow(total_runs)){
  tmp <- dd %>%
    filter(mod_run %in% good_fits_gauss$mod_run,
           n == total_runs$n[i],
           betas == total_runs$betas[i],
           model == total_runs$model[i])
  rows <- sample(1:nrow(tmp), 43, replace = FALSE)

  gauss_dat <- bind_rows(gauss_dat, tmp[rows,])
}


# Plot model results:


# dat %>%
#   pivot_longer(cols = c('fr_true_pos', 'fr_false_pos', 'beta_true_pos', 'beta_false_pos', 'rmse_forecast'),
#                values_to = 'value', names_to = 'metric') %>%
#   # group_by(betas, n, model) %>%
#   # summarize(med_rmse = median(rmse_forecast, na.rm = TRUE),
#   #           lower = quantile(rmse_forecast, 0.025, na.rm = TRUE),
#   #           upper = quantile(rmse_forecast, 0.975, na.rm = TRUE))
# ggplot(aes( y = value, fill = model)) +
#   geom_violin(aes(x = model), alpha = 0.5, width = 1.4) +
#   geom_boxplot(aes(x = model), width = 0.1, position = position_nudge(0.2)) +
#   facet_wrap(metric ~ ., ncol = 1, scales = 'free') +
#   labs(x = "", y = "Prediction RMSE") +
#   coord_flip() +
#   theme_minimal() +
#   theme(legend.position = "bottom")
#   # scale_fill_manual(values = c("Stepwise AIC" = "blue", "Horseshoe" = "red"))


# Generate a grayscale palette with distinct shades

modify_alpha <- function(hex, alpha) {
  alpha_hex <- sprintf("%02X", round(alpha * 255))
  paste0(substr(hex, 1, 7), alpha_hex)
}

# Convert hexadecimal color to RGB
rgb_values <- col2rgb(mod_cols)


dark_factor <- 0.75
dark_cols <- rgb(rgb_values[1,1]*dark_factor,
                 rgb_values[2,1]*dark_factor,
                 rgb_values[3,1]*dark_factor,
                 maxColorValue = 255)
dark_cols <- c(dark_cols,
               rgb(rgb_values[1,2]*dark_factor,
                   rgb_values[2,2]*dark_factor,
                   rgb_values[3,2]*dark_factor,
                   maxColorValue = 255))

legend_data <- data.frame(
  Model = rep("Horseshoe \nPrior", 3),

  rmse_forecast = rep(0, 3),
  n = c(365, 730, 1095),
  betas = c(0,2,4),
  color = barcols#rep('white', 3)
)



png(file = "Manuscript/Figures/Fourier_seasonality_Arima_comparison.png",
    width = 5, height = 5, units = 'in', res = 300)

  par(mar=c(4,5,1,2))
  par(mfrow=c(2,1))
  # Define the data for each group
  arima_data <- list(dat$rmse_forecast[which(dat$n == 365 & dat$model == 'Seasonal \nAuto Arima')],
                     dat$rmse_forecast[which(dat$n == 730 & dat$model == 'Seasonal \nAuto Arima')],
                     dat$rmse_forecast[which(dat$n == 1095 & dat$model == 'Seasonal \nAuto Arima')])

  horseshoe_data <- list(dat$rmse_forecast[which(dat$n == 365 & dat$model == 'Horseshoe \nPrior')],
                         dat$rmse_forecast[which(dat$n == 730 & dat$model == 'Horseshoe \nPrior')],
                         dat$rmse_forecast[which(dat$n == 1095 & dat$model == 'Horseshoe \nPrior')])

  # Manually specify the positions for the violins
  x_positions <- c(1, 2, 3, 5, 6, 7)

  # Combine the data for vioplot
  vioplot::vioplot(arima_data[[1]], arima_data[[2]], arima_data[[3]],
                   horseshoe_data[[1]], horseshoe_data[[2]], horseshoe_data[[3]],
                   at = x_positions,
                   xlab = "", ylab = "", horizontal = TRUE, las = 1,
                   names = rep(c('1-year', '2-year', '3-year'), 2),
                   col = c(modify_alpha(mod_cols, 0.4)[1],
                           modify_alpha(mod_cols, 0.8)[1],
                           mod_cols[1],
                           modify_alpha(mod_cols, 0.4)[2],
                           modify_alpha(mod_cols, 0.8)[2],
                           mod_cols[2]),
                   pchMed = 20,
                   border = rep(dark_cols, each = 3),
                   rectCol = rep(dark_cols, each = 3),
                   lineCol = rep(dark_cols, each = 3),
                   colMed = rep(dark_cols, each = 3))

  axis(2, at = c(2, 6), labels = c("Arima", "Horseshoe"),
       las = 3, line = 3, tick = FALSE)  # Cluster labels

  # title(ylab = "Density of \nPrediction RMSE",line=6,cex.lab=1)
  title(xlab="Prediction RMSE",line=2.5,cex.lab=1)

  arima_data <- list(dat$rmse_beta[which(dat$n == 365 & dat$model == 'Seasonal \nAuto Arima')],
                     dat$rmse_beta[which(dat$n == 730 & dat$model == 'Seasonal \nAuto Arima')],
                     dat$rmse_beta[which(dat$n == 1095 & dat$model == 'Seasonal \nAuto Arima')])

  horseshoe_data <- list(dat$rmse_beta[which(dat$n == 365 & dat$model == 'Horseshoe \nPrior')],
                         dat$rmse_beta[which(dat$n == 730 & dat$model == 'Horseshoe \nPrior')],
                         dat$rmse_beta[which(dat$n == 1095 & dat$model == 'Horseshoe \nPrior')])

  # Combine the data for vioplot
  vioplot::vioplot(arima_data[[1]], arima_data[[2]], arima_data[[3]],
                   horseshoe_data[[1]], horseshoe_data[[2]], horseshoe_data[[3]],
                   at = x_positions,
                   xlab = "", ylab = "", horizontal = TRUE, las = 1,
                   names = rep(c('1-year', '2-year', '3-year'), 2),
                   col = c(modify_alpha(mod_cols, 0.4)[1],
                           modify_alpha(mod_cols, 0.8)[1],
                           mod_cols[1],
                           modify_alpha(mod_cols, 0.4)[2],
                           modify_alpha(mod_cols, 0.8)[2],
                           mod_cols[2]),
                   pchMed = 20,
                   border = rep(dark_cols, each = 3),
                   rectCol = rep(dark_cols, each = 3),
                   lineCol = rep(dark_cols, each = 3),
                   colMed = rep(dark_cols, each = 3))


  axis(2, at = c(2, 6), labels = c("Arima", "Horseshoe"),
       las = 3, line = 3, tick = FALSE)  # Cluster labels

  # title(ylab="Density of \nCovariate Effect RMSE",line=6,cex.lab=1)
  title(xlab="Covariate Effect RMSE",line= 2.5, cex.lab=1)

dev.off()



png(file = "Manuscript/Figures/Fourier_seasonality_Arima_comparison.png",
    width = 5, height = 5, units = 'in', res = 300)

  par(mar=c(4,5,1,2))
  par(mfrow=c(2,1))
  # Define the data for each group
  arima_data <- list(dat$rmse_forecast[which(dat$n == 365 & dat$model == 'Seasonal \nAuto Arima')],
                     dat$rmse_forecast[which(dat$n == 730 & dat$model == 'Seasonal \nAuto Arima')],
                     dat$rmse_forecast[which(dat$n == 1095 & dat$model == 'Seasonal \nAuto Arima')])

  horseshoe_data <- list(dat$rmse_forecast[which(dat$n == 365 & dat$model == 'Horseshoe \nPrior')],
                         dat$rmse_forecast[which(dat$n == 730 & dat$model == 'Horseshoe \nPrior')],
                         dat$rmse_forecast[which(dat$n == 1095 & dat$model == 'Horseshoe \nPrior')])

  # Manually specify the positions for the violins
  x_positions <- c(1, 2, 3, 5, 6, 7)

  # Combine the data for vioplot
  vioplot::vioplot(arima_data[[1]], arima_data[[2]], arima_data[[3]],
                   horseshoe_data[[1]], horseshoe_data[[2]], horseshoe_data[[3]],
                   at = x_positions,
                   xlab = "", ylab = "", horizontal = TRUE, las = 1,
                   names = rep(c('1-year', '2-year', '3-year'), 2),
                   col = c(modify_alpha(mod_cols, 0.4)[1],
                           modify_alpha(mod_cols, 0.8)[1],
                           mod_cols[1],
                           modify_alpha(mod_cols, 0.4)[2],
                           modify_alpha(mod_cols, 0.8)[2],
                           mod_cols[2]),
                   pchMed = 20,
                   border = rep(dark_cols, each = 3),
                   rectCol = rep(dark_cols, each = 3),
                   lineCol = rep(dark_cols, each = 3),
                   colMed = rep(dark_cols, each = 3))

  axis(2, at = c(2, 6), labels = c("Arima", "Horseshoe"),
       las = 3, line = 3, tick = FALSE)  # Cluster labels

  # title(ylab = "Density of \nPrediction RMSE",line=6,cex.lab=1)
  title(xlab="Prediction RMSE",line=2.5,cex.lab=1)

  arima_data <- list(dat$rmse_beta[which(dat$n == 365 & dat$model == 'Seasonal \nAuto Arima')],
                     dat$rmse_beta[which(dat$n == 730 & dat$model == 'Seasonal \nAuto Arima')],
                     dat$rmse_beta[which(dat$n == 1095 & dat$model == 'Seasonal \nAuto Arima')])

  horseshoe_data <- list(dat$rmse_beta[which(dat$n == 365 & dat$model == 'Horseshoe \nPrior')],
                         dat$rmse_beta[which(dat$n == 730 & dat$model == 'Horseshoe \nPrior')],
                         dat$rmse_beta[which(dat$n == 1095 & dat$model == 'Horseshoe \nPrior')])

  # Combine the data for vioplot
  vioplot::vioplot(arima_data[[1]], arima_data[[2]], arima_data[[3]],
                   horseshoe_data[[1]], horseshoe_data[[2]], horseshoe_data[[3]],
                   at = x_positions,
                   xlab = "", ylab = "", horizontal = TRUE, las = 1,
                   names = rep(c('1-year', '2-year', '3-year'), 2),
                   col = c(modify_alpha(mod_cols, 0.4)[1],
                           modify_alpha(mod_cols, 0.8)[1],
                           mod_cols[1],
                           modify_alpha(mod_cols, 0.4)[2],
                           modify_alpha(mod_cols, 0.8)[2],
                           mod_cols[2]),
                   pchMed = 20,
                   border = rep(dark_cols, each = 3),
                   rectCol = rep(dark_cols, each = 3),
                   lineCol = rep(dark_cols, each = 3),
                   colMed = rep(dark_cols, each = 3))


  axis(2, at = c(2, 6), labels = c("Arima", "Horseshoe"),
       las = 3, line = 3, tick = FALSE)  # Cluster labels

  # title(ylab="Density of \nCovariate Effect RMSE",line=6,cex.lab=1)
  title(xlab="Covariate Effect RMSE",line= 2.5, cex.lab=1)

dev.off()

png(file = "Manuscript/Figures/Fourier_seasonality_Arima_comparison_with_Gauss.png",
    width = 5, height = 7.5, units = 'in', res = 300)

  par(mar=c(4,5,1,2))
  par(mfrow=c(2,1))
  # Define the data for each group
  arima_data <- list(gauss_dat$rmse_forecast[which(gauss_dat$n == 365 & gauss_dat$model == 'arima')],
                     gauss_dat$rmse_forecast[which(gauss_dat$n == 730 & gauss_dat$model == 'arima')],
                     gauss_dat$rmse_forecast[which(gauss_dat$n == 1095 & gauss_dat$model == 'arima')])

  gauss_data <- list(gauss_dat$rmse_forecast[which(gauss_dat$n == 365 & gauss_dat$model == 'gauss')],
                     gauss_dat$rmse_forecast[which(gauss_dat$n == 730 & gauss_dat$model == 'gauss')],
                     gauss_dat$rmse_forecast[which(gauss_dat$n == 1095 & gauss_dat$model == 'gauss')])

  horseshoe_data <- list(gauss_dat$rmse_forecast[which(gauss_dat$n == 365 & gauss_dat$model == 'hs')],
                         gauss_dat$rmse_forecast[which(gauss_dat$n == 730 & gauss_dat$model == 'hs')],
                         gauss_dat$rmse_forecast[which(gauss_dat$n == 1095 & gauss_dat$model == 'hs')])

  # Manually specify the positions for the violins
  x_positions <- c(1, 2, 3, 5, 6, 7, 9, 10, 11)

  # Combine the data for vioplot
  vioplot::vioplot(arima_data[[1]], arima_data[[2]], arima_data[[3]],
                   gauss_data[[1]], gauss_data[[2]], gauss_data[[3]],
                   horseshoe_data[[1]], horseshoe_data[[2]], horseshoe_data[[3]],
                   at = x_positions,
                   xlab = "", ylab = "", horizontal = TRUE, las = 1,
                   names = rep(c('1-year', '2-year', '3-year'), 3),
                   col = c(modify_alpha(mod_cols, 0.4)[1],
                           modify_alpha(mod_cols, 0.8)[1],
                           mod_cols[1],
                           modify_alpha('#A9A9A9', 0.4),
                           modify_alpha('#A9A9A9', 0.8),
                           '#A9A9A9',
                           modify_alpha(mod_cols, 0.4)[2],
                           modify_alpha(mod_cols, 0.8)[2],
                           mod_cols[2]),
                   pchMed = 20,
                   border = rep(c(dark_cols[1], '#000000', dark_cols[2]), each = 3),
                   rectCol = rep(c(dark_cols[1], '#000000', dark_cols[2]), each = 3),
                   lineCol = rep(c(dark_cols[1], '#000000', dark_cols[2]), each = 3),
                   colMed = rep(c(dark_cols[1], '#000000', dark_cols[2]), each = 3))

  axis(2, at = c(2, 6, 10), labels = c("Arima", "Gaussian", "Horseshoe"),
       las = 3, line = 3, tick = FALSE)  # Cluster labels

  title(ylab = "Time series length",line=6,cex.lab=1)
  title(xlab="Prediction RMSE",line=2.5,cex.lab=1)

  arima_data <- list(gauss_dat$rmse_forecast[which(gauss_dat$betas == 0 & gauss_dat$model == 'arima')],
                     gauss_dat$rmse_forecast[which(gauss_dat$betas == 2 & gauss_dat$model == 'arima')],
                     gauss_dat$rmse_forecast[which(gauss_dat$betas == 4 & gauss_dat$model == 'arima')])

  gauss_data <- list(gauss_dat$rmse_forecast[which(gauss_dat$betas == 0 & gauss_dat$model == 'gauss')],
                     gauss_dat$rmse_forecast[which(gauss_dat$betas == 2 & gauss_dat$model == 'gauss')],
                     gauss_dat$rmse_forecast[which(gauss_dat$betas == 4 & gauss_dat$model == 'gauss')])

  horseshoe_data <- list(gauss_dat$rmse_forecast[which(gauss_dat$betas == 0 & gauss_dat$model == 'hs')],
                         gauss_dat$rmse_forecast[which(gauss_dat$betas == 2 & gauss_dat$model == 'hs')],
                         gauss_dat$rmse_forecast[which(gauss_dat$betas == 4 & gauss_dat$model == 'hs')])

  # Combine the data for vioplot
  vioplot::vioplot(arima_data[[1]], arima_data[[2]], arima_data[[3]],
                   gauss_data[[1]], gauss_data[[2]], gauss_data[[3]],
                   horseshoe_data[[1]], horseshoe_data[[2]], horseshoe_data[[3]],
                   at = x_positions,
                   xlab = "", ylab = "", horizontal = TRUE, las = 1,
                   names = rep(c('0 of 5', '2 of 5', '4 of 5'), 3),
                   col = c(modify_alpha(mod_cols, 0.4)[1],
                           modify_alpha(mod_cols, 0.8)[1],
                           mod_cols[1],
                           modify_alpha('#A9A9A9', 0.4),
                           modify_alpha('#A9A9A9', 0.8),
                           '#A9A9A9',
                           modify_alpha(mod_cols, 0.4)[2],
                           modify_alpha(mod_cols, 0.8)[2],
                           mod_cols[2]),
                   pchMed = 20,
                   border = rep(c(dark_cols[1], '#000000', dark_cols[2]), each = 3),
                   rectCol = rep(c(dark_cols[1], '#000000', dark_cols[2]), each = 3),
                   lineCol = rep(c(dark_cols[1], '#000000', dark_cols[2]), each = 3),
                   colMed = rep(c(dark_cols[1], '#000000', dark_cols[2]), each = 3))


  axis(2, at = c(2, 6, 10), labels = c("Arima", "Gaussian", "Horseshoe"),
       las = 3, line = 3, tick = FALSE)  # Cluster labels

  title(ylab="Number of known covariates",line=6,cex.lab=1)
  title(xlab="Prediction RMSE",line= 2.5, cex.lab=1)

dev.off()

png(file = "Manuscript/Figures/Fourier_seasonality_Arima_comparison.png",
    width = 10, height = 6, units = 'in', res = 300)

  par(mar=c(4,7,1,2),
      oma = c(0,1,2,0))
  par(mfrow=c(2,2))
  # Define the data for each group
  arima_data <- list(dat$rmse_forecast[which(dat$n == 365 & dat$model == 'Seasonal \nAuto Arima')],
                     dat$rmse_forecast[which(dat$n == 730 & dat$model == 'Seasonal \nAuto Arima')],
                     dat$rmse_forecast[which(dat$n == 1095 & dat$model == 'Seasonal \nAuto Arima')])

  horseshoe_data <- list(dat$rmse_forecast[which(dat$n == 365 & dat$model == 'Horseshoe \nPrior')],
                         dat$rmse_forecast[which(dat$n == 730 & dat$model == 'Horseshoe \nPrior')],
                         dat$rmse_forecast[which(dat$n == 1095 & dat$model == 'Horseshoe \nPrior')])

  # Manually specify the positions for the violins
  x_positions <- c(1, 2, 3, 5, 6, 7)

  # Combine the data for vioplot
  vioplot::vioplot(arima_data[[1]], arima_data[[2]], arima_data[[3]],
                   horseshoe_data[[1]], horseshoe_data[[2]], horseshoe_data[[3]],
                   at = x_positions,
                   xlab = "", ylab = "", horizontal = TRUE, las = 1,
                   names = rep(c('1-year', '2-year', '3-year'), 2),
                   col = c(modify_alpha(mod_cols, 0.4)[1],
                           modify_alpha(mod_cols, 0.8)[1],
                           mod_cols[1],
                           modify_alpha(mod_cols, 0.4)[2],
                           modify_alpha(mod_cols, 0.8)[2],
                           mod_cols[2]),
                   pchMed = 20,
                   border = rep(dark_cols, each = 3),
                   rectCol = rep(dark_cols, each = 3),
                   lineCol = rep(dark_cols, each = 3),
                   colMed = rep(dark_cols, each = 3))

  axis(2, at = c(2, 6), labels = c("Arima", "Horseshoe"),
       las = 3, line = 3, tick = FALSE)  # Cluster labels

  title(ylab="Time series length",line=6,cex.lab=1)
  title(xlab = "Prediction RMSE",line=2.5,cex.lab=1)
  # title(main = "Time series length", line = 0, cex.lab = 1, outer = TRUE, adj = 0.25)

  arima_data <- list(dat$rmse_forecast[which(dat$betas == 0 & dat$model == 'Seasonal \nAuto Arima')],
                     dat$rmse_forecast[which(dat$betas == 2 & dat$model == 'Seasonal \nAuto Arima')],
                     dat$rmse_forecast[which(dat$betas == 4 & dat$model == 'Seasonal \nAuto Arima')])

  horseshoe_data <- list(dat$rmse_forecast[which(dat$betas == 0 & dat$model == 'Horseshoe \nPrior')],
                         dat$rmse_forecast[which(dat$betas == 2 & dat$model == 'Horseshoe \nPrior')],
                         dat$rmse_forecast[which(dat$betas == 4 & dat$model == 'Horseshoe \nPrior')])

  # Combine the data for vioplot
  vioplot::vioplot(arima_data[[1]], arima_data[[2]], arima_data[[3]],
                   horseshoe_data[[1]], horseshoe_data[[2]], horseshoe_data[[3]],
                   at = x_positions,
                   xlab = "", ylab = "", horizontal = TRUE, las = 1,
                   names = rep(c('0 of 5', '2 of 5', '4 of 5'), 2),
                   col = c(modify_alpha(mod_cols, 0.4)[1],
                           modify_alpha(mod_cols, 0.8)[1],
                           mod_cols[1],
                           modify_alpha(mod_cols, 0.4)[2],
                           modify_alpha(mod_cols, 0.8)[2],
                           mod_cols[2]),
                   pchMed = 20,
                   border = rep(dark_cols, each = 3),
                   rectCol = rep(dark_cols, each = 3),
                   lineCol = rep(dark_cols, each = 3),
                   colMed = rep(dark_cols, each = 3))

  axis(2, at = c(2, 6), labels = c("Arima", "Horseshoe"),
       las = 3, line = 3, tick = FALSE)  # Cluster labels

  title(xlab="Prediction RMSE",line=2.5,cex.lab=1)
  title(ylab="Number of known covariates",line=6,cex.lab=1)

  arima_data <- list(dat$rmse_beta[which(dat$n == 365 & dat$model == 'Seasonal \nAuto Arima')],
                     dat$rmse_beta[which(dat$n == 730 & dat$model == 'Seasonal \nAuto Arima')],
                     dat$rmse_beta[which(dat$n == 1095 & dat$model == 'Seasonal \nAuto Arima')])

  horseshoe_data <- list(dat$rmse_beta[which(dat$n == 365 & dat$model == 'Horseshoe \nPrior')],
                         dat$rmse_beta[which(dat$n == 730 & dat$model == 'Horseshoe \nPrior')],
                         dat$rmse_beta[which(dat$n == 1095 & dat$model == 'Horseshoe \nPrior')])

  # Combine the data for vioplot
  vioplot::vioplot(arima_data[[1]], arima_data[[2]], arima_data[[3]],
                   horseshoe_data[[1]], horseshoe_data[[2]], horseshoe_data[[3]],
                   at = x_positions,
                   xlab = "", ylab = "", horizontal = TRUE, las = 1,
                   names = rep(c('1-year', '2-year', '3-year'), 2),
                   col = c(modify_alpha(mod_cols, 0.4)[1],
                           modify_alpha(mod_cols, 0.8)[1],
                           mod_cols[1],
                           modify_alpha(mod_cols, 0.4)[2],
                           modify_alpha(mod_cols, 0.8)[2],
                           mod_cols[2]),
                   pchMed = 20,
                   border = rep(dark_cols, each = 3),
                   rectCol = rep(dark_cols, each = 3),
                   lineCol = rep(dark_cols, each = 3),
                   colMed = rep(dark_cols, each = 3))


  axis(2, at = c(2, 6), labels = c("Arima", "Horseshoe"),
       las = 3, line = 3, tick = FALSE)  # Cluster labels


  # title(main = "Number of known covariates", line = 0, cex.lab = 1, outer = TRUE, adj = 0.86)
  title(ylab="Time Series Length",line=6,cex.lab=1)
  title(xlab="Covariate Effect RMSE",line= 2.5, cex.lab=1)

  arima_data <- list(dat$rmse_beta[which(dat$betas == 0 & dat$model == 'Seasonal \nAuto Arima')],
                     dat$rmse_beta[which(dat$betas == 2 & dat$model == 'Seasonal \nAuto Arima')],
                     dat$rmse_beta[which(dat$betas == 4 & dat$model == 'Seasonal \nAuto Arima')])

  horseshoe_data <- list(dat$rmse_beta[which(dat$betas == 0 & dat$model == 'Horseshoe \nPrior')],
                         dat$rmse_beta[which(dat$betas == 2 & dat$model == 'Horseshoe \nPrior')],
                         dat$rmse_beta[which(dat$betas == 4 & dat$model == 'Horseshoe \nPrior')])

  # Combine the data for vioplot
  vioplot::vioplot(arima_data[[1]], arima_data[[2]], arima_data[[3]],
                   horseshoe_data[[1]], horseshoe_data[[2]], horseshoe_data[[3]],
                   at = x_positions,
                   xlab = "", ylab = "", horizontal = TRUE, las = 1,
                   names = rep(c('0 of 5', '2 of 5', '4 of 5'), 2),
                   col = c(modify_alpha(mod_cols, 0.4)[1],
                           modify_alpha(mod_cols, 0.8)[1],
                           mod_cols[1],
                           modify_alpha(mod_cols, 0.4)[2],
                           modify_alpha(mod_cols, 0.8)[2],
                           mod_cols[2]),
                   pchMed = 20,
                   border = rep(dark_cols, each = 3),
                   rectCol = rep(dark_cols, each = 3),
                   lineCol = rep(dark_cols, each = 3),
                   colMed = rep(dark_cols, each = 3))

  axis(2, at = c(2, 6), labels = c("Arima", "Horseshoe"),
       las = 3, line = 3, tick = FALSE)  # Cluster labels

  title(xlab="Covariate Effect RMSE",line=2.5,cex.lab=1)
  title(ylab="Number of known covariates",line=6,cex.lab=1)


dev.off()

png(file = "Manuscript/Figures/Fourier_seasonality_Arima_comparison.png",
    width = 10, height = 6, units = 'in', res = 300)

  par(mar=c(4,7,1,2),
      oma = c(0,1,2,0))
  par(mfrow=c(2,2))
  # Define the data for each group
  arima_data <- list(dat$rmse_forecast[which(dat$n == 365 & dat$model == 'Seasonal \nAuto Arima')],
                     dat$rmse_forecast[which(dat$n == 730 & dat$model == 'Seasonal \nAuto Arima')],
                     dat$rmse_forecast[which(dat$n == 1095 & dat$model == 'Seasonal \nAuto Arima')])

  horseshoe_data <- list(dat$rmse_forecast[which(dat$n == 365 & dat$model == 'Horseshoe \nPrior')],
                         dat$rmse_forecast[which(dat$n == 730 & dat$model == 'Horseshoe \nPrior')],
                         dat$rmse_forecast[which(dat$n == 1095 & dat$model == 'Horseshoe \nPrior')])

  # Manually specify the positions for the violins
  x_positions <- c(1, 2, 3, 5, 6, 7)

  # Combine the data for vioplot
  vioplot::vioplot(arima_data[[1]], arima_data[[2]], arima_data[[3]],
                   horseshoe_data[[1]], horseshoe_data[[2]], horseshoe_data[[3]],
                   at = x_positions,
                   xlab = "", ylab = "", horizontal = TRUE, las = 1,
                   names = rep(c('1-year', '2-year', '3-year'), 2),
                   col = c(modify_alpha(mod_cols, 0.4)[1],
                           modify_alpha(mod_cols, 0.8)[1],
                           mod_cols[1],
                           modify_alpha(mod_cols, 0.4)[2],
                           modify_alpha(mod_cols, 0.8)[2],
                           mod_cols[2]),
                   pchMed = 20,
                   border = rep(dark_cols, each = 3),
                   rectCol = rep(dark_cols, each = 3),
                   lineCol = rep(dark_cols, each = 3),
                   colMed = rep(dark_cols, each = 3))

  axis(2, at = c(2, 6), labels = c("Arima", "Horseshoe"),
       las = 3, line = 3, tick = FALSE)  # Cluster labels

  title(ylab="Time series length",line=6,cex.lab=1)
  title(xlab = "Prediction RMSE",line=2.5,cex.lab=1)
  # title(main = "Time series length", line = 0, cex.lab = 1, outer = TRUE, adj = 0.25)

  arima_data <- list(dat$rmse_forecast[which(dat$betas == 0 & dat$model == 'Seasonal \nAuto Arima')],
                     dat$rmse_forecast[which(dat$betas == 2 & dat$model == 'Seasonal \nAuto Arima')],
                     dat$rmse_forecast[which(dat$betas == 4 & dat$model == 'Seasonal \nAuto Arima')])

  horseshoe_data <- list(dat$rmse_forecast[which(dat$betas == 0 & dat$model == 'Horseshoe \nPrior')],
                         dat$rmse_forecast[which(dat$betas == 2 & dat$model == 'Horseshoe \nPrior')],
                         dat$rmse_forecast[which(dat$betas == 4 & dat$model == 'Horseshoe \nPrior')])

  # Combine the data for vioplot
  vioplot::vioplot(arima_data[[1]], arima_data[[2]], arima_data[[3]],
                   horseshoe_data[[1]], horseshoe_data[[2]], horseshoe_data[[3]],
                   at = x_positions,
                   xlab = "", ylab = "", horizontal = TRUE, las = 1,
                   names = rep(c('0 of 5', '2 of 5', '4 of 5'), 2),
                   col = c(modify_alpha(mod_cols, 0.4)[1],
                           modify_alpha(mod_cols, 0.8)[1],
                           mod_cols[1],
                           modify_alpha(mod_cols, 0.4)[2],
                           modify_alpha(mod_cols, 0.8)[2],
                           mod_cols[2]),
                   pchMed = 20,
                   border = rep(dark_cols, each = 3),
                   rectCol = rep(dark_cols, each = 3),
                   lineCol = rep(dark_cols, each = 3),
                   colMed = rep(dark_cols, each = 3))

  axis(2, at = c(2, 6), labels = c("Arima", "Horseshoe"),
       las = 3, line = 3, tick = FALSE)  # Cluster labels

  title(xlab="Prediction RMSE",line=2.5,cex.lab=1)
  title(ylab="Number of known covariates",line=6,cex.lab=1)

  arima_data <- list(dat$rmse_beta[which(dat$n == 365 & dat$model == 'Seasonal \nAuto Arima')],
                     dat$rmse_beta[which(dat$n == 730 & dat$model == 'Seasonal \nAuto Arima')],
                     dat$rmse_beta[which(dat$n == 1095 & dat$model == 'Seasonal \nAuto Arima')])

  horseshoe_data <- list(dat$rmse_beta[which(dat$n == 365 & dat$model == 'Horseshoe \nPrior')],
                         dat$rmse_beta[which(dat$n == 730 & dat$model == 'Horseshoe \nPrior')],
                         dat$rmse_beta[which(dat$n == 1095 & dat$model == 'Horseshoe \nPrior')])

  # Combine the data for vioplot
  vioplot::vioplot(arima_data[[1]], arima_data[[2]], arima_data[[3]],
                   horseshoe_data[[1]], horseshoe_data[[2]], horseshoe_data[[3]],
                   at = x_positions,
                   xlab = "", ylab = "", horizontal = TRUE, las = 1,
                   names = rep(c('1-year', '2-year', '3-year'), 2),
                   col = c(modify_alpha(mod_cols, 0.4)[1],
                           modify_alpha(mod_cols, 0.8)[1],
                           mod_cols[1],
                           modify_alpha(mod_cols, 0.4)[2],
                           modify_alpha(mod_cols, 0.8)[2],
                           mod_cols[2]),
                   pchMed = 20,
                   border = rep(dark_cols, each = 3),
                   rectCol = rep(dark_cols, each = 3),
                   lineCol = rep(dark_cols, each = 3),
                   colMed = rep(dark_cols, each = 3))


  axis(2, at = c(2, 6), labels = c("Arima", "Horseshoe"),
       las = 3, line = 3, tick = FALSE)  # Cluster labels


  # title(main = "Number of known covariates", line = 0, cex.lab = 1, outer = TRUE, adj = 0.86)
  title(ylab="Time Series Length",line=6,cex.lab=1)
  title(xlab="Covariate Effect RMSE",line= 2.5, cex.lab=1)

  arima_data <- list(dat$rmse_beta[which(dat$betas == 0 & dat$model == 'Seasonal \nAuto Arima')],
                     dat$rmse_beta[which(dat$betas == 2 & dat$model == 'Seasonal \nAuto Arima')],
                     dat$rmse_beta[which(dat$betas == 4 & dat$model == 'Seasonal \nAuto Arima')])

  horseshoe_data <- list(dat$rmse_beta[which(dat$betas == 0 & dat$model == 'Horseshoe \nPrior')],
                         dat$rmse_beta[which(dat$betas == 2 & dat$model == 'Horseshoe \nPrior')],
                         dat$rmse_beta[which(dat$betas == 4 & dat$model == 'Horseshoe \nPrior')])

  # Combine the data for vioplot
  vioplot::vioplot(arima_data[[1]], arima_data[[2]], arima_data[[3]],
                   horseshoe_data[[1]], horseshoe_data[[2]], horseshoe_data[[3]],
                   at = x_positions,
                   xlab = "", ylab = "", horizontal = TRUE, las = 1,
                   names = rep(c('0 of 5', '2 of 5', '4 of 5'), 2),
                   col = c(modify_alpha(mod_cols, 0.4)[1],
                           modify_alpha(mod_cols, 0.8)[1],
                           mod_cols[1],
                           modify_alpha(mod_cols, 0.4)[2],
                           modify_alpha(mod_cols, 0.8)[2],
                           mod_cols[2]),
                   pchMed = 20,
                   border = rep(dark_cols, each = 3),
                   rectCol = rep(dark_cols, each = 3),
                   lineCol = rep(dark_cols, each = 3),
                   colMed = rep(dark_cols, each = 3))

  axis(2, at = c(2, 6), labels = c("Arima", "Horseshoe"),
       las = 3, line = 3, tick = FALSE)  # Cluster labels

  title(xlab="Covariate Effect RMSE",line=2.5,cex.lab=1)
  title(ylab="Number of known covariates",line=6,cex.lab=1)


dev.off()

png(file = "Manuscript/Figures/Fourier_seasonality_Arima_comparison_log.png",
    width = 5, height = 5, units = 'in', res = 300)

  par(mar=c(4,9,1,2))
  par(mfrow=c(2,1))
  # Define the data for each group
  arima_data <- list(log10(dat$rmse_forecast[which(dat$n == 365 & dat$model == 'Seasonal \nAuto Arima')]),
                     log10(dat$rmse_forecast[which(dat$n == 730 & dat$model == 'Seasonal \nAuto Arima')]),
                     log10(dat$rmse_forecast[which(dat$n == 1095 & dat$model == 'Seasonal \nAuto Arima')]))

  horseshoe_data <- list(log10(dat$rmse_forecast[which(dat$n == 365 & dat$model == 'Horseshoe \nPrior')]),
                         log10(dat$rmse_forecast[which(dat$n == 730 & dat$model == 'Horseshoe \nPrior')]),
                         log10(dat$rmse_forecast[which(dat$n == 1095 & dat$model == 'Horseshoe \nPrior')]))

  # Manually specify the positions for the violins
  x_positions <- c(1, 2, 3, 5, 6, 7)

  # Combine the data for vioplot
  vioplot::vioplot(arima_data[[1]], arima_data[[2]], arima_data[[3]],
                   horseshoe_data[[1]], horseshoe_data[[2]], horseshoe_data[[3]],
                   at = x_positions, xaxt = 'n', yaxt = 'n',
                   ylim = c(0, log10(25)),
                   xlab = "", ylab = "", horizontal = TRUE, las = 1,
                   names = rep(c('1-year', '2-year', '3-year'), 2),
                   col = c(modify_alpha(mod_cols, 0.4)[1],
                           modify_alpha(mod_cols, 0.8)[1],
                           mod_cols[1],
                           modify_alpha(mod_cols, 0.4)[2],
                           modify_alpha(mod_cols, 0.8)[2],
                           mod_cols[2]),
                   pchMed = 20,
                   border = rep(dark_cols, each = 3),
                   rectCol = rep(dark_cols, each = 3),
                   lineCol = rep(dark_cols, each = 3),
                   colMed = rep(dark_cols, each = 3))
  axis(2, at = c(2, 6), labels = c("Arima", "Horseshoe"),
       las = 3, line = 3, tick = FALSE)  # Cluster labels

  axis(2, at = x_positions, labels = rep(c('1-year', '2-year', '3-year'), 2), las = 2)
  axis(1, at = log10(c(1, 2, 5, 10, 25)), labels = c(1, 2, 5, 10, 25))  # Cluster labels

  title(ylab = "Density of \nPrediction RMSE",line=6,cex.lab=1)
  title(xlab="Prediction RMSE",line=2.5,cex.lab=1)

  arima_data <- list(log10(dat$rmse_beta[which(dat$n == 365 & dat$model == 'Seasonal \nAuto Arima')]),
                     log10(dat$rmse_beta[which(dat$n == 730 & dat$model == 'Seasonal \nAuto Arima')]),
                     log10(dat$rmse_beta[which(dat$n == 1095 & dat$model == 'Seasonal \nAuto Arima')]))

  horseshoe_data <- list(log10(dat$rmse_beta[which(dat$n == 365 & dat$model == 'Horseshoe \nPrior')]),
                         log10(dat$rmse_beta[which(dat$n == 730 & dat$model == 'Horseshoe \nPrior')]),
                         log10(dat$rmse_beta[which(dat$n == 1095 & dat$model == 'Horseshoe \nPrior')]))

  # Combine the data for vioplot
  vioplot::vioplot(arima_data[[1]], arima_data[[2]], arima_data[[3]],
                   horseshoe_data[[1]], horseshoe_data[[2]], horseshoe_data[[3]],
                   at = x_positions, yaxt = 'n', xaxt = 'n',
                   # ylim = c(0, 50),
                   xlab = "", ylab = "", horizontal = TRUE, las = 1,
                   col = c(modify_alpha(mod_cols, 0.4)[1],
                           modify_alpha(mod_cols, 0.8)[1],
                           mod_cols[1],
                           modify_alpha(mod_cols, 0.4)[2],
                           modify_alpha(mod_cols, 0.8)[2],
                           mod_cols[2]),
                   pchMed = 20,
                   border = rep(dark_cols, each = 3),
                   rectCol = rep(dark_cols, each = 3),
                   lineCol = rep(dark_cols, each = 3),
                   colMed = rep(dark_cols, each = 3))


  axis(2, at = c(2, 6), labels = c("Arima", "Horseshoe"),
       las = 3, line = 3, tick = FALSE)  # Cluster labels
  axis(2, at = x_positions, labels = rep(c('1-year', '2-year', '3-year'), 2), las = 2)
  axis(1, at = seq(-2, 2), labels = c(0.01, 0.1, 1, 10, 100))  # Cluster labels

  title(ylab="Density of \nCovariate Effect RMSE",line=6,cex.lab=1)
  title(xlab="Covariate Effect RMSE",line= 2.5, cex.lab=1)

dev.off()
# plot with bars for time series length:
png('Manuscript/Figures/Fourier_seasonality_Arima_Comparison_Nbeta.png',
    width = 5, height = 10, units = 'in', res = 300)

ggplot(dat_long, aes(x = Model, y = value)) +
  geom_violin(aes(fill = Model), alpha = 0.5, width = 1, show.legend = FALSE) +
  geom_boxplot(data = dat_long[which(dat_long$betas == 0),],
               col = barcols[1], fill = barcols[1], width = 0.04,
               position = position_nudge(0.08), outlier.shape = NA) +
  geom_boxplot(data = dat_long[which(dat_long$betas == 2),],
               col = barcols[2], fill = barcols[2], width = 0.04, outlier.shape = NA) +
  geom_boxplot(data = dat_long[which(dat_long$betas == 4),],
               col = barcols[3], fill = barcols[3], width = 0.04,
               position = position_nudge(-0.08), outlier.shape = NA) +
  # ylim(0,12)+
  labs(x = "", y = "") +
  coord_flip() +
  scale_fill_manual(values = mod_cols) +
  theme_bw() +
  geom_boxplot(data = legend_data, aes(y = rmse_forecast,
                                       color = color),# fill = color),
               width = 0.1) +
  facet_wrap(metric~., ncol = 1, scales = 'free')+
  scale_color_manual(name = "Known Covariates (out of 5)",
                     values = c("grey50" = barcols[1], "grey25" = barcols[2],
                                "black" = barcols[3]),
                     labels = c("0", "2", "4")) +
  guides(color = guide_legend(override.aes = list(
    fill = c(barcols[1], barcols[2], barcols[3]),
    color = c(barcols[1], barcols[2], barcols[3]),
    size = 0.5, width = 0.1, outlier.shape = NA)))+
  geom_boxplot(data = legend_data, aes(y = rmse_forecast),
               color = 'white', fill = 'white', width = 1) +
  theme(legend.position = 'top',
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

dev.off()



ggplot(dat_long, aes(x = Model, y = rmse_forecast)) +
  geom_violin(aes(fill = Model), alpha = 0.5, width = 1.15, show.legend = FALSE) +
  geom_boxplot(data = dat_long[which(dat_long$betas == 0),],
               col = barcols[1], fill = barcols[1], width = 0.025,
               position = position_nudge(0.05), outlier.shape = NA) +
  geom_boxplot(data = dat_long[which(dat_long$betas == 2),],
               col = barcols[2], fill = barcols[2], width = 0.025, outlier.shape = NA) +
  geom_boxplot(data = dat_long[which(dat_long$betas == 4),],
               col = barcols[3], fill = barcols[3], width = 0.025,
               position = position_nudge(-0.05), outlier.shape = NA) +
  ylim(0,12)+
  labs(x = "", y = "Prediction RMSE") +
  coord_flip() +
  scale_fill_manual(values = mod_cols) +
  theme_bw() +
  geom_boxplot(data = legend_data, aes(y = rmse_forecast,
                                       color = color),# fill = color),
               width = 0.1) +
  scale_color_manual(name = "Known Covariates (out of 5)",
                     values = c("grey50" = barcols[1], "grey25" = barcols[2],
                                "black" = barcols[3]),
                     labels = c("0", "2", "4")) +
  guides(color = guide_legend(override.aes = list(
    fill = c(barcols[1], barcols[2], barcols[3]),
    color = c(barcols[1], barcols[2], barcols[3]),
    size = 0.5, width = 0.1, outlier.shape = NA)))+
  geom_boxplot(data = legend_data, aes(y = rmse_forecast),
               color = 'white', fill = 'white', width = 1) +
  theme(legend.position = 'top',
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())













ggplot(datr, aes(n, rmse_forecast)) +
  geom_violin(aes(group = cut_width(n, 5)), alpha = 0.5,
              fill = mod_cols[1]) +
  geom_violin(data = datnr, aes(group = cut_width(n, 5)), alpha = 0.5,
              fill = mod_cols[2]) +
  geom_smooth(aes(lty = siglab), col = mod_cols[1]) +
  geom_smooth(data = datnr, aes(lty = siglab), col = mod_cols[2]) +
  scale_y_log10() +
  facet_wrap(siglab~., ncol = 1)+
  ylab('Model Forecast RMSE')+
  xlab('Time series length')+
  theme_classic()

# make a dummy plot to get a legend
dummy_df <- data.frame(Model = c('Regularized', 'Not Regularized'),
                       y = c(1,2))
dummy_plot <- ggplot(dummy_df, aes(Model, y, fill = Model)) +
  geom_bar(stat = 'identity', width = 1, alpha = 0.5) +
  scale_fill_manual(values = mod_cols,
                    name = "") +
  theme(legend.direction = 'horizontal')

# lgnd <- ggpubr::get_legend(dummy_plot)
# pp <- ggpubr::ggarrange(p, legend.grob = lgnd)

# png(width = 10.5, height = 9, units = 'in', type = 'cairo', res = 300,
#     filename = 'Figures/ARp_sim_fits2.png')
#     ggpubr::annotate_figure(pp, top = ggpubr::text_grob('AR-p simulation model fits'))
# dev.off()

dat %>%
  select(-sigma, -prior) %>%
  pivot_longer(cols = starts_with(c('beta_', 'phi_')),
               values_to = 'value',
               names_to = c('parameter', 'rate'),
               names_pattern = '(beta|phi)_([a-z_]+$)') %>%
  pivot_wider(values_from = 'value', names_from = 'rate') %>%
ggplot(aes(n, true_pos, col = model)) +
  geom_point() +
  geom_line() +
  facet_grid(betas~parameter) +
  scale_color_manual(values = mod_cols) +
  theme_bw()

dd2 <- dat %>%
  pivot_longer(cols = starts_with(c('beta', 'phi')),
               values_to = 'value',
               names_to = c('parameter', 'rate'),
               names_pattern = '(beta|phi)_([a-z_]+$)') %>%
  pivot_wider(values_from = 'value', names_from = 'rate')

rates <-
  dd2 %>%
  pivot_longer(cols = ends_with('pos'),
               values_to = 'rate',
               names_to = 'category') %>%
  mutate(category = case_when(category == 'true_pos' ~ 'TPR',
                              category == 'false_pos' ~ 'FPR'),
         category = factor(category, levels = c('TPR', 'FPR'))) %>%
  group_by(prior, category, siglab, parameter, n) %>%
  summarize(rate_high = quantile(rate, 0.75, na.rm = T),
            rate_low = quantile(rate, 0.25, na.rm = T),
            rate = median(rate, na.rm = T)) %>%
  ungroup() %>%
  arrange(prior, category, siglab, parameter, n)

r <- rates %>%
  select(-rate_high, -rate_low) %>%
  tidyr::nest(data = c(rate, n)) %>%
  dplyr::mutate(
    m = purrr::map(data, loess,
                   formula = rate ~ n, span = 0.5),
    fitted = purrr::map(m, `[[`, 'fitted')
  )

r <- r %>%
  select(-m) %>%
  tidyr::unnest(cols = c(data, fitted)) %>%
  rename(rate_fit = fitted)

rh <- rates %>%
  select(-rate, -rate_low) %>%
  tidyr::nest(data = c(rate_high, n)) %>%
  dplyr::mutate(
    m = purrr::map(data, loess,
                   formula = rate_high ~ n, span = 0.5),
    fitted = purrr::map(m, `[[`, 'fitted')
  )

rh <- rh %>%
  select(-m) %>%
  tidyr::unnest(cols = c(data,fitted)) %>%
  rename(rate_high_fit = fitted)

rl <- rates %>%
  select(-rate, -rate_high) %>%
  tidyr::nest(data = c(rate_low, n)) %>%
  dplyr::mutate(
    m = purrr::map(data, loess,
                   formula = rate_low ~ n, span = 0.5),
    fitted = purrr::map(m, `[[`, 'fitted')
  )

rl <- rl %>%
  select(-m) %>%
  tidyr::unnest(cols = c(data, fitted)) %>%
  rename(rate_low_fit = fitted)

rr <- inner_join(r, rh) %>% inner_join(rl) %>%
  mutate(across(starts_with('rate'),
                ~case_when(. > 1 ~ 1,
                           . < 0 ~ 0,
                           TRUE ~ .)))

tpr <- rates %>%
  rename(Category = category) %>%
  ggplot(aes(n, rate, color = prior, lty = Category))+
  geom_ribbon(aes(ymin = rate_low, ymax = rate_high, fill = prior),
              alpha = 0.25, color = NA)+
  geom_line(size = 0.75)+
  facet_grid(parameter~siglab) +
  scale_color_manual('Prior', values = mod_cols) +
  scale_fill_manual('Prior', values = mod_cols) +
  ylab('Rate')+
  xlab('Time series length') +
  theme_classic()+
  theme(panel.border = element_rect(fill = NA),
        panel.spacing = unit(0, 'line'),
        legend.position = 'top')
# tpr <- rr %>%
#   rename(Category = category) %>%
#   ggplot(aes(n, rate_fit, color = prior, lty = Category))+
#   geom_ribbon(aes(ymin = rate_low_fit, ymax = rate_high_fit, fill = prior),
#               alpha = 0.25, color = NA)+
#   geom_line(size = 0.75)+
#   facet_grid(parameter~siglab) +
#   scale_color_manual('Prior', values = mod_cols) +
#   scale_fill_manual('Prior', values = mod_cols) +
#   ylab('Rate')+
#   xlab('Time series length') +
#   theme_classic()+
#   theme(panel.border = element_rect(fill = NA),
#         panel.spacing = unit(0, 'line'),
#         legend.position = 'top')
png('Manuscript/Figures/ARp_err_TPR_FPR.png',
    width = 6.5, height = 5, units = 'in', res = 300)
tpr
dev.off()

ggpubr::ggarrange(tpr, frmse, common.legend = TRUE,
                  nrow = 2, align = 'v', heights = c(2,1))

dat %>%
  group_by(n, model, betas) %>%
  summarize(true_se = sd(fr_true_pos),
            false_se = sd(fr_false_pos),
            true_pos = mean(fr_true_pos),
            false_pos = mean(fr_false_pos)
            ) %>%
  pivot_longer(cols = starts_with(c('true', 'false')),
               names_to = c('TF', 'val'),
               names_sep = '_', values_to = 'percent') %>%
  pivot_wider(values_from = 'percent', names_from = 'val') %>%
  mutate(Rate = case_when(TF == 'true' ~ 'true pos',
                          TF == 'false' ~ 'false pos'),
         Rate = factor(Rate, levels = c('true pos', 'false pos')))%>%
  ggplot(aes(n, pos, col = model))+
  geom_line(size = 1.2)+
  geom_ribbon(aes(ymin = pos - se, ymax = pos + se, fill = model),
              col = NA, alpha = 0.2)+
  labs(x = 'Time series length',
       y = 'Detection rate')+
  scale_color_manual('Prior', values = mod_cols) +
  scale_fill_manual('Prior', values = mod_cols) +
  facet_grid(Rate~betas, scales = 'free_y')+
  theme_classic()

dat %>%
  pivot_longer(cols = starts_with('rmse'), values_to = 'rmse',
               names_to = 'par') %>%
  filter(par != 'rmse_sigma') %>%
  ggplot(aes(factor(n), rmse, fill = model))+
  geom_boxplot(alpha = 0.5) +
  scale_fill_manual('Prior', values = mod_cols) +
  facet_wrap(.~par, ncol = 1, scale = 'free')+
  scale_y_log10()+
  labs(y = 'RMSE', x = 'timeseries length')+
  theme_classic()
dat %>%
  ggplot(aes(factor(n), rmse_phi, fill = model))+
  geom_boxplot(alpha = 0.5) +
  scale_fill_manual('Prior', values = mod_cols) +
  scale_y_log10()
dat %>%
  ggplot(aes(factor(n), rmse_forecast, fill = model))+
  geom_boxplot(alpha = 0.5) +
  scale_fill_manual('Prior', values = mod_cols) +
  # scale_y_log10()+
  labs(y = 'Forecast RMSE', x = 'timeseries length')+
  ylim(2,12)+
  theme_classic()

