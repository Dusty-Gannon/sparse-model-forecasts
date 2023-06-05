# plot figure for individual model run:

out <- readRDS('Data/aquatic_sim_data/test/simdat_run12.rds')

mod_cols <- c("#a52a2aff", "#33406fff")

beta <- out$model_pars$beta
phi <- out$model_pars$phi
beta_r <- out$mod_fit_r$par_ests$beta_hat
beta_nr <- out$mod_fit_nr$par_ests$beta_hat
phi_r <- out$mod_fit_r$par_ests$phi_hat
phi_nr <- out$mod_fit_nr$par_ests$phi_hat

forecast_r <- out$mod_fit_r$forecast[51:150,]
forecast_nr <- out$mod_fit_nr$forecast[51:150,]

ratio_beta = -max(beta_r$high, beta_nr$high, beta)/min(beta_r$low, beta_nr$low, beta)
ratio_phi = -max(phi_r$high, phi_nr$high, phi)/min(phi_r$low, phi_nr$low, phi)

ratio = min(ratio_phi, ratio_beta)


png('Manuscript/Figures/AR_param_recovery.png', width = 10, height = 3,
    units = 'in', res = 300)
  par( mar = c(2,2,1,2), oma = c(2,2,1,0))
  layout(matrix(c(1,1,2,2,2,2,2), nrow = 1, ncol = 7, byrow = TRUE))

  plot(seq(1:20), phi, pch = 2, cex = 1.2, xlab = 'Phi', ylab = 'Estimate',
       ylim = c(-max(phi_r$high, phi_nr$high, phi)/ratio,
                max(phi_r$high, phi_nr$high, phi)))
  points(seq(0.85, 19.85, by = 1), phi_r[,1],
         pch = 19, col = mod_cols[2])
  arrows(seq(0.85, 19.85, by = 1), phi_r[,2],
         seq(0.85, 19.85, by = 1), phi_r[,3],
         length = 0, col = mod_cols[2])
  points(seq(1.15, 20.15, by = 1), phi_nr[,1], pch = 19, col = mod_cols[1])
  arrows(seq(1.15, 20.15, by = 1), phi_nr[,2],
         seq(1.15, 20.15, by = 1), phi_nr[,3],
         length = 0, col = mod_cols[1])
  mtext('Estimate', side = 2, line = 2, cex = 0.9)
  mtext('Phi 1-20', side = 1, line = 2.2, cex = 0.9)
  plot(seq(1:50), beta[1:50], pch = 2, cex = 1.2, xlab = 'Beta', ylab = 'Value',
       ylim = c(-max(beta_r$high, beta_nr$high, beta)/ratio,
                max(beta_r$high, beta_nr$high, beta)))
  points(seq(0.85, 49.85, by = 1), beta_r[1:50,1], pch = 19, col = mod_cols[2])
  arrows(seq(0.85, 49.85, by = 1), beta_r[1:50,2],
         seq(0.85, 49.85, by = 1), beta_r[1:50,3],
         length = 0, col = mod_cols[2])
  points(seq(1.15, 50.15, by = 1), beta_nr[1:50,1], pch = 19, col = mod_cols[1])
  arrows(seq(1.15, 50.15, by = 1), beta_nr[1:50,2],
         seq(1.15, 50.15, by = 1), beta_nr[1:50,3],
         length = 0, col = mod_cols[1])
  mtext('Beta 1-50', side = 1, line = 2.2, cex = 0.9)

  legend('topright',
         legend = c('True value', 'Horseshoe', 'Gaussian'),
         pch = c(2, 19, 19), col = c('black', mod_cols[2], mod_cols[1]),
         bty = 'n', cex = 1.2)

dev.off()





# version with forecast

par( mar = c(5,3,1,2), oma = c(0,2,1,0))
layout(matrix(c(1,1,1,1,1,1,1,2,2,3,3,3,3,3), nrow = 2, ncol = 7, byrow = TRUE))

plot(forecast_r$time, forecast_r$y, type = 'l',
     ylim = c(min(forecast_nr$low), max(forecast_nr$high)),
     ylab = 'GPP', xlab = 'Time' )
polygon(x = c(forecast_nr$time, rev(forecast_nr$time)),
        y = c(forecast_nr$low, rev(forecast_nr$high)),
        col = adjustcolor(mod_cols[1], alpha.f = 0.5), border = NA)
polygon(x = c(forecast_r$time, rev(forecast_r$time)),
        y = c(forecast_r$low, rev(forecast_r$high)),
        col = adjustcolor(mod_cols[2], alpha.f = 0.5), border = NA)
lines(forecast_r$time, forecast_r$estim, col = mod_cols[2])
lines(forecast_nr$time, forecast_r$estim, col = mod_cols[1])
abline(v = 100)
mtext('GPP', 2, line = 2)

legend(100, 30,
       legend = c('Regularized', 'Not regularized'),
       title = "Model Forecast:",
       fill = c(mod_cols[2], mod_cols[1]),
       bty = 'n')



plot(seq(1:20), phi, pch = 2, xlab = 'Phi', ylab = '')#, ylim = c(-5.5,5.8))
points(seq(0.85, 19.85, by = 1), phi_r[,1], pch = 19, col = mod_cols[2])
arrows(seq(0.85, 19.85, by = 1), phi_r[,2],
       seq(0.85, 19.85, by = 1), phi_r[,3],
       length = 0, col = mod_cols[2])
points(seq(1.15, 20.15, by = 1), phi_nr[,1], pch = 19, col = mod_cols[1])
arrows(seq(1.15, 20.15, by = 1), phi_nr[,2],
       seq(1.15, 20.15, by = 1), phi_nr[,3],
       length = 0, col = mod_cols[1])
mtext('Estimate', 2, line = 2)

plot(seq(1:51), beta, pch = 2, xlab = 'Beta', ylab = 'Value', ylim = c(-5.5,5.8))
points(seq(0.85, 50.85, by = 1), beta_r[,1], pch = 19, col = mod_cols[2])
arrows(seq(0.85, 50.85, by = 1), beta_r[,2],
       seq(0.85, 50.85, by = 1), beta_r[,3],
       length = 0, col = mod_cols[2])
points(seq(1.15, 51.15, by = 1), beta_nr[,1], pch = 19, col = mod_cols[1])
arrows(seq(1.15, 51.15, by = 1), beta_nr[,2],
       seq(1.15, 51.15, by = 1), beta_nr[,3],
       length = 0, col = mod_cols[1])

legend('topright',
       legend = c('true value', 'regularized', 'not regularized'),
       pch = c(2, 19, 19), col = c('black', mod_cols[2], mod_cols[1]),
       bty = 'n')
