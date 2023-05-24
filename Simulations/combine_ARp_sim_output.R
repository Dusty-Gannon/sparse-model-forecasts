# concatenate files from batch SLURM runs

library(tidyverse)
# args <- 'ARp_sims_test'
setwd("/project/modelscape/analyses/sponges/")
args <- commandArgs(trailingOnly = TRUE)

# filelist <- paste0('Data/aquatic_sim_data/test/',
#                    list.files(paste0('Data/aquatic_sim_data/test/'))
#                    )
filelist <- paste0('Data/aquatic_sim_data/', args[1], '/',
                   list.files(paste0('Data/aquatic_sim_data/', args[1], '/'))
                   )

#out <- lapply( filelist, function(x) c(readRDS(x)))

#saveRDS(out, paste0('Data/aquatic_sim_data/', args[1], '.rds'))

# Calculate true positives:
calculate_true_pos_rate <- function(ests, value, par = 'beta'){
    min_pos = 0
    if(par == 'beta') min_pos = 1.96 * 0.05
    ests <- ests %>%
        mutate(value = value,
             zero = case_when(low <= 0 & 0 <= high ~ TRUE,
                              TRUE ~FALSE),
             captures_est = case_when(low <= value & value <= high ~ TRUE,
                                      TRUE ~FALSE),
             estimate = case_when(abs(value) > min_pos ~ 'Positive',
                                  TRUE ~ 'Negative'))
    true_pos = sum(ests$estimate == 'Positive' & ests$captures_est)/
      sum(ests$estimate == 'Positive')
    true_neg = sum(ests$estimate == 'Negative' & ests$zero)/
      sum(ests$estimate == 'Negative')
    false_pos = sum(ests$estimate == 'Negative' & !ests$zero)/
      sum(ests$estimate == 'Negative')
    false_neg = sum(ests$estimate == 'Positive' & ests$zero)/
      sum(ests$estimate == 'Positive')

    return(data.frame(true_pos = true_pos,
                      true_neg = true_neg,
                      false_pos = false_pos,
                      false_neg = false_pos)
    )
}

summarize_pos_rate <- function(out){
    beta <- out$model_pars$beta
    phi <- out$model_pars$phi
    beta_r <- out$mod_fit_r$par_ests$beta_hat
    beta_nr <- out$mod_fit_nr$par_ests$beta_hat
    phi_r <- out$mod_fit_r$par_ests$phi_hat
    phi_nr <- out$mod_fit_nr$par_ests$phi_hat

    pos_rate_beta <- bind_rows(calculate_true_pos_rate(beta_r, beta),
                               calculate_true_pos_rate(beta_nr, beta))
    colnames(pos_rate_beta) <- paste0('beta_', colnames(pos_rate_beta))
    pos_rate_phi <- bind_rows(calculate_true_pos_rate(phi_r, phi, par = 'phi'),
                              calculate_true_pos_rate(phi_nr, phi, par = 'phi'))
    colnames(pos_rate_phi) <- paste0('phi_', colnames(pos_rate_phi))

    return(cbind(pos_rate_beta, pos_rate_phi))

}

out <- readRDS(filelist[1])
beta <- out$model_pars$beta
phi <- out$model_pars$phi
beta_r <- out$mod_fit_r$par_ests$beta_hat
beta_nr <- out$mod_fit_nr$par_ests$beta_hat
phi_r <- out$mod_fit_r$par_ests$phi_hat
phi_nr <- out$mod_fit_nr$par_ests$phi_hat

forecast_r <- out$mod_fit_r$forecast[51:150,]
forecast_nr <- out$mod_fit_nr$forecast[51:150,]
    par( mar = c(5,3,1,2), oma = c(0,2,1,0))
    layout(matrix(c(1,1,1,1,1,1,1,2,2,3,3,3,3,3), nrow = 2, ncol = 7, byrow = TRUE))

    plot(forecast_r$time, forecast_r$y, type = 'l',
         ylim = c(min(forecast_nr$low), max(forecast_nr$high)),
         ylab = 'GPP', xlab = 'Time' )
    polygon(x = c(forecast_nr$time, rev(forecast_nr$time)),
            y = c(forecast_nr$low, rev(forecast_nr$high)),
            col = adjustcolor('grey50', alpha.f = 0.5), border = NA)
    polygon(x = c(forecast_r$time, rev(forecast_r$time)),
            y = c(forecast_r$low, rev(forecast_r$high)),
            col = adjustcolor('brown3', alpha.f = 0.5), border = NA)
    lines(forecast_r$time, forecast_r$estim, col = 'brown3')
    lines(forecast_nr$time, forecast_r$estim, col = 'grey50')
    abline(v = 100)
    mtext('GPP', 2, line = 2)

    legend(100, 30,
           legend = c('Regularized', 'Not regularized'),
           title = "Model Forecast:",
           fill = c('brown3', 'grey50'),
           bty = 'n')



    plot(seq(1:20), phi, pch = 2, xlab = 'Phi', ylab = '')#, ylim = c(-5.5,5.8))
    points(seq(0.85, 19.85, by = 1), phi_r[,1], pch = 19, col = 'brown3')
    arrows(seq(0.85, 19.85, by = 1), phi_r[,2],
           seq(0.85, 19.85, by = 1), phi_r[,3],
           length = 0, col = 'brown3')
    points(seq(1.15, 20.15, by = 1), phi_nr[,1], pch = 19, col = 'grey50')
    arrows(seq(1.15, 20.15, by = 1), phi_nr[,2],
           seq(1.15, 20.15, by = 1), phi_nr[,3],
           length = 0, col = 'grey50')
    mtext('Estimate', 2, line = 2)

    plot(seq(1:51), beta, pch = 2, xlab = 'Beta', ylab = 'Value', ylim = c(-5.5,5.8))
    points(seq(0.85, 50.85, by = 1), beta_r[,1], pch = 19, col = 'brown3')
    arrows(seq(0.85, 50.85, by = 1), beta_r[,2],
           seq(0.85, 50.85, by = 1), beta_r[,3],
           length = 0, col = 'brown3')
    points(seq(1.15, 51.15, by = 1), beta_nr[,1], pch = 19, col = 'grey50')
    arrows(seq(1.15, 51.15, by = 1), beta_nr[,2],
           seq(1.15, 51.15, by = 1), beta_nr[,3],
           length = 0, col = 'grey50')

    legend('topright',
           legend = c('true value', 'regularized', 'not regularized'),
           pch = c(2, 19, 19), col = c('black', 'brown3', 'grey50'),
           bty = 'n')


df <- data.frame()
for(i in 1:length(filelist)){

    out <- readRDS(filelist[i])

    up_r <- paste(out$mod_fit_r$bad_fits$par_conv$pars, collapse = ', ')
    if(is.null(up_r)) up_r <- NA_character_
    up_nr <- paste(out$mod_fit_nr$bad_fits$par_conv$pars, collapse = ', ')
    if(is.null(up_nr)) up_nr <- NA_character_


    dd <- data.frame(
        n = rep(out$model_pars$n, 2),
        sigma = rep(out$model_pars$sigma_e, 2),
        model = c('reg', 'not_reg'),
        rmse_forecast = c(mean(out$mod_fit_r$rmse$rmse_forecast),
                          mean(out$mod_fit_nr$rmse$rmse_forecast)),
        rmse_beta = c(mean(out$mod_fit_r$rmse$rmse_beta),
                      mean(out$mod_fit_nr$rmse$rmse_beta)),
        rmse_phi = c(mean(out$mod_fit_r$rmse$rmse_phi),
                      mean(out$mod_fit_nr$rmse$rmse_phi)),
        rmse_sigma = c(mean(out$mod_fit_r$rmse$rmse_sigma),
                      mean(out$mod_fit_nr$rmse$rmse_sigma)),
        unconverged_pars = c(up_r, up_nr),
        divergent_trans = c(out$mod_fit_r$bad_fits$divergent,
                            out$mod_fit_nr$bad_fits$divergent)
    )

    dd <- cbind(dd, summarize_pos_rate(out))


    df <- rbind(df,dd)
}

write.csv(df, paste0('Data/aquatic_sim_data/summary_files/', args[1], '_condensed.csv'),
	  row.names = FALSE)
