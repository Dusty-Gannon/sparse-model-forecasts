# concatenate files from batch SLURM runs

library(tidyverse)
# args <- 'ARp_sims_test'
devtools::load_all()
setwd("/project/modelscape/analyses/sponges/")
args <- commandArgs(trailingOnly = TRUE)
# outdir <- 'test_seasonal'
outdir <- args[1]

filelist <- paste0('Data/aquatic_sim_data/', outdir, '/',
                   list.files(paste0('Data/aquatic_sim_data/', outdir, '/'))
                   )
# df <- read_csv('Data/aquatic_sim_data/ARp_err_sims_9_18_condensed.csv')
df <- data.frame()
for(i in 1:length(filelist)){

    #print(filelist[i])
    out <- readRDS(filelist[i])

    up_hs <- paste(out$mod_fits$hs_fit$bad_fits$par_conv$pars, collapse = ', ')
    if(is.null(up_hs)) up_hs <- NA_character_
    up_gauss <- paste(out$mod_fits$gauss_fit$bad_fits$par_conv$pars, collapse = ', ')
    if(is.null(up_gauss)) up_gauss <- NA_character_
    # up_flat <- paste(out$mod_fits$flat_fit$bad_fits$par_conv$pars, collapse = ', ')
    # if(is.null(up_flat)) up_flat <- NA_character_


    dd <- data.frame(
        n = rep(out$model_pars$n, 2),
        sigma = rep(out$model_pars$sd, 2),
        betas = rep(out$model_pars$beta_select, 2),
        model = c('hs', 'gauss'),
        rmse_forecast = c(mean(out$mod_fits$hs_fit$rmse$rmse_forecast),
                          mean(out$mod_fits$gauss_fit$rmse$rmse_forecast)),
        rmse_beta = c(mean(out$mod_fits$hs_fit$rmse$rmse_beta),
                      mean(out$mod_fits$gauss_fit$rmse$rmse_beta)),
        rmse_phi = c(mean(out$mod_fits$hs_fit$rmse$rmse_phi),
                      mean(out$mod_fits$gauss_fit$rmse$rmse_phi)),
        rmse_sigma = c(mean(out$mod_fits$hs_fit$rmse$rmse_sigma),
                      mean(out$mod_fits$gauss_fit$rmse$rmse_sigma)),
        # unconverged_pars = c(up_hs, up_gauss),
        divergent_trans = c(out$mod_fits$hs_fit$bad_fits$divergent,
                            out$mod_fits$gauss_fit$bad_fits$divergent)
    )

    tpr <- summarize_pos_rate(out$mod_fits, out$model_pars,
                              threshold = 0.9, fr = TRUE) %>%
      select(-model)

    dd <- bind_cols(dd, tpr)
    dd$mod_run <- rep(filelist[i], 2)

    print(paste0('simulation ', filelist[i]))
    arima_fit <- fit_arima_model(out$model_pars)
    arima_seasonal <- fit_seasonal_arima_model(out$model_pars)

    ar <- data.frame(
      n = out$model_pars$n,
      sigma = out$model_pars$sd,
      betas = out$model_pars$beta_select,
      model = 'arima',
      rmse_forecast = arima_seasonal$rmse)

    tpr_arima <- summarize_arima_pos_rate(arima_seasonal$fit_ar,
                                          out$mod_fits$hs_fit,
                                          out$model_pars)

    ar <- bind_cols(ar, tpr_arima)

    dd <- bind_rows(dd, ar)
    # ar_for <- arima_seasonal$ar_forecast

    df <- rbind(df,dd)

    # hs_for <- t(apply(out$mod_fits$hs_fit$forecast, 2, quantile,
    #                   probs = c(0.025, 0.5, 0.975)))
    # colnames(hs_for) <- c('hs_lower', 'hs_med', 'hs_upper')
    # ar_for <- arima_fit$ar_forecast
    # ar_for <- data.frame(ar_lower = c(ar_for$fitted, as.numeric(ar_for$lower[,2])),
    #                      ar_med = c(ar_for$fitted, ar_for$mean),
    #                      ar_upper = c(ar_for$fitted, as.numeric(ar_for$upper[,2])))
    # fcst <- data.frame(y = as.numeric(out$model_pars$y)) %>%
    #   bind_cols(hs_for, ar_for)
    #
    #
    # fcst[(nrow(fcst)-199):nrow(fcst),] %>%
    #   mutate(time = 1:200) %>%
    #   pivot_longer(cols = starts_with(c('hs', 'ar')),
    #                names_to = c('model', 'stat'), values_to = 'val',
    #                names_sep = '_') %>%
    #   pivot_wider(names_from = 'stat', values_from = 'val') %>%
    # ggplot(aes(time, y)) +
    #   geom_ribbon(aes(ymin = lower, ymax = upper, fill = model), alpha = 0.3) +
    #   geom_line(aes(y = med, col = model)) +
    #   geom_line()+
    #   geom_vline(xintercept = 100, lty = 2) +
    #   theme_bw()

}

write.csv(df, paste0('Data/aquatic_sim_data/summary_files/', outdir, '_condensed.csv'),
	  row.names = FALSE)

#glimpse(df)
#mod_cols <- c('grey', "#33406fff", "#a52a2aff")

#png('Figures/AR_fourier_model_convergence.png')
#df %>%# select(-unconverged_pars) %>%
#  mutate(converged = case_when(divergent_trans < 20 ~ 1,
#                               TRUE ~ 0),
#         Prior = case_when(model == 'flat' ~ 'Flat',
#                           model == 'gauss' ~ 'Gaussian',
#                           model == 'hs' ~ 'Horseshoe')) %>%
#  group_by(n, Prior) %>%
#  summarize(converged = mean(converged),
#            count = n()) %>%
#  ggplot(aes(n/365, converged, col = Prior)) +
#  geom_line(size = 1) +
#  theme_classic() + xlab('Time series length (y)') + ylab('Fraction of models that converged')
#
#dev.off()
#  mutate(percent_conv = 1 - divergent_trans/2000) %>%
#  ggplot(aes(factor(n), percent_conv, fill = model)) +
#  geom_boxplot()
#df %>% select(-unconverged_pars) %>%
#  filter(divergent_trans <= 20, betas <5) %>%
#  group_by(model, n, sigma, betas) %>%
#  summarize(beta_true_pos = median(beta_true_pos, na.rm = T),
#            rmse_forecast = median(rmse_forecast, na.rm = T),
#            count = n()) %>%
#  ggplot(aes(x = n/365, y=log(rmse_forecast), col = model, lty = factor(sigma))) +
#  geom_point() + geom_line() +
#  facet_wrap(.~factor(betas)) + theme_classic() +
#  xlab('Time series length (y)')
#
#png('Figures/seasonal_AR_test_comparisonA.png', width = 6.5, height = 3, units = 'in', res = 300)
#df %>% select(-unconverged_pars) %>%
#  filter(divergent_trans <= 400, betas <5) %>%
#  ggplot(aes(x = factor(betas), y=log(rmse_forecast), fill = model)) +
#  geom_boxplot() + xlab('Number of observed covariates')+
#  theme_classic()
#dev.off()
#
#df %>% select(-unconverged_pars) %>%
#  filter(divergent_trans <= 400, betas <5, model != 'flat') %>%
#  group_by(model, n, sigma, betas) %>%
#  # group_by(model, n, betas) %>%
#  summarize(beta_true_pos = median(beta_true_pos, na.rm = T),
#            rmse_forecast = median(rmse_forecast, na.rm = T),
#            count = n()) %>%
#  # ggplot(aes(x = n/365, y=rmse_forecast, col = model, lty = factor(sigma))) +
#  ggplot(aes(x = n/365, y=rmse_forecast, col = model)) +
#  geom_point() + geom_line() +
#  facet_wrap(.~factor(betas)) + theme_classic() +
#  xlab('Time series length (y)')
#  # xlab('sd of random innovations')
#
#df %>% select(-unconverged_pars) %>%
#  filter(divergent_trans <= 400, betas <5, model != 'flat') %>%
#  mutate(ts_length = as.factor(n/365),
#         sigma = as.factor(sigma),
#         betas = as.factor(betas)) %>%
#  pivot_longer(cols = c('model', 'sigma', 'betas'),
#               names_to = 'variable', values_to = 'value') %>%
#  ggplot(aes(x = ts_length, y = rmse_forecast, col = model)) +
#  geom_point() + geom_line() +
#  facet_wrap(.~factor(betas)) + theme_classic() +
#  xlab('Time series length (y)')
#  # xlab('sd of random innovations')
#
#png('Figures/seasonal_AR_test_comparisonB.png', width = 8, height = 8, units = 'in', res = 300)
#  ggpubr::ggarrange(A,B,C, nrow = 3, common.legend = T)
#dev.off()
