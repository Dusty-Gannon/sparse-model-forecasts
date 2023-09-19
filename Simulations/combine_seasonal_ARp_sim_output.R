# concatenate files from batch SLURM runs

library(tidyverse)
# args <- 'ARp_sims_test'
devtools::load_all()
setwd("/project/modelscape/analyses/sponges/")
args <- commandArgs(trailingOnly = TRUE)
outdir <- args[1]

filelist <- paste0('Data/aquatic_sim_data/', outdir, '/',
                   list.files(paste0('Data/aquatic_sim_data/', outdir, '/'))
                   )

df <- data.frame()
for(i in 1:length(filelist)){

    #print(filelist[i])
    out <- readRDS(filelist[i])

    up_hs <- paste(out$mod_fits$hs_fit$bad_fits$par_conv$pars, collapse = ', ')
    if(is.null(up_hs)) up_hs <- NA_character_
    up_gauss <- paste(out$mod_fits$gauss_fit$bad_fits$par_conv$pars, collapse = ', ')
    if(is.null(up_gauss)) up_gauss <- NA_character_
    up_flat <- paste(out$mod_fits$flat_fit$bad_fits$par_conv$pars, collapse = ', ')
    if(is.null(up_flat)) up_flat <- NA_character_


    dd <- data.frame(
        n = rep(out$model_pars$n, 3),
        sigma = rep(out$model_pars$sd, 3),
        betas = rep(out$model_pars$beta_select, 3),
        model = c('hs', 'gauss', 'flat'),
        rmse_forecast = c(mean(out$mod_fits$hs_fit$rmse$rmse_forecast),
                          mean(out$mod_fits$gauss_fit$rmse$rmse_forecast),
                          mean(out$mod_fits$flat_fit$rmse$rmse_forecast)),
        rmse_beta = c(mean(out$mod_fits$hs_fit$rmse$rmse_beta),
                      mean(out$mod_fits$gauss_fit$rmse$rmse_beta),
                      mean(out$mod_fits$flat_fit$rmse$rmse_beta)),
        rmse_phi = c(mean(out$mod_fits$hs_fit$rmse$rmse_phi),
                      mean(out$mod_fits$gauss_fit$rmse$rmse_phi),
                      mean(out$mod_fits$flat_fit$rmse$rmse_phi)),
        rmse_sigma = c(mean(out$mod_fits$hs_fit$rmse$rmse_sigma),
                      mean(out$mod_fits$gauss_fit$rmse$rmse_sigma),
                      mean(out$mod_fits$flat_fit$rmse$rmse_sigma)),
        unconverged_pars = c(up_hs, up_gauss, up_flat),
        divergent_trans = c(out$mod_fits$hs_fit$bad_fits$divergent,
                            out$mod_fits$gauss_fit$bad_fits$divergent,
                            out$mod_fits$flat_fit$bad_fits$divergent)
    )

    tpr <- summarize_pos_rate(out$mod_fits, out$model_pars,
                              threshold = 0.9)

    dd <- bind_cols(dd, tpr)

    df <- rbind(df,dd)
}

write.csv(df, paste0('Data/aquatic_sim_data/summary_files/', outdir, '_condensed.csv'),
	  row.names = FALSE)

# df %>% select(-unconverged_pars) %>%
#   filter(divergent_trans <= 400, betas <5) %>%
#   group_by(model, n, betas) %>%
#   summarize(beta_true_pos = median(beta_true_pos, na.rm = T),
#             beta_false_pos = median(beta_false_pos, na.rm = T),
#             count = n()) %>%
#   ggplot(aes(x = n, y=beta_true_pos, col = model)) +
#   geom_point() + geom_line() +
#   geom_point(aes(y = beta_false_pos)) +
#   geom_line(aes(y = beta_false_pos), lty = 2) +
#   facet_wrap(.~factor(betas))

