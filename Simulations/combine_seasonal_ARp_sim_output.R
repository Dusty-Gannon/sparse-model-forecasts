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

    up_hs <- paste(out$fits$hs_fit$bad_fits$par_conv$pars, collapse = ', ')
    if(is.null(up_hs)) up_hs <- NA_character_
    up_gauss <- paste(out$fits$gauss_fit$bad_fits$par_conv$pars, collapse = ', ')
    if(is.null(up_gauss)) up_gauss <- NA_character_
    up_flat <- paste(out$fits$flat_fit$bad_fits$par_conv$pars, collapse = ', ')
    if(is.null(up_flat)) up_flat <- NA_character_


    dd <- data.frame(
        n = rep(out$model_pars$n, 3),
        sigma = rep(out$model_pars$sd, 3),
        model = c('hs', 'gauss', 'flat'),
        rmse_forecast = c(mean(out$fits$hs_fit$rmse$rmse_forecast),
                          mean(out$fits$gauss_fit$rmse$rmse_forecast),
                          mean(out$fits$flat_fit$rmse$rmse_forecast)),
        rmse_beta = c(mean(out$fits$hs_fit$rmse$rmse_beta),
                      mean(out$fits$gauss_fit$rmse$rmse_beta),
                      mean(out$fits$flat_fit$rmse$rmse_beta)),
        rmse_phi = c(mean(out$fits$hs_fit$rmse$rmse_phi),
                      mean(out$fits$gauss_fit$rmse$rmse_phi),
                      mean(out$fits$flat_fit$rmse$rmse_phi)),
        rmse_sigma = c(mean(out$fits$hs_fit$rmse$rmse_sigma),
                      mean(out$fits$gauss_fit$rmse$rmse_sigma),
                      mean(out$fits$flat_fit$rmse$rmse_sigma)),
        unconverged_pars = c(up_hs, up_gauss, up_flat),
        divergent_trans = c(out$fits$hs_fit$bad_fits$divergent,
                            out$fits$gauss_fit$bad_fits$divergent,
                            out$fits$flat_fit$bad_fits$divergent)
    )

    tpr <- summarize_pos_rate(out$fits, out$model_pars,
                              threshold = 0.9)

    dd <- left_join(dd, tpr, by = 'model')

    df <- rbind(df,dd)
}

write.csv(df, paste0('Data/aquatic_sim_data/summary_files/', outdir, '_condensed.csv'),
	  row.names = FALSE)


