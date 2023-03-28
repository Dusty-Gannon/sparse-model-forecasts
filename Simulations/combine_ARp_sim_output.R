# concatenate files from batch SLURM runs

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


df <- data.frame()
for(i in 1:length(filelist)){

 out <- readRDS(filelist[i])

 up_r <- paste(out$mod_fit_r$bad_fits$par_cov$pars)
 if(is.null(up_r)) up_r <- NA_character_
 up_nr <- c(out$mod_fit_nr$bad_fits$par_cov$pars)
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



  df <- rbind(df,dd)
}

write.csv(df, paste0('Data/aquatic_sim_data/summary_files/', args[1], '_condensed.csv'),
	  row.names = FALSE)
