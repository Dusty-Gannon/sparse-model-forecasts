# concatenate files from batch SLURM runs

# args <- 'ARp_sims_test'
setwd("/project/modelscape/analyses/sponges/")
args <- commandArgs(trailingOnly = TRUE)

filelist <- paste0('Data/aquatic_sim_data/', args[1], '/',
                   list.files(paste0('Data/aquatic_sim_data/', args[1], '/'))
                   )

out <- lapply( filelist, function(x) c(readRDS(x)))

saveRDS(out, paste0('Data/aquatic_sim_data/', args[1], '_', args[2], 'steps_',
		    args[3], 'sigma.rds'))


df <- data.frame()
for(i in 1:length(out)){
  dd <- data.frame(
    n = rep(out[[i]]$model_pars$n, 2),
    sigma = rep(out[[i]]$model_pars$sigma_e, 2),
    model = c('reg', 'not_reg'),
    rmse_forecast = c(mean(out[[i]]$mod_fit_r$rmse$rmse_forecast),
                      mean(out[[i]]$mod_fit_nr$rmse$rmse_forecast)),
    rmse_beta = c(mean(out[[i]]$mod_fit_r$rmse$rmse_beta),
                  mean(out[[i]]$mod_fit_nr$rmse$rmse_beta)),
    rmse_phi = c(mean(out[[i]]$mod_fit_r$rmse$rmse_phi),
                  mean(out[[i]]$mod_fit_nr$rmse$rmse_phi)),
    rmse_sigma = c(mean(out[[i]]$mod_fit_r$rmse$rmse_sigma),
                  mean(out[[i]]$mod_fit_nr$rmse$rmse_sigma))
  )

  df <- rbind(df,dd)
}

write.csv(df, paste0('Data/aquatic_sim_data/', args[1], '_', args[2], 'steps_',
		     args[3], 'sigma_condensed.csv'),
	  row.names = FALSE)
