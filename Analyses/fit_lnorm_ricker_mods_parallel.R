#################################################################
# Fitting Ricker models with log-normal demographic stochasticity
#################################################################


##### Doc setup #####

  # libraries
  library(rstan)
  library(here)
  library(parallel)
  devtools::load_all()


##################################
# Fit models and summarize results
##################################

  # command line arguments
  args <- commandArgs(trailingOnly = T)
  datfile <- args[1]
  tsteps <- as.numeric(args[2]):as.numeric(args[3])
  sim_ids <- as.numeric(args[4]):as.numeric(args[5])

  # test if disturbance is specified
  if(length(args) > 6 & args[7] == "dist"){
    dist <- T
  } else{
    dist <- F
  }

  # print some things for the log file
  print(
    paste0("Starting on sims ", sim_ids[1], " to ", sim_ids[length(sim_ids)])
  )

  # load simulated data
  datfp <- paste0("Data/terrestrial_sim_data/lnorm_ricker/", datfile)
  sims <- readRDS(here::here(datfp))[sim_ids]

  # compile stan model
  growth_mod <- rstan::stan_model(
    here::here("Stan/pop_growth_rate_FHS.stan")
  )

  # apply the above wrapper to each simulated dataset
  cl <- makeCluster(length(sim_ids))
    clusterEvalQ(cl, expr = {
      library(here); library(rstan); devtools::load_all()
      # src_files <- list.files(here("R/"), pattern = "*.R", full.names = T);
      # sapply(src_files, source, .GlobalEnv);
    })

    results <- parLapply(
      cl,
      sims,
      fun = fit_n_summarize,
      stan_mod = growth_mod,
      tsteps = tsteps,
      pip = 0.9,
      dist = dist
    )

  stopCluster(cl)

  # save the results
  outfp <- paste0("Data/terrestrial_sim_data/lnorm_ricker/", args[6])
  saveRDS(results, file = here(outfp))








