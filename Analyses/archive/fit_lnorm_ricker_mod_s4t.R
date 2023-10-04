#################################################################
# Fitting Ricker models with log-normal demographic stochasticity
# and data on replicate communities but shorter time scales
#################################################################


##### Doc setup #####

  # libraries
  library(rstan)
  library(here)
  devtools::load_all()


##################################
# Fit models and summarize results
##################################

  # command line arguments
  args <- commandArgs(trailingOnly = T)
  datfile <- args[1]
  id <- args[2]
  tsteps <- as.numeric(args[3]):as.numeric(args[4])

  # test if disturbance is specified
  if(length(args) > 5 & args[6] == "dist"){
    dist <- T
  } else{
    dist <- F
  }

  # print some things for the log file
  print(paste0(
    "Arguments: ", args
  ))
  print(
    paste0("Starting on sim result:", datfile)
  )

  # load simulated data
  sim <- readRDS(here::here(datfile))

  # compile stan model
  growth_mod <- rstan::stan_model(
    here::here("Stan/indrep_pop_growth_rate_FHS.stan")
  )

  # apply the wrapper function
  results <- fit_n_summarize(
    sim,
    stan_mod = growth_mod,
    tsteps = tsteps,
    dist = dist,
    sp = T
  )

  # add details about the input params
  results <- c(
    results,
    list(
      sim_params = sim$sim_params
    )
  )

  # save the results
  outfp <- paste0(
    args[5],
    "_", id, ".rds"
  )
  saveRDS(results, file = here(outfp))








