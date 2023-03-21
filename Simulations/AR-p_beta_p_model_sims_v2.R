
# libraries
library(rstan)
#library(here)
devtools::load_all()
rstan_options(auto_write = TRUE)
setwd("/project/modelscape/analyses/sponges/")
# setwd("C:/Users/alice.carter/git/sponges/")

#### Simulate data and run regularized and non-regularized AR-p_beta models ####
# store data inputs and model fits in a log file
   nsims <- 300
   sigmas <- c(0.1, 1, 10)
   #lengths <- 100
   lengths <- c(65, 125, 185, 245, 370)

   sim_df <- data.frame(
     sigma = rep(sigmas, each = nsims * length(lengths)),
     length = rep(rep(lengths, each = nsims), length(sigmas))
   )
  
  ARp_beta_sims <- function(input_pars){

    model_pars <- simulate_AR_p_beta_p_timeseries(input_pars, draw_beta = 'zero', draw_phi = 'near_zero')
    fits <- fit_ARp_beta_model(model_pars, iter = 4000)
    sim_list <- unpack_ARp_fit(fits)
    return(sim_list)

  }


  args <- commandArgs(trailingOnly = TRUE)

  #number of simulations to run:
  i <- as.numeric(args[2])
  i <- 30*(i-1) + 1

  for(j in 0:29){

  nsteps <- sim_df$length[i+j]
  sigma <- sim_df$sigma[i+j]

  # add new lines to error and out file identifying which model this is
  write(paste('model number = ', i+j, ', nsteps = ', nsteps, ', sigma = ', sigma),
        stderr())
  
# simulate AR-p data (will update in future to have variable inputs)
  input_pars <- list(
    n = nsteps,   # length of time series
    p  = 16,      # number of AR lags to consider
    beta_p = 5,   # number of beta lags to consider in lagged covariate (beta_1)
    beta_n = 45,  # number of additional covariates to include
    b0 = 0.5,     # intercept
    n_phi = 3,    # number of non-zero autoregressive terms
    n_lags = 2,   # number of non-zero lags in covariate
    n_beta = 2,   # number of non-zero covariate parameters
    non_zero_coef_guess = 5, # guess for the number of non-zero coefficients
    sigma_e = sigma,  #standard deviation of the innovations
    holdout = 30
  )


  sim_dat <- ARp_beta_sims(input_pars)

  # Save the results
  fname <- paste0("/simdat_", i+j, ".rds")
  fpath <- paste0("Data/aquatic_sim_data/", args[1], fname)
  saveRDS(sim_dat, file = paste0("/project/modelscape/analyses/sponges/", fpath))

  }
