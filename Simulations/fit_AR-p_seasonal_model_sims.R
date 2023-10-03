
# libraries
library(rstan)
library(dplyr)
#library(here)
devtools::load_all()
#rstan_options(auto_write = TRUE)
setwd("/project/modelscape/analyses/sponges/")

#### Simulate data and run regularized and non-regularized AR-p_beta models ####
# store data inputs and model fits in a log file

# compile stan models:
ar_err_gauss <- stan_model("Stan/AR-p_err3_Gauss_DG.stan")
# ar_err_flat <- stan_model("Stan/AR-p_err3_Flat_DG.stan")
ar_err_hs <- stan_model("Stan/AR-p_err3_FHS_DG.stan")

beta <- c(0, 2, 4) # how many of the important betas did you measure?
sigmas <- c(0.5, 2, 5)
lengths <- 365 * (1:3)

s_df <- expand.grid(beta, sigmas, lengths)
colnames(s_df) <- c('beta', 'sigma', 'length')

reps <- 5
sim_df <- s_df
for(i in 1:(reps-1)){
  sim_df <- rbind(sim_df, s_df)
}

mods_per_node <- 9
array_size <- nrow(sim_df)/mods_per_node

ARp_beta_sims <- function(input_pars){

    model_pars <- simulate_seasonal_AR_p_timeseries(input_pars)
    out <- fit_seasonal_ARp_models(model_pars, fit_flat = FALSE)
    sim_list <- unpack_seasonal_ARp_fits(fits = out$fits, model_pars = out$model_pars)

    return(sim_list)

}


# outdir <- 'test_seasonal'; array_num <- 1; i = 1
args <- commandArgs(trailingOnly = TRUE)
outdir <- args[1]
array_num <- as.numeric(args[2])

for(i in 1:mods_per_node){

    mod_num = mods_per_node * (array_num - 1) + i

    nsteps <- sim_df$length[mod_num]
    sigma <- sim_df$sigma[mod_num]
    beta <- sim_df$beta[mod_num]

    # add new lines to error and out file identifying which model this is
    write(paste('model number = ', mod_num, ', nsteps = ', nsteps,
                ', sigma = ', sigma, ', beta = ', beta),
          stderr())

    # simulate AR-p data (will update in future to have variable inputs)
    input_pars <- list(
      n = nsteps,        # length of time series
      sd = sigma,        # standard deviation of the innovations
      phi = c(0.5, 0.1), # default vector of AR terms
      beta_n = 5,        # number of covariates to include
      beta_p = 0,        # number of lags of covariate 1 to include
      beta_sig = NULL,   # the number of significant betas. NULL defaults to all significant.
      beta_select = beta,# the number of beta terms to keep in the model
      K = 100,           # number of fourier components to include
      holdout = 100
    )

    sim_dat <- ARp_beta_sims(input_pars)

    # Save the results
    fname <- paste0("/simdat_run", mod_num, ".rds")
    fpath <- paste0("Data/aquatic_sim_data/", outdir, fname)
    saveRDS(sim_dat, file = paste0("/project/modelscape/analyses/sponges/", fpath))
    # saveRDS(sim_dat, file =  fpath)

}
