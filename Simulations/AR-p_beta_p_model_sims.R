
# libraries
library(rstan)
#library(here)
devtools::load_all()
#rstan_options(auto_write = TRUE)
setwd("/project/modelscape/analyses/sponges/")

#### Simulate data and run regularized and non-regularized AR-p_beta models ####
# store data inputs and model fits in a log file


sigmas <- c(0.5, 2, 5)
# lengths <- c(60, 75, 120, 180, 240)
lengths <- c(rep(60, 12), rep(65, 12), rep(70, 11), rep(75, 11),
             rep(80, 10), rep(85, 10), rep(90, 8), rep(95, 8),
             rep(100, 2), rep(105, 2), rep(110, 2),
             seq(115, 200, by = 5))
sim_df <- data.frame(
    sigma = rep(sigmas, each = length(lengths)),
    length = rep(lengths, length(sigmas))
)

full_sim_df <- data.frame()
for(i in 1:60){
  full_sim_df <- rbind(full_sim_df, sim_df)
}
mods_per_node <- 20

array_size <- nrow(full_sim_df)/mods_per_node
nsims <- array_size

# fits <- readRDS('Data/aquatic_sim_data/test_fit.rds')

ARp_beta_sims <- function(input_pars){

    model_pars <- simulate_AR_p_beta_p_timeseries(input_pars, draw_beta = 'near_zero',
                                                  draw_phi = 'zero')
    fits <- fit_ARp_beta_model(model_pars, iter = 6000, warmup = 4000)
    sim_list <- unpack_ARp_fit(fits = fits)

    return(sim_list)

}


args <- commandArgs(trailingOnly = TRUE)

#number of simulations to run:
i <- as.numeric(args[2]) # which array instance is this
#K = nsims * length(sigmas) * length(lengths) / array_size # how many models/batch
K = mods_per_node

#j <- K*(i-1) + 1 # where in the sim dataframe to start based on which batch this is

for(k in 1:K){

    mod_num = (i-1)*K + k

    nsteps <- full_sim_df$length[mod_num]
    sigma <- full_sim_df$sigma[mod_num]

    # add new lines to error and out file identifying which model this is
    write(paste('model number = ', mod_num, ', nsteps = ', nsteps, ', sigma = ', sigma),
          stderr())

    # simulate AR-p data (will update in future to have variable inputs)
    input_pars <- list(
        n = nsteps,   # length of time series
        p  = 10,      # number of AR lags to consider
        beta_p = 5,   # number of beta lags to consider in lagged covariate (beta_1)
        beta_n = 45,  # number of additional covariates to include
        b0 = 0,       # intercept
        n_phi = 3,    # number of non-zero autoregressive terms
        n_lags = 2,   # number of non-zero lags in covariate
        n_beta = 2,   # number of non-zero covariate parameters
        non_zero_coef_guess = 5, # guess for the number of non-zero coefficients
        sigma_e = sigma,  #standard deviation of the innovations
        holdout = 50
    )

    sim_dat <- ARp_beta_sims(input_pars)

    # Save the results
    fname <- paste0("/simdat_run", mod_num, ".rds")
    fpath <- paste0("Data/aquatic_sim_data/", args[1], fname)
    saveRDS(sim_dat, file = paste0("/project/modelscape/analyses/sponges/", fpath))

}
