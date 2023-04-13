
# libraries
library(rstan)
#library(here)
devtools::load_all()
#rstan_options(auto_write = TRUE)
setwd("/project/modelscape/analyses/sponges/")
# setwd("C:/Users/alice.carter/git/sponges/")

#### Simulate data and run regularized and non-regularized AR-p_beta models ####
# store data inputs and model fits in a log file

array_size <- 300
nsims <- 300
sigmas <- c(0.1, 1, 10)
lengths <- c(50, 100, 150, 250, 300)

sim_df <- data.frame(
    sigma = rep(sigmas, each = nsims * length(lengths)),
    length = rep(rep(lengths, each = nsims), length(sigmas))
)

sim_df[sample(1:nrow(sim_df), nrow(sim_df), replace = FALSE), ]

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
K = nsims * length(sigmas) * length(lengths) / array_size # how many models/batch

j <- K*(i-1) + 1 # where in the sim dataframe to start based on which batch this is

for(k in 0:K){

    nsteps <- sim_df$length[j+k]
    sigma <- sim_df$sigma[j+k]

    # add new lines to error and out file identifying which model this is
    write(paste('model number = ', k+j, ', nsteps = ', nsteps, ', sigma = ', sigma),
          stderr())

    # simulate AR-p data (will update in future to have variable inputs)
    input_pars <- list(
        n = nsteps,   # length of time series
        p  = 16,      # number of AR lags to consider
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
    fname <- paste0("/simdat_run", j+k, ".rds")
    fpath <- paste0("Data/aquatic_sim_data/", args[1], fname)
    saveRDS(sim_dat, file = paste0("/project/modelscape/analyses/sponges/", fpath))

}
