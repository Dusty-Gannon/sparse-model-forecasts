#################################################################
# Fitting Ricker models with log-normal demographic stochasticity
#################################################################


##### Doc setup #####

  # libraries
  library(rstan)
  library(here)
  library(parallel)
  devtools::load_all()



###########################################
# Functions to fit and summarize the models
###########################################


##### Function to fit growth models #####

  fit_growth_models <- function(N, stan_mod, tsteps){

    # compile data to feed into Stan
    N_het <- t(N[-1, tsteps])
    N_het_std <- scale(N_het)
    N_foc <- as.double(N[1, tsteps])
    n <- length(tsteps)
    tau_0 <- tau0(
      y = log(N_foc[2:n] / N_foc[1:(n - 1)]),
      m0 = min(5, ncol(N_het) - 1),
      M = ncol(N_het),
      N = n - 1,
      fam = "gaussian"
    )

    datlist <- list(
      N = n,
      P = ncol(N_het),
      y = N_foc,
      X_beta = N_het_std,
      error_scl = 0.5,
      tau0 = tau_0,
      slab_scl = 0.25,
      slab_df = 6
    )

    mfit <- sampling(
      stan_mod,
      data = datlist,
      cores = 3, chains = 3,
      control = list(adapt_delta = 0.99, max_treedepth = 15)
    )

    beta_post <- rstan::extract(mfit, pars = "beta")$beta
    beta_post %*% solve(diag(apply(N_het, 2, sd)))

    return(beta_post)

  }


##### Function to summarize results from posterior betas #####
  conf_mat_summaries <- function(beta_post, A_mat, pip = 0.9){

    nz_estims <- which(apply(beta_post, 2, function(x){
      mean(x < 0) >= pip | mean(x > 0) >= pip
    }))
    z_estims <- which(!(1:ncol(beta_post) %in% nz_estims))
    nz_true <- which(A_mat[1, -1] != 0)
    z_true <-  which(!(1:ncol(beta_post) %in% nz_true))

    conf_mat <- matrix(
      c(
        sum(nz_estims %in% nz_true),
        sum(nz_estims %in% z_true),
        sum(z_estims %in% nz_true),
        sum(z_estims %in% z_true)
      ),
      nrow = 2, ncol = 2
    )

    rmse <- RMSE_bayes(
      as.double(A_mat[1, -1]),
      ppreds = beta_post
    )

    return(
      list(
        conf_mat = conf_mat,
        rmse = rmse
      )
    )

  }

##### Wrapper function to fit and summarize #####
  fit_3_summaries <- function(X, stan_mod, tsteps = 51:100, pip = 0.9){

    N_list <- list(N_full = X$N_full, N_cor = X$N_cor, N_eqv = X$N_eqv)

    beta_posts <- lapply(
      N_list,
      FUN = fit_growth_models,
      stan_mod = stan_mod,
      tsteps = tsteps
    )

    summaries <- lapply(
      beta_posts,
      conf_mat_summaries,
      A_mat = X$sim_params$A_mat,
      pip = pip
    )

    return(summaries)

  }




##################################
# Fit models and summarize results
##################################

  # command line arguments
  args <- commandArgs(trailingOnly = T)
  datfile <- args[1]

  # load simulated data
  datfp <- paste0("Data/terrestrial_sim_data/lnorm_ricker/", datfile)
  sims <- readRDS(here(datfp))

  # compile stan model
  growth_mod <- rstan::stan_model(
    here("Stan/pop_growth_rate_FHS.stan")
  )

  # apply the above wrapper in parallel
  results <- lapply(
    sims,
    fun = fit_3_summaries,
    stan_mod = growth_mod,
    tsteps = 41:100,
    pip = 0.9
  )

  # save the results
  saveRDS(results, file = here(
    "Data/terrestrial_sim_data/lnorm_ricker/vr_tests.rds"
  ))





















