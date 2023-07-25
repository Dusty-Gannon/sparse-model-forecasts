# document setup
library(here)
library(tidyverse)
library(parallel)
devtools::load_all()


###################################################################
# This script simulates competitive communities using
# a Ricker population model for each competing species.
# Demographic stochasticity enters in the exponent,
# resulting in log-normally distributed errors. We assess, after
# 50 time steps, which species are coexisting, then select a focal
# species and conduct targeted thinning of non-focal species
###################################################################


###### Wrapper function to simulate competitive communities #####
simulate_communities <- function(sim_params){

  # generate initial communities prior to thinning
  if(sim_params$thin_freq == 0){
    N_pre <- with(sim_params, {
      ricker_spts_lnorm(
        nsp = nsp,
        N_0 = N_0,
        lambdas = lambdas,
        A_mat = A_mat,
        sigmas = sigmas,
        steps = tot_steps,
        thin_freq = 0,
        n_sites = n_sites,
        sp_cov = sp_cov
      )
    })$N

    # remove species that were extinct most of the time
    time_extinct <- apply(
      N_pre, 3,
      function(X){
        apply(
          X, 1,
          function(x){
            mean(x == 0)
          }
        )
      }
    )
    sp2rm <- which(
      apply(
        time_extinct, 1,
        function(x){
          sum(x > 0.5)
        }
      ) > 0.5 * sim_params$n_sites
    )
    N <- N_pre[-sp2rm, , ]

    # create new order to label focal
    foc <- which.min(
      apply(N, 1, FUN = function(X){
        mean(X == 0)
      })
    )
    if(foc == 1){
      new_ord <- 1:dim(N)[1]
    } else if(foc == dim(N)[1]){
      new_ord <- c(
        foc,
        1:(dim(N)[1] - 1)
      )
    } else{
      new_ord <- c(
        foc,
        1:(foc - 1),
        (foc + 1):dim(N)[1]
      )
    }
    N <- N[new_ord, , ]

    return(list(
      N = N,
      sim_params = list(
        lambdas = sim_params$lambdas[-sp2rm][new_ord],
        A_mat = sim_params$A_mat[-sp2rm, -sp2rm][new_ord, new_ord],
        sigmas = sim_params$sigmas[-sp2rm][new_ord],
        start_thin = NA,
        thin_freq = 0,
        prop_cthin = 0
      )
    ))

  }

  # conditional for whether the community should be thinned
  if(sim_params$thin_freq > 0){
    N_pre <- with(sim_params, {
      ricker_spts_lnorm(
        nsp = nsp,
        N_0 = N_0,
        lambdas = lambdas,
        A_mat = A_mat,
        sigmas = sigmas,
        steps = init_steps,
        thin_freq = 0,
        n_sites = n_sites,
        sp_cov = sp_cov
      )
    })$N

    # name species for easier tracking
    rownames(N_pre) <- 1:sim_params$nsp

    # get species that are not extinct after warmup
    time_extinct <- apply(
      N_pre, 3,
      function(X){
        apply(
          X, 1,
          function(x){
            mean(x == 0)
          }
        )
      }
    )
    sp2rm <- which(
      apply(
        time_extinct, 1,
        function(x){
          sum(x > 0.5)
        }
      ) > 0.5 * sim_params$n_sites
    )
    N_pre_final <- N_pre[-sp2rm, , ]
    nesp <- (1:sim_params$nsp)[!(1:sim_params$nsp %in% sp2rm)]
    foc <- nesp[1]

    # get number of species to thin
    nthin <- ceiling(sim_params$prop_cthin * (nrow(N_pre_final) - 1))

    # get thinning order based on rank abundances
    temp_means <-
      apply(
        N_pre_final[-1, round(sim_params$init_steps/2):sim_params$init_steps, ],
        1, mean
      )
    if(sim_params$target_thin){
      thin_order <- (2:nrow(N_pre_final))[
        order(temp_means, decreasing = T)
      ][1:nthin]
      thin_sp_names <- rownames(N_pre_final)[thin_order]
      thin_levels <- (sort(temp_means, decreasing = T) * sim_params$thin_factor)[1:nthin]
    } else{
      thin_order <- sample(2:nrow(N_pre_final), size = nthin)
      thin_levels <- temp_means[thin_order] * sim_params$thin_factor
    }

    # combine original list with newly created variables
    param_list_r2 <- c(
      list(
        N_pre_f = N_pre_final,
        nesp = nesp,
        ext = sp2rm,
        thin_order = thin_order,
        thin_levels = thin_levels
      ),
      sim_params
    )

    # now conduct thinning experiments
    N_exp <- with(param_list_r2, {
      ricker_spts_lnorm(
        nsp = dim(N_pre_f)[1],
        N_0 = N_pre_f[, init_steps, ],
        lambdas = lambdas[nesp],
        A_mat = A_mat[nesp, nesp],
        sigmas = sigmas[nesp],
        steps = tot_steps - init_steps + 1,
        n_sites = n_sites,
        sp_cov = sp_cov,
        thin_freq = thin_freq,
        thin_levels = thin_levels,
        thin_order = thin_order
      )
    })$N

    rownames(N_exp) <- rownames(N_pre_final)

    # find and remove any species that went extinct after thinning
    time_extinct2 <- apply(
      N_exp, 3,
      function(X){
        apply(
          X, 1,
          function(x){
            mean(x == 0)
          }
        )
      }
    )
    sp2keep <- which(
      apply(
        time_extinct2, 1,
        function(x){
          sum(x < 0.5)
        }
      ) > 0.5 * sim_params$n_sites
    )

    # ensure that the unthinned focal species is included
    if(!(1 %in% sp2keep)){
      sp2keep <- c(1, sp2keep)
    }

    # make some return objects
    N <- abind::abind(
      N_pre_final[sp2keep, , ], N_exp[sp2keep, -1, ],
      along = 2
    )
    thinned_sp <- which(dimnames(N)[[1]] %in% thin_sp_names)

    return(list(
      N = N,
      thinned_sp = thinned_sp,
      sim_params = list(
        lambdas = sim_params$lambdas[nesp][sp2keep],
        A_mat = sim_params$A_mat[nesp, nesp][sp2keep, sp2keep],
        sigmas = sim_params$sigmas[nesp][sp2keep],
        start_thin = sim_params$init_steps + 1,
        thin_freq = sim_params$thin_freq,
        prop_cthin = sim_params$prop_cthin
      )
    ))
  }
}



##### Defining ranges for parameter generation #####

  # load arguments from command line
  args <- commandArgs(trailingOnly = TRUE)

  # some general simulation conditions
  thin_freq <- 1:10
  prop_cthin <- seq(0.1, 1, by = 0.1)
  reps <- 10
  target_reps <- 1
#  round <- 2

  # complete the factorial table
  freq_prop_df <- expand.grid(thin_freq, prop_cthin)
  names(freq_prop_df) <- c("thin_freq", "prop_cthin")
  # add a row for no thinning at all
  freq_prop_df <- rbind(c(0,0), freq_prop_df)
  freq_prop_df <- freq_prop_df[rep(1:nrow(freq_prop_df), reps), ]

  nsp <- 60; init_steps <- 20; tot_steps = 100
  n_sites <- 20; num_ngs <- 5; sigma_rng <- c(0.1, 0.3);
  alpha_rng <- c(0.01, 0.05); lambda_rng <- c(1.2, 1.6);
  ng_range <- c(0.1, 0.5); rho <- 0; mean_init_abund <- 50
  thin_factor = 0.1; target_thin = T

##### Running the simulations #####

  # generate list of parameters
  params <- purrr::map(
    1:reps,
    ~ generate_sim_params_thin(
      thin_freq = 1,
      prop_cthin = 0,
      nsp = nsp, init_steps = init_steps,
      tot_steps = tot_steps, num_ngs = num_ngs,
      sigma_rng = sigma_rng, alpha_rng = alpha_rng,
      lambda_rng = lambda_rng, ng_range = ng_range,
      rho = rho, mean_init_abund = mean_init_abund,
      thin_factor = thin_factor,
      target_thin = target_thin,
      spatial = TRUE,
      n_sites = n_sites
    )
  )

  # now expand the parameter list to repeat each set of parameters for each combination of
  # thinning frequency and proportion of cummunity that gets thinned
  params <- params[rep(1:length(params), each = nrow(unique(freq_prop_df)))]
  for(i in 1:length(params)){
    params[[i]]$thin_freq <- freq_prop_df$thin_freq[i]
    params[[i]]$prop_cthin <- freq_prop_df$prop_cthin[i]
  }

  # simulate in parallel
  cl <- makeCluster(reps)

  # load functions on each node
  clusterEvalQ(cl, expr = {
    devtools::load_all()
  })

  # simulate
  sims <- parLapply(
    cl = cl,
    params,
    fun = simulate_communities
  )

  stopCluster(cl)

  # collect a list of replicated sim params for which the focal
  # did not go extinct in any treatment combo
  r <- 1
  good_reps <- 0
  n_trtmnts <- nrow(unique(freq_prop_df))
  sims_final <- vector(mode = "list")
  while(good_reps < target_reps & r <= reps){

    # get sim ids with same starting parameters
    ids <- (n_trtmnts * (r - 1) + 1):(n_trtmnts * r)

    # loop through and check for issues, stop if one is
    # encountered and move on
    prblms <- 0
    i <- 1
    while(prblms == 0 & i <= length(ids)){
      if(dim(sims[[ids[i]]]$N)[1] > 0){
        prblms <- prblms +
          sum(is.infinite(sims[[ids[i]]]$N[1, ,])) +
          (mean(sims[[ids[i]]]$N[1,,] == 0) > 0.05)
      } else {
        prblms <- prblms + 1
      }
      i <- i + 1
    }

    if(prblms == 0){
      # add the sims if there were no issues
      sims_final <- c(sims_final, sims[ids])
      good_reps <- good_reps + 1
    }
    r <- r + 1

  }

  # remove sims
  rm(sims)


##### Save simulations #####

  if(length(sims_final) > 0){
    fname <- paste0(
      "lnorm_ricker_thin_freq_x_nsp_s4t",
      "_rep_", args[1], ".rds"
    )

    fp <- paste0(
      args[2],
      fname
    )

    saveRDS(sims_final, file = here(fp))
  }



