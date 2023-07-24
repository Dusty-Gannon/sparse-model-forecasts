# document setup
library(here)
library(tidyverse)
library(parallel)
devtools::load_all()


###################################################################
# This script simulates competitive communities using
# a Ricker population model for each competing species.
# Demographic stochasticity enters in the exponent,
# resulting in log-normally distributed errors. Stochastic
# disturbances can be introduced with the frequency of disturbance,
# disturbance intensity, and the proportion of the community the
# disturbance impacts controlled by simulation parameters.
# see ?ricker_ts_lnorm()
###################################################################


###### Wrapper function to simulate competitive communities #####
simulate_communities <- function(sim_params){

    sims_pre <- with(sim_params, {
      ricker_spts_lnorm(
        N_0 = N_0,
        lambdas = lambdas,
        A_mat = A_mat,
        sigmas = sigmas,
        steps = steps,
        sp_cov = sp_cov,
        n_sites = n_sites,
        dist_prob = dist_prob,
        dist_int = dist_int,
        prop_cdist = prop_cdist
      )
    })

    # remove species that were extinct most of the time
    time_extinct <- apply(
      sims_pre$N, 3,
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

    N <- sims_pre$N[-sp2rm, , ]

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

    # format and reorder disturbances
    if(sim_params$dist_prob > 0){
      disturbances <- lapply(
        sims_pre$disturbances,
        FUN = function(X){
          X[-sp2rm, ][new_ord, ]
        }
      )
    } else{
      disturbances <- NULL
    }

    if(is.null(nrow(N))){
      return(list(N = NULL))
    } else{
      return(list(
        N = N,
        disturbances = disturbances,
        sim_params = list(
          lambdas = sim_params$lambdas[-sp2rm][new_ord],
          A_mat = sim_params$A_mat[-sp2rm, -sp2rm][new_ord, new_ord],
          sigmas = sim_params$sigmas[-sp2rm][new_ord],
          dist_prob = sim_params$dist_prob,
          dist_int = sim_params$dist_int,
          prop_cdist = sim_params$prop_cdist
        )
      ))
    }

}




##### Getting command line args #####

args <- commandArgs(trailingOnly = T)


##### Defining ranges for parameter generation #####

  # some general simulation conditions
  # set.seed(6528)
  dist_prob <- seq(0.05, 0.1, length.out = 10)
  prop_cdist <- seq(0.1, 1, by = 0.1)
  dist_int <- c(0.5, 0.8)
  reps <- 10
  target_reps <- 1

  # complete the factorial table
  trt_df <- expand.grid(dist_prob, prop_cdist, dist_int)
  names(trt_df) <- c("dist_prob", "prop_cdist", "dist_int")
  # add a row for no disturbance
  trt_df <- rbind(c(0, 0, 0), trt_df)
  # store number of treatments
  n_trtmnts <- nrow(trt_df)
  # now expand df
  trt_df <- trt_df[rep(1:nrow(trt_df), reps), ]

  nsp <- 60; steps <- 100
  n_sites <- 20
  num_ngs <- 5; sigma_rng <- c(0.1, 0.3);
  alpha_rng <- c(0.01, 0.05); lambda_rng <- c(1.2, 1.6);
  ng_range <- c(0.1, 0.5); rho <- 0; mean_init_abund <- 50

##### Running the simulations #####

  # generate list of parameters
  params <- purrr::map(
    1:reps,
    ~ generate_sim_params_dist(
        nsp = nsp, steps = steps, num_ngs = num_ngs,
        sigma_rng = sigma_rng, alpha_rng = alpha_rng,
        lambda_rng = lambda_rng, ng_range = ng_range,
        rho = rho, mean_init_abund = mean_init_abund,
        dist_prob = 0, dist_int = 0, prop_cdist = 0,
        dist_min_thresh = 1, spatial = T, n_sites = 20
    )
  )

# now expand the parameter list to repeat each set of parameters for each combination of
  # disturbance parameters
  params <- params[rep(1:length(params), each = n_trtmnts)]
  for(i in 1:length(params)){
    params[[i]]$dist_prob <- trt_df$dist_prob[i]
    params[[i]]$prop_cdist <- trt_df$prop_cdist[i]
    params[[i]]$dist_int <- trt_df$dist_int[i]
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


##### Get a list of successful sims #####

  # collect a list of replicated sim params for which the focal
  # did not go extinct in any treatment combo
  r <- 1
  good_reps <- 0
  sims_final <- vector(mode = "list")
  while(good_reps < target_reps & r <= reps){

    # get sim ids with same starting parameters
    ids <- (n_trtmnts * (r - 1) + 1):(n_trtmnts * r)

    # loop through and check for issues, stop if one is
    # encountered and move on
    prblms <- 0
    i <- 1
    while(prblms == 0 & i <= length(ids)){
      if(!is.null(sims[[ids[i]]]$N)){
        if(dim(sims[[ids[i]]]$N)[1] > 1){
          prblms <- prblms +
            sum(is.infinite(sims[[ids[i]]]$N[1, , ])) +
            (mean(sims[[ids[i]]]$N[1,,] == 0) > 0.05)
        }
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

  # # a little cleanup of species that were extinct over the fitting period
  # for(i in 1:length(sims_final)){
  #   N <- sims_final[[i]]$N
  #   ext3 <- which(apply(
  #     scale(t(N[, 51:250])),
  #     2,
  #     FUN = function(x){
  #       sum(is.nan(x))
  #     }
  #   ) > 0)
  #   if(length(ext3) > 0){
  #     sims_final[[i]]$N <- N[-ext3, ]
  #     sims_final[[i]]$sim_params$lambdas <- sims_final[[i]]$sim_params$lambdas[-ext3]
  #     sims_final[[i]]$sim_params$sigmas <- sims_final[[i]]$sim_params$sigmas[-ext3]
  #     sims_final[[i]]$sim_params$A_mat <- sims_final[[i]]$sim_params$A_mat[-ext3, -ext3]
  #     if(!is.null(sims_final[[i]]$dist_foc)){
  #       sims_final[[i]]$dist_foc <- sims_final[[i]]$dist_foc[-ext3, ][1, ]
  #       # convert to 0 and 1
  #       sims_final[[i]]$dist_foc[sims_final[[i]]$dist_foc != 0] <- 1
  #     }
  #   } else if(!is.null(sims_final[[i]]$dist_foc)){
  #     sims_final[[i]]$dist_foc <- sims_final[[i]]$dist_foc[1, ]
  #     # convert to 0 and 1
  #     sims_final[[i]]$dist_foc[sims_final[[i]]$dist_foc != 0] <- 1
  #   }
  #
  # }

  # remove sims
  rm(sims)

  # sims_final_sort <- sims_final[
  #   order(
  #     purrr::map_dbl(sims_final, ~.x$sim_params$dist_int),
  #     purrr::map_dbl(sims_final, ~.x$sim_params$prop_cdist),
  #     purrr::map_dbl(sims_final, ~.x$sim_params$dist_prob)
  #   )
  # ]



##### Save simulations #####

  if(length(sims_final) > 0){
    fname <- paste0(
      "lnorm_ricker_dist_freq_x_nsp_x_int",
      "_rep_", args[1], ".rds"
    )

    fp <- paste0(
      args[2],
      fname
    )

    saveRDS(sims_final, file = here(fp))
  }


