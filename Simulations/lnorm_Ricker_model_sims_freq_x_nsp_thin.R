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
      ricker_ts_lnorm(
        N_0 = N_0,
        lambdas = lambdas,
        A_mat = A_mat,
        sigmas = sigmas,
        steps = tot_steps,
        thin_freq = 0
      )
    })

    # which species were extinct for more than 70% of the time?
    ext <- (1:nrow(N_pre))[which(apply(N_pre, 1, function(x){
      mean(x == 0)
    }) > 0.7)]

    N <- N_pre[-ext, ]

    return(list(
      N = N,
      sim_params = list(
        lambdas = sim_params$lambdas[-ext],
        A_mat = sim_params$A_mat[-ext, -ext],
        sigmas = sim_params$sigmas[-ext],
        start_thin = NA,
        thin_freq = 0,
        prop_cthin = 0
      )
    ))

  }

  # conditional for whether the community should be thinned
  if(sim_params$thin_freq > 0){
    N_pre <- with(sim_params, {
      ricker_ts_lnorm(
        N_0 = N_0,
        lambdas = lambdas,
        A_mat = A_mat,
        sigmas = sigmas,
        steps = init_steps,
        thin_freq = 0
      )
    })

    # name species for easier tracking
    rownames(N_pre) <- 1:sim_params$nsp

    # get species that are not extinct after warmup
    nesp <- which(!(N_pre[, sim_params$init_steps] == 0))
    ext <- which(N_pre[, sim_params$init_steps] == 0)
    N_pre_final <- N_pre[nesp, ]
    foc <- nesp[1]

    # get number of species to thin
    nthin <- ceiling(sim_params$prop_cthin * (nrow(N_pre_final) - 1))

    # get thinning order based on rank abundances
    temp_means <- apply(N_pre_final[-1, round(sim_params$init_steps/2):sim_params$init_steps], 1, mean)
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
        N_pre = N_pre_final,
        nesp = nesp,
        ext = ext,
        thin_order = thin_order,
        thin_levels = thin_levels
      ),
      sim_params
    )

    # now conduct thinning experiments
    N_exp <- with(param_list_r2, {
      ricker_ts_lnorm(
        N_0 = N_pre[, init_steps],
        lambdas = lambdas[nesp],
        A_mat = A_mat[nesp, nesp],
        sigmas = sigmas[nesp],
        steps = tot_steps - init_steps + 1,
        thin_freq = thin_freq,
        thin_levels = thin_levels,
        thin_order = thin_order
      )
    })

    rownames(N_exp) <- rownames(N_pre_final)

    # find and remove any species that went extinct in the first 50 steps after thinning
    ext2 <- (2:nrow(N_exp))[which(apply(N_exp[-1, ], 1, function(x){
      mean(x == 0)
    }) > 0.7)]

    # make some return objects
    N = cbind(N_pre_final[-ext2, -sim_params$init_steps], N_exp[-ext2, ])
    thinned_sp = which(rownames(N) %in% thin_sp_names)

    return(list(
      N = N,
      thinned_sp = thinned_sp,
      sim_params = list(
        lambdas = sim_params$lambdas[nesp][-ext2],
        A_mat = sim_params$A_mat[nesp, nesp][-ext2, -ext2],
        sigmas = sim_params$sigmas[nesp][-ext2],
        start_thin = sim_params$init_steps + 1,
        thin_freq = sim_params$thin_freq,
        prop_cthin = sim_params$prop_cthin
      )
    ))
  }
}



##### Defining ranges for parameter generation #####

  # some general simulation conditions
  set.seed(5254)
  thin_freq <- 1:10
  prop_cthin <- seq(0.1, 1, by = 0.1)
  reps <- 1000
  target_reps <- 500

  # complete the factorial table
  freq_prop_df <- expand.grid(thin_freq, prop_cthin)
  names(freq_prop_df) <- c("thin_freq", "prop_cthin")
  # add a row for no thinning at all
  freq_prop_df <- rbind(c(0,0), freq_prop_df)
  freq_prop_df <- freq_prop_df[rep(1:nrow(freq_prop_df), each = reps), ]

  nsp <- 60; init_steps <- 50; tot_steps = 500
  num_ngs <- 5; sigma_rng <- c(0.1, 0.3);
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
      target_thin = target_thin
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
  cl <- makeCluster(20)

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

# ##### Quick checks #####
 # sapply(sims, function(x){sum(x$N_full[, steps] > 0)})
 #
 #  plot_comm <- function(N){
 #    library(ggplot2)
 #    nsp <- nrow(N)
 #    steps <- ncol(N)
 #
 #    df <- data.frame(
 #      t = rep(1:steps, each = nsp),
 #      sp = as.factor(rep(1:nsp, steps)),
 #      N = as.vector(N)
 #    )
 #
 #    return(
 #      ggplot(data = df, aes(x = t, y = N, color = sp)) +
 #        geom_line() +
 #        theme_classic()
 #    )
 #  }
 #
 # plot_comm(N)

##### Save the simulated datasets #####

  # find and remove any cases where the focal went extinct
  to_keep <- which(
    {
      as.vector(sapply(sims, function(x){sum(is.na(x$N))})) +
      as.vector(sapply(sims, function(x){
          if(nrow(x$N) == 0){return(1)} else{
            sum(is.nan(x$N[1, ]) | x$N[1, ] == 0)
          }
      }))
    } == 0
  )

  # sample from good sims
  grps <- length(sims) / reps
  sub_ids <- vector(mode = "double")
  for(i in 1:grps){

    good_ids <- to_keep[to_keep %in% ((i - 1) * reps + 1):(i * reps)]
    sub_ids <- c(
      sub_ids,
      sample(good_ids, target_reps)
    )
  }

  sub_ids <- sub_ids[order(sub_ids)]

  sims_final <- sims[sub_ids]


  # sanity check
  test_df <- data.frame(
    thin_freq = purrr::map_dbl(
      sims_final,
      ~ .x$sim_params$thin_freq
    ),
    prop_cthin = purrr::map_dbl(
      sims_final,
      ~ .x$sim_params$prop_cthin
    )
  )

##### Save simulations #####

  fname <- paste0(
    "lnorm_ricker_thin_freq_x_nsp_ordered",
    "_S", num_ngs,
    "_s", nsp - num_ngs, ".rds"
  )

  fp <- paste0(
    "Data/terrestrial_sim_data/lnorm_ricker/",
    fname
  )

  saveRDS(sims_final, file = here(fp))



