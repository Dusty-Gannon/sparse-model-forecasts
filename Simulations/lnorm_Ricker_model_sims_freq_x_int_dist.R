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

    sims_pre <- with(sim_params, {
      ricker_ts_lnorm(
        N_0 = N_0,
        lambdas = lambdas,
        A_mat = A_mat,
        sigmas = sigmas,
        steps = steps,
        dist_prob = dist_prob,
        dist_int = dist_int,
        prop_cdist = prop_cdist
      )
    })

    # which species were extinct for more than 70% of the time?
    if(sim_params$dist_prob > 0){
      ext <- which(apply(sims_pre$N, 1, function(x){
        mean(x == 0)
      }) > 0.7)
      N <- sims_pre$N[-ext, ]

      # when was the focal disturbed?
      dist_foc <- sapply(
        sims_pre$disturbances,
        FUN = function(x, ext){
          as.numeric(x[-ext][1] != 0)
        },
        ext = ext
      )
    } else{
      ext <- which(apply(sims_pre, 1, function(x){
        mean(x == 0)
      }) > 0.7)
      N <- sims_pre[-ext, ]

      dist_foc <- NA
    }

    return(list(
      N = N,
      dist_foc = dist_foc,
      sim_params = list(
        lambdas = sim_params$lambdas[-ext],
        A_mat = sim_params$A_mat[-ext, -ext],
        sigmas = sim_params$sigmas[-ext],
        dist_prob = sim_params$dist_prob,
        dist_int = sim_params$dist_int,
        prop_cdist = sim_param$prop_cdist
      )
    ))

}




##### Defining ranges for parameter generation #####

  # some general simulation conditions
  set.seed(6528)
  dist_prob <- seq(0.05, 0.5, length.out = 10)
  prop_cdist <- seq(0.1, 1, by = 0.1)
  dist_int <- c(0.5, 0.9)
  reps <- 10
  target_reps <- 100

  # complete the factorial table
  trt_df <- expand.grid(dist_prob, prop_cdist, dist_int)
  names(trt_df) <- c("dist_prob", "prop_cdist", "dist_int")
  # add a row for no disturbance
  trt_df <- rbind(c(0, 0, 0), trt_df)
  trt_df <- trt_df[rep(1:nrow(trt_df), reps), ]

  nsp <- 60; steps <- 300
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
        dist_prob = 0, dist_int = 0, prop_cdist = 0
    )
  )

# now expand the parameter list to repeat each set of parameters for each combination of
  # thinning frequency and proportion of cummunity that gets thinned
  params <- params[rep(1:length(params), each = nrow(unique(trt_df)))]
  for(i in 1:length(params)){
    params[[i]]$dist_prob <- trt_df$dist_prob[i]
    params[[i]]$prop_cdist <- trt_df$prop_cdist[i]
    params[[i]]$dist_int <- trt_df$dist_int[i]
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



