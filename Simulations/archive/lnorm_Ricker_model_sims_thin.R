# document setup
library(here)
library(tidyverse)
devtools::load_all()


###################################################################
# This script simulates competitive communities using
# a Ricker population model for each competing species.
# Demographic stochasticity enters in the exponent,
# resulting in log-normally distributed errors. We assess, after
# 100 time steps, which species are coexisting, then select a focal
# species and conduct targeted thinning of non-focal species
###################################################################


###### Wrapper function to simulate competitive communities #####
simulate_communities <- function(sim_params){

  # generate initial communities prior to thinning
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

  # get species that are not extinct after warmup
  nesp <- which(!(N_pre[, sim_params$init_steps] == 0))
  ext <- which(N_pre[, sim_params$init_steps] == 0)
  N_pre_final <- N_pre[nesp, ]
  foc <- nesp[1]

  # get thinning order based on rank abundances
  if(sim_params$target_thin){
    thin_order <- (2:nrow(N_pre_final))[
      order(apply(N_pre_final[-1, round(sim_params$init_steps/2):sim_params$init_steps], 1, mean), decreasing = T)
    ]
  } else{
    thin_order <- sample(2:nrow(N_pre_final))
  }

  param_list_r2 <- c(
    list(
      N_pre = N_pre_final,
      nesp = nesp,
      ext = ext,
      thin_order = thin_order
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
      steps = tot_steps - init_steps,
      thin_freq = thin_freq,
      thin_factor = thin_factor,
      thin_order = thin_order
    )
  })

  # find and remove any species that went extinct in the first 50 steps after thinning
  ext2 <- (2:nrow(N_exp))[which(apply(N_exp[-1, ], 1, function(x){
    mean(x == 0)
  }) > 0.7)]


  return(list(
    N = cbind(N_pre_final[-ext2, ], N_exp[-ext2, ]),
    sim_params = list(
      lambdas = sim_params$lambdas[nesp][-ext2],
      A_mat = sim_params$A_mat[nesp, nesp][-ext2, -ext2],
      sigmas = sim_params$sigmas[nesp][-ext2],
      start_thin = sim_params$init_steps + 1,
      thin_freq = sim_params$thin_freq
    )
  ))

}


##### Defining ranges for parameter generation #####

  # some general simulation conditions
  set.seed(5254)
  reps_per_thin_freq <- 500
  thin_freq <- rep(c(1, 2, 5, 10), reps_per_thin_freq)
  target_reps <- 100

  nsp <- 60; init_steps <- 50; tot_steps = 500
  num_ngs <- 5; sigma_rng <- c(0.1, 0.3);
  alpha_rng <- c(0.01, 0.05); lambda_rng <- c(1.2, 1.6);
  ng_range <- c(0.1, 0.5); rho <- 0; mean_init_abund <- 50
  thin_factor = 0.1; target_thin = T

##### Running the simulations #####

  # generate list of parameters
  params <- purrr::map(
    thin_freq,
    ~ generate_sim_params_thin(
      thin_freq = .x,
      nsp = nsp, init_steps = init_steps,
      tot_steps = tot_steps, num_ngs = num_ngs,
      sigma_rng = sigma_rng, alpha_rng = alpha_rng,
      lambda_rng = lambda_rng, ng_range = ng_range,
      rho = rho, mean_init_abund = mean_init_abund,
      thin_factor = thin_factor,
      target_thin = target_thin
    )
  )

  sims <- lapply(
    params,
    FUN = simulate_communities
  )

# ##### Quick checks #####
#  sapply(sims, function(x){sum(x$N_full[, steps] > 0)})
#
#   plot_comm <- function(N){
#     library(ggplot2)
#     nsp <- nrow(N)
#     steps <- ncol(N)
#
#     df <- data.frame(
#       t = rep(1:steps, each = nsp),
#       sp = as.factor(rep(1:nsp, steps)),
#       N = as.vector(N)
#     )
#
#     return(
#       ggplot(data = df, aes(x = t, y = N, color = sp)) +
#         geom_line() +
#         theme_classic()
#     )
#   }
#
#  plot_comm(N)

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
  grps <- length(sims) / reps_per_thin_freq
  sub_ids <- vector(mode = "double")
  for(i in 1:grps){

    good_ids <- to_keep[to_keep %in% ((i - 1) * reps_per_thin_freq + 1):(i * reps_per_thin_freq)]
    sub_ids <- c(
      sub_ids,
      sample(good_ids, target_reps)
    )
  }

  sub_ids <- sub_ids[order(sub_ids)]

  sims_final <- sims[sub_ids]


##### Save simulations #####

  fname <- paste0(
    "lnorm_ricker_thin_sims_ordered",
    "_S", num_ngs,
    "_s", nsp - num_ngs, ".rds"
  )

  fp <- paste0(
    "Data/terrestrial_sim_data/lnorm_ricker/",
    fname
  )

  saveRDS(sims_final, file = here(fp))



