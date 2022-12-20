# document setup
library(here)
library(parallel)


#########################################################
# This function simulates a competitive community using
# a Ricker population model for each competing species.
# environmental variation affects growth rates as well as
# the strength of competition among some species.
#########################################################

simulate_Ricker_comp_communities <- function(X, sim_params){

  # store number of species
  nsp <- sim_params$S + sim_params$s

  # track 2 environmental variables
  # one stationary, one increasing through time
  phi <- 0.4                               # autocorrelation across time steps
  env <- matrix(data = 0, nrow = sim_params$steps, ncol = 2)
  env[1, 1] <- rnorm(1, sd = sim_params$sigma_env[1]/sqrt(1 - phi^2))

  # generate autocorrelated but increasing time series
  mu_env2 <- seq(-0.5, 0.5, length.out = sim_params$steps)
  env[1, 2] <- rnorm(1, mean = mu_env2[1], sd = sim_params$sigma_env[2])

  # generate matrix of parameters that determine competition based
  #  on number of non-generics and environmental covariates
  B <- ricker_comp_array(
    nsp = nsp, num_ngs = sim_params$S,
    env = c(0.1),
    num_dynamic = 2, generic = 0.001
  )

  # max per-capita fecundity
  lambda_max <- c(
    runif(nsp, min = 1.2, max = 2)
  )

  # environmental optima
  sp_optims <- runif(nsp, min = -0.2, max = 0.2)

  # track abundances through time
  N_ts <- matrix(data = 0, nrow = nsp, ncol = sim_params$steps)
  N_ts[, 1] <- rpois(nsp, lambda = 10)

  # step through the evolution of the community
  for(t in 1:(sim_params$steps-1)){

    # determine competition
    env_comp_t <- c(1, 1, env[t, 2])
    A_mat_t <- t(sapply(
      1:nsp,
      function(x, A, vec){
        A[, , x] %*% vec
      },
      A = B$B,
      vec = env_comp_t
    ))

    # growth rates for each step for each population
    lambda_t <- gauss_env_effect(env[t, 1], lambda_max, optims = sp_optims)

    # determine species abundance in the next step
    N_ts[, t + 1] <- ricker_step(
      N_ts[, t],
      lambdas = lambda_t,
      A_mat = A_mat_t,
      stochastic = F
    )

    # change environment for next iteration
    env[t + 1, 1] <- phi * env[t, 1] + rnorm(1, sd = sim_params$sigma_env[1])
    env[t + 1, 2] <- mu_env2[t + 1] + phi * (env[t, 2] - mu_env2[1]) +
      rnorm(1, sd = sim_params$sigma_env[2])

  }

  # abundance dataframe (long)
  df_long <- data.frame(
    t = rep(1:sim_params$steps, each = nsp),
    species = as.factor(rep(1:nsp, sim_params$steps)),
    abundance = as.vector(N_ts)
  )

  ggplot2::ggplot(data = df_long, ggplot2::aes(x= t, y = abundance, color = species))+
    ggplot2::geom_line()

  # get cover and fecundity values for a randomly selected sub-plot
  ul_corner <- sample(1:(sim_params$cells - sim_params$sub_cells), size = 1)
  rc_ids <- ul_corner:(ul_corner + sim_params$sub_cells - 1)

  # subset the lattice
  ts_sub <- lapply(
    ts_all,
    FUN = function(x, rc_ids){x[rc_ids, rc_ids]},
    rc_ids = rc_ids
  )
  ts_fec_sub <- lapply(
    ts_fecundity,
    FUN = function(x, rc_ids){x[rc_ids, rc_ids]},
    rc_ids = rc_ids
  )

  # initialize dataframes
  subplot_cover <- data.frame(
    t = rep(1:sim_params$steps, nsp),
    species = rep(1:nsp, each = sim_params$steps)
  )
  subplot_fec <- data.frame(
    t = rep(1:(sim_params$steps - 1), nsp),
    species = rep(1:nsp, each = sim_params$steps - 1)
  )

  # fill in the cover values
  cover <- vector(mode = "double")
  subcover <- vector(mode = "double")
  for(i in 1:nsp){
    cover <- c(
      cover,
      sapply(ts_all, function(x){sum(x == i)/(nrow(x) * ncol(x))})
    )
    subcover <- c(
      subcover,
      sapply(ts_sub, function(x){sum(x == i)/(nrow(x) * ncol(x))})
    )
  }
  cover_df$cover <- cover
  subplot_cover$cover <- subcover
  cover_df$species <- as.factor(cover_df$species)
  subplot_cover$species <- as.factor(subplot_cover$species)

  # fill in fecundity values
  subfec <- vector(mode = "double")
  for(i in 1:nsp){
    subfec <- c(
      subfec,
      sapply(
        1:length(ts_fecundity),
        function(x, mat, id_mat, s){
          sum((id_mat[[x]] == s) * mat[[x]])
        },
        mat = ts_fec_sub,
        id_mat = ts_sub,
        s = i
      )
    )
  }
  subplot_fec$fecundity <- subfec
  subplot_fec$species <- as.factor(subplot_fec$species)

  # print message about the number of species coexisting
  #  at the end of the simulation
  message(
    paste(
      "Number of species coexisting:",
      sum(cover_df[cover_df$t == sim_params$steps, ]$cover > 0)
    )
  )

  # return objects
  return(
    list(
      params = c(
        sim_params,
        A_mat = A_mat,
        lambda_max = lambda_max,
        sp_optims = sp_optims,
        Pr_death = Pr_death
      ),
      cover = cover_df,
      subplot_cover = subplot_cover,
      subplot_fecundity = subplot_fec,
      env = env
    )
  )
}


#### Run the simulations and store datasets ####

  args <- commandArgs(trailingOnly = T)
  # define number of sims
  nsims <- as.numeric(args[1])

  # set sim parameters
  sim_params <- list(
    steps = as.numeric(args[2]),
    S = 2, s = 48,
    n_annuals = 5,
    sigma_env = as.numeric(args[3])
  )

  # make the cluster
  cl <- makeCluster(20)

  # load necessary functions
  clusterEvalQ(
    cl = cl,
    expr = {
      library(here)
      src_files <- list.files(here("R/"), pattern = "*.R", full.names = T)
      sapply(src_files, source, .GlobalEnv)
    }
  )

  # run the simulations in parallel
  sim_dat <- parLapply(
    cl = cl,
    X = 1:nsims,
    fun = simulate_comp_communities,
    sim_params = sim_params
  )

  # stop the cluster and save the results
  stopCluster(cl)
  fname <- paste0("simdat_", args[1], "reps_", args[2], "steps_S2s48_5ann_env", args[3], ".rds")
  fpath <- paste0("Data/terrestrial_sim_data/", fname)
  saveRDS(sim_dat, file = here(fpath))






