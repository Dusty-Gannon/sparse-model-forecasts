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
    num_env = 1, generic = 0.0001,
    num_dynamic = sim_params$num_dynamic,
    ng_strength = c(15, 20),
    intra_strength = c(90, 100)
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

    # set sensititivy of growth rates to changes in environment
    tau <- 1
    # growth rates for each step for each population
    lambda_t <- gauss_env_effect(env[t, 1], lambda_max, optims = sp_optims, tau = tau)

    # determine species abundance in the next step
    N_ts[, t + 1] <- ricker_step(
      N_ts[, t],
      lambdas = lambda_t,
      A_mat = A_mat_t,
      stochastic = T
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
    abundance = as.vector(N_ts),
    env1 = rep(env[, 1], each = nsp),
    env2 = rep(env[, 2], each = nsp)
  )

  # ggplot(data = df_long, aes(x = t, y = abundance, color = species))+
  #   geom_line()

  df_wide <- as.data.frame(
    cbind(
      1:sim_params$steps,
      env,
      t(N_ts)
    )
  )

  # make sure columns have distinguishable names
  names(df_wide) <- c(
    "t",
    paste0("v", 1:ncol(env)),
    paste0("s", 1:nsp)
  )

  return(list(
    B_mat = B$B,
    dynamic_sp = B$dynamic_spids,
    lambda_max = lambda_max,
    tau = tau,
    sp_optims = sp_optims,
    abundance_long = df_long,
    abundance_wide = df_wide,
    sim_params = sim_params
  ))

}


#### Run the simulations and store datasets ####

  # bring in command line arguments
  args <- commandArgs(trailingOnly = T)

  # define number of sims
  nsims <- as.numeric(args[1])

  # set sim parameters
  sim_params <- list(
    steps = as.numeric(args[2]),
    S = 2, s = 28,
    num_dynamic = 10,
    sigma_env = c(0.1, 0.05)
  )

  # make the cluster
  cl <- makeCluster(4)

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
    fun = simulate_Ricker_comp_communities,
    sim_params = sim_params
  )

  # stop the cluster and save the results
  stopCluster(cl)
  fname <- paste0(
    "simdat_500reps_", sim_params$steps, "steps_S",
    sim_params$S, "s", sim_params$s, "_",
    sim_params$num_dynamic, "dyn", "_",
    length(sim_params$sigma_env), "env", ".rds"
  )
  fpath <- paste0(here("Data/terrestrial_sim_data/"), fname)
  saveRDS(sim_dat, file = fpath)


# save plots from a randomly selected 50 simulations
  samps <- sample(1:nsims, size = 25)
  sub_sims <- sim_dat[samps]

  plots <- lapply(
    sub_sims,
    FUN = function(X){
      ggplot2::ggplot(data = X$abundance_long, ggplot2::aes(x = t, y = abundance, color = species))+
        ggplot2::geom_line()+
        ggplot2::theme(legend.position = "none")
    }
  )

  fname_plots <- paste0(
    "plots_", nsims, "reps_",
    sim_params$steps, "steps_",
    "S", sim_params$S, "s", sim_params$s, "_",
    sim_params$num_dynamic, "dyn_",
    length(sim_params$sigma_env), "env.png"
  )

  fpath_plots <- paste0("Data/terrestrial_sim_data/", fname_plots)

  png(
    filename = here(fpath_plots),
    width = 3000, height = 3000,
    units = "px", res = 300
  )
    gridExtra::grid.arrange(
      grobs = plots,
      nrow = 5, ncol = 5
    )
  dev.off()


