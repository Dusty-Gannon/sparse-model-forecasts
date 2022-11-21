# document setup
library(here)
library(parallel)

simulate_comp_communities <- function(X, sim_params){

  # store number of species
  nsp <- sim_params$S + sim_params$s

  # vector of intra-specific competitive effects
  alpha <- runif(nsp, min=3.5, max = 4)

  # competition matrix
  A_mat <- comp_matrix2(
    n_sp = nsp,
    alpha = alpha,
    rho = 0.01,
    num_ngs = sim_params$S
  )

  # max per-capita fecundity
  lambda_max <- c(
    runif(sim_params$n_annuals, min = 45, max = 50),
    runif(nsp - sim_params$n_annuals, min = 10, max = 15)
  )

  # environmental optima
  sp_optims <- runif(nsp, min = -0.2, max = 0.2)

  # probability an individual dies in a given time point for each species
  Pr_death <- c(
    rep(1, sim_params$n_annuals),
    runif(nsp - sim_params$n_annuals, min = 0.05, max = 0.5)
  )

  # dispersal rates for each species
  disp_rate <- c(
    runif(sim_params$n_annuals, 0.12, 0.2),
    runif(nsp - sim_params$n_annuals, 0.9, 1.1)
  )

  # initialize lattice
  X <- matrix(
    data = sample(1:nsp, sim_params$cells * sim_params$cells, replace = T),
    nrow = sim_params$cells, ncol = sim_params$cells
  )

  # count the neighbors for each cell given the radius and
  #  lattice size
  n_neighbors <- count_neighbors(M = sim_params$cells, J = sim_params$cells, r = sim_params$nbrhood_radius)

  # track evolution of community
  ts_all <- vector(mode = "list", length = sim_params$steps)
  ts_all[[1]] <- X

  # track environment
  sigma_env <- 0.1    # variation in the environment
  phi <- 0.4          # autocorrelation across time steps
  env <- vector(mode = "double", length = sim_params$steps)
  env[1] <- rnorm(1, sd = sigma_env/sqrt(1 - phi^2))


  # step through the evolution of the community
  for(t in 1:(sim_params$steps-1)){

    # count neighbors of each species around each cell
    nbrs_t <- kernel_count(ts_all[[t]], r = sim_params$nbrhood_radius, sp_list = 1:nsp)

    # get fecundity matrix by applying the fecundidty_ll() command to each row of X
    Fecundity <- t(sapply(
      1:sim_params$cells,
      FUN = function(x, lattice, lambda_t, alpha, nbrhood, n){
        fecundity_ll(
          foc_sp = lattice[x, ],
          lambda_max = lambda_max,
          optims = sp_optims,
          env_t = env[t],
          alpha = alpha,
          nbrhood = nbrhood[x, ],
          n = n[x, , ]
        )
      },
      lattice = ts_all[[t]],
      lambda_t = lambda,
      alpha = A_mat,
      nbrhood = n_neighbors,
      n = nbrs_t
    ))

    # calculate seed rain into each cell by each species
    sr_t <- seed_rain_array(
      F_mat = Fecundity,
      X = ts_all[[t]],
      d_max = 3,
      rate = disp_rate,
      nsp = nsp
    )

    # kill off and replace some adults with seedlings
    ts_all[[t+1]] <- die_replace(
      ts_all[[t]],
      prob_death = Pr_death,
      seed_rain = sr_t
    )

    # change environment for next iteration
    env[t + 1] <- phi * env[t] + rnorm(1, sd = sigma_env)

  }

  # get percent cover dataframe
  cover_df <- data.frame(
    t = rep(1:sim_params$steps, nsp),
    species = rep(1:nsp, each = sim_params$steps)
  )

  # fill in the cover values
  cover <- vector(mode = "double")
  for(i in 1:nsp){
    cover <- c(
      cover,
      sapply(ts_all, function(x){sum(x == i)/(nrow(x) * ncol(x))})
    )
  }
  cover_df$cover <- cover
  cover_df$species <- as.factor(cover_df$species)

  # which species died out early on?
  extinct <- which(cover_df[cover_df$t == round(sim_params$steps/2), ]$cover == 0)
  cover_df <- cover_df[-which(cover_df$species %in% extinct), ]

  # return params for select species
  params <- list(
    A_mat = A_mat[-extinct, -extinct],
    lambda = lambda_max[-extinct],
    pr_death = Pr_death[-extinct],
    disp_rate = disp_rate[-extinct],
    env_optims = sp_optims[-extinct]
  )
  message(
    paste(
      "Number of species coexisting:",
      sum(cover_df[cover_df$t == sim_params$steps, ]$cover > 0)
    )
  )

  # return objects
  return(
    list(
      params = params,
      cover = cover_df,
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
    cells = 50,
    nbrhood_radius = 3,
    steps = as.numeric(args[2]),
    S = 2, s = 48,
    n_annuals = 5
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
  fname <- paste0("simdat_", args[1], "reps_", args[2], "steps_S2s48_5ann.rds")
  fpath <- paste0("Data/terrestrial_sim_data/", fname)
  saveRDS(sim_dat, file = here(fpath))






