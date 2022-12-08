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

  # track fecundity through time
  ts_fecundity <- vector(mode = "list", length = sim_params$steps - 1)

  # track environment
  sigma_env <- sim_params$sigma_env    # variation in the environment
  phi <- 0.4                           # autocorrelation across time steps
  env <- vector(mode = "double", length = sim_params$steps)
  env[1] <- rnorm(1, sd = sigma_env/sqrt(1 - phi^2))


  # step through the evolution of the community
  for(t in 1:(sim_params$steps-1)){

    # count neighbors of each species around each cell
    nbrs_t <- kernel_count(ts_all[[t]], r = sim_params$nbrhood_radius, sp_list = 1:nsp)

    # get fecundity matrix by applying the fecundity_ll() command to each row of X
    ts_fecundity[[t]] <- t(sapply(
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
      F_mat = ts_fecundity[[t]],
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
      params = sim_params,
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
    cells = 50,
    sub_cells = 9,
    nbrhood_radius = 4,
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






