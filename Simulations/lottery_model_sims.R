# document setup
library(here)
library(ggplot2)

devtools::load_all()

# setup parameters

  # lattice size along one margin
  cells <- 10

  # neighborhood size
  nbrhood_radius <- 1

  # time steps
  steps <- 100

  # species
  S <- 1 # number of strong competitors
  s <- 1 # number weak competitors
  nsp <- S + s + 1

  # vector of intra-specific competitive effects
  alpha <- runif(nsp, min = 0.2, max = 1)

  # competition matrix
  A_mat <- comp_matrix(nsp, alpha = alpha, rho = 0.05, num_ngs = S)

  # average per-capita fecundity
  lambda <- rgamma(nsp, shape = 10, rate = 1)

  # probability an individual dies in a given time point for each species
  Pr_death <- c(0.3, 0.2, 0.1)

  # initialize lattice
  X <- matrix(
    data = sample(1:nsp, cells * cells, replace = T),
    nrow = cells, ncol = cells
  )

  # count the neighbors for each cell given the radius and
  #  lattice size
  n_neighbors <- count_neighbors(M = cells, J = cells, r = nbrhood_radius)

  # track evolution of community
  ts_all <- vector(mode = "list", length = steps)
  ts_all[[1]] <- X

  # step through the evolution of the community
  for(t in 1:(steps-1)){

    # count neighbors of each species around each cell
    nbrs_t <- kernel_count(ts_all[[t]], r = nbrhood_radius, sp_list = 1:nsp)

    # get fecundity matrix by applying the fecundidty_ll() command to each row of X
    Fecundity <- t(sapply(
      1:cells,
      FUN = function(x, lattice, lambda_t, alpha, nbrhood, n){
        fecundity_ll(
          foc_sp = lattice[x, ],
          lambda_t = lambda_t,
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
      rate = c(1,0.5,3)
    )

    # kill off and replace some adults with seedlings
    ts_all[[t+1]] <- die_replace(
      ts_all[[t]],
      prob_death = Pr_death,
      seed_rain = sr_t
    )

  }

# get percent cover dataframe
  cover_df <- data.frame(
    t = rep(1:steps, nsp),
    species = rep(1:nsp, each = steps),
    cover = c(
      sapply(ts_all, FUN = function(x){sum(x == 1)}),
      sapply(ts_all, FUN = function(x){sum(x == 2)}),
      sapply(ts_all, FUN = function(x){sum(x == 3)})
    )
  )
  cover_df$species <- as.factor(cover_df$species)

  ggplot(data = cover_df, aes(x = t, y = cover, color = species))+
    geom_line()+
    theme_classic()+
    scale_color_manual(values = c("blue", "red", "black"))









