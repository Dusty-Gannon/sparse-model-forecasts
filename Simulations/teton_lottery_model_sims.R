# document setup
library(here)

# load user-defined functions (teton way, not devtools way)
src_files <- list.files(here("R/"), pattern = "*.R", full.names = T)
sapply(src_files, source, .GlobalEnv)

# setup parameters

# lattice size along one margin
cells <- 50

# neighborhood size
nbrhood_radius <- 2

# time steps
steps <- 800

# species
S <- 2 # number of strong competitors
s <- 47 # number weak competitors
nsp <- S + s + 1

# vector of intra-specific competitive effects
#  alpha <- rgamma(nsp, shape = 3, rate = 1)
alpha <- c(
  sort(runif(S+1, min = 3.5, max = 4), decreasing = T),
  sort(runif(s, min = 4, max = 5), decreasing = T)
)

# competition matrix
A_mat <- comp_matrix(
  n_sp = nsp,
  alpha = alpha,
  rho = c(0.01, 0.6),
  num_ngs = S, num_regs = 2
)

# average per-capita fecundity
# lambda <- c(
#   rgamma(S + 1, shape = 40, rate = 1),
#   rgamma(s, shape = 10, rate = 1)
# )
lambda_max <- c(
  sort(runif(S + 1, min = 45, max = 50)),
  sort(runif(s, min = 10, max = 15))
)
sp_optims <- runif(nsp, min = -0.2, max = 0.2)

# probability an individual dies in a given time point for each species
# Pr_death <- c(
#   rep(1, 1 + S),
#   runif(s, min = 0.1, max = 0.5)
# )
Pr_death <- c(
  rep(1, S + 1),
  runif(s, min = 0.05, max = 0.5)
)

# dispersal rates for each species
# disp_rate <- c(
#   runif(S + 1),
#   runif(s, min = 0, max = 3)
# )
disp_rate <- c(
  c(0.1, sort(runif(S, 0.2, 0.3))),
  sort(runif(s, 0.5, 1.1))
)

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

# track environment
sigma_env <- 0.1    # variation in the environment
phi <- 0.1          # autocorrelation across time steps
env <- vector(mode = "double", length = steps)
env[1] <- rnorm(1, sd = sigma_env/sqrt(1 - phi^2))


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
  t = rep(1:steps, nsp),
  species = rep(1:nsp, each = steps)
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

# how many species coexist?
sum(cover_df[cover_df$t == steps, ]$cover > 0)





