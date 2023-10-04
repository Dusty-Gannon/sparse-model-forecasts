

library(rstan)
devtools::load_all()

# generate some parameters

params <- generate_sim_params_dist(
  nsp = 10, steps = 50,
  num_ngs = 3, rho = 0.1,
  spatial = T, n_sites = 5,
  dist_prob = 0.1, dist_int = 0.6,
  prop_cdist = 0.9
)

# simulate a community
sim <- do.call(ricker_spts_lnorm, params)

# for ease, add small value to the 0s
N <- sim$N
N[N == 0] <- 0.5

# construct data for stan model
Ntm1 <- matrix(
  nrow = params$nsp * params$n_sites,
  ncol = params$steps - 1
)

for(t in 1:(params$steps - 1)){
  Ntm1[, t] <- as.vector(N[, t, ])
}

y <- matrix(
  nrow = nrow(Ntm1),
  ncol = ncol(Ntm1)
)

for(t in 2:params$steps){
  n_t <- as.vector(N[, t, ])
  n_tm1 <- as.vector(N[, t - 1, ])
  y[, t - 1] = log(n_t / n_tm1)
}

# construct vectors of disturbance indicators
X_beta <- Reduce(rbind, sim$disturbances)[, -1]

dat <- list(
  S = params$nsp,
  N = params$steps - 1,
  K = params$n_sites,
  Ntm1 = Ntm1,
  y = y,
  X_beta = X_beta,
  tau0 = 0.0001,
  slab_df = 6,
  slab_scl = 0.1
)

gp_comm <- stan_model("Stan/regGP_comm_growth_FHS.stan")

test <- sampling(gp_comm, data = dat, chains = 1)

##########
# Testing version 2

params <- generate_sim_params_dist(
  nsp = 20, steps = 100,
  num_ngs = 3, rho = 0.01,
  dist_prob = 0.2, dist_int = 0.6,
  prop_cdist = 0.5
)

# simulate a community
sim <- with(params, {
  ricker_ts_lnorm(
    N_0 = N_0,
    lambdas = lambdas,
    A_mat,
    sigmas,
    steps,
    dist_prob,
    dist_int,
    prop_cdist
  )
})

# for ease, add small value to the 0s
N <- sim$N
N[N == 0] <- 0.5

# construct data for stan model
Ntm1 <- N[, -params$steps]

ymat <- matrix(
  nrow = nrow(Ntm1),
  ncol = ncol(Ntm1)
)

for(t in 2:params$steps){
  n_t <- as.vector(N[, t])
  n_tm1 <- as.vector(N[, t - 1])
  ymat[, t - 1] = log(n_t / n_tm1)
}

y <- as.vector(ymat)

# construct disturbance covariate
X_beta <- Reduce(cbind, sim$disturbances)[, -1]
X_beta[X_beta > 0] <- 1

dat2 <- list(
  S = params$nsp,
  N = params$steps - 1,
  Ntm1 = Ntm1,
  y = y,
  P = 1,
  X_beta = matrix(
    data = as.vector(X_beta),
    ncol = 1
  ),
  sp_id = rep(1:params$nsp, params$steps - 1),
  tau0 = 0.0001,
  slab_df = 6,
  slab_scl = 0.1
)


gp_comm2 <- stan_model("Stan/regGP_comm_growth_FHS_v2.stan")

test <- sampling(gp_comm2, data = dat2, chains = 1)







