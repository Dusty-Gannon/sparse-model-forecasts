

library(rstan)
devtools::load_all()

# generate some parameters

params <- generate_sim_params_dist(
  nsp = 10, steps = 50,
  num_ngs = 2, rho = 0.1,
  spatial = T, n_sites = 2
)

# simulate a community
N <- do.call(ricker_spts_lnorm, params)$N

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

dat <- list(
  S = params$nsp,
  N = params$steps - 1,
  K = params$n_sites,
  Ntm1 = Ntm1,
  y = y
)

gp_comm <- stan_model("Stan/regGP_comm_growth.stan")

test <- sampling(gp_comm, data = dat, chains = 1)
