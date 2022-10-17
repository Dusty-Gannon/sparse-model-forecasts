# libraries
library(here)
library(rjags)
library(coda)
devtools::load_all()

post_mode <- function(x){
  id <- which.max(density(x)$y)
  return(density(x)$x[id])
}

# simulating from a GLARMA(1,1) process

# set parameters
# time series length
n <- 300

# intercept
alpha <- 1

# nz slopes
beta_nz <- c(1, -1, 0.5)

# combine
beta <- c(
  alpha,
  beta_nz,
  rep(0, 46)
)

# AR and MA parameters
phi <- 0.6
theta <- 0.5

# model matrix
X <- matrix(
  rnorm(n = n * (length(beta) - 1)),
  nrow = n
)
X <- cbind(
  rep(1, n),
  X
)

# initiate the process
y <- vector(mode = "double", length = n)
y_star <- vector(mode = "double", length = n)
mu <- vector(mode = "double", length = n)
z <- vector(mode = "double", length = n)
y[1] <- rpois(1, lambda = exp(alpha))
y_star[1] <- max(y[1], 0.1)
mu[1] <- y_star[1]
z[1] <- 0

# continue the process

for(t in 2:n){

  # arma term
  z_t <- phi * (log(y_star[t-1]) - X[t, ] %*% beta) + theta * log(y_star[t-1]/mu[t-1])
  mu[t] <- exp(X[t, ] %*% beta + z_t)
  y[t] <- rpois(1, lambda = mu[t])
  y_star[t] <- max(y[t], 0.1)

}


# compile data to fit model in JAGS
datlist <- list(
  N = n,
  P0 = length(alpha),
  P = length(beta) - length(alpha),
  y = y,
  y_star = y_star,
  X_alpha = as.matrix(X[, 1]),
  X_beta = X[, -1],
  nu0 = 0.01,
  a = c(1, 4)
)

# # initial values
# nchains <- 3
#
# inits <- vector(mode = "list")
# for(i in 1:nchains){
#   inits[[i]] <- list(
#     alpha = rnorm(length(alpha)),
#     beta = rnorm(length(beta) - length(alpha)),
#     phi = runif(1),
#     theta = runif(1),
#     rho = runif(1),
#     tau_inv = rgamma(1, shape = 3, rate = 1)
#   )
# }
#
# test <- jags.model(
#   here("JAGS/Pois_GLARMA-1-1_SpSl.txt"),
#   data = datlist,
#   inits = inits,
#   n.chains = nchains
# )
#
# update(test, n.iter=5000)
# samples <- coda.samples(test, variable.names=c("alpha", "beta", "phi", "theta", "rho"),
#                         n.iter=2000)
#
# samps_all <- rbind(samples[[1]], samples[[2]], samples[[3]])


# stan model fitting
  pois_glarma11_spsl <- stan_model(here("Stan/Pois_GLARMA-1-1_SpSl.stan"))

  mfit <- sampling(
    pois_glarma11_spsl,
    data = datlist,
    cores = 3,
    chains = 3,
    control = list(adapt_delta = 0.9)
  )

# which variables would get included?
  rho_post <- rstan::extract(mfit, pars = "rho")$rho
  post_modes_rho <- apply(rho_post, 2, post_mode)

  sum(post_modes_rho > 0.5)





