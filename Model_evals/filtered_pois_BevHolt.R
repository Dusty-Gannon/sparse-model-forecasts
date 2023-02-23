# libraries
library(here)
library(rjags)
library(coda)
library(ggplot2)
devtools::load_all()

post_mode <- function(x){
  id <- which.max(density(x)$y)
  return(density(x)$x[id])
}

# simulating a Beverton-Holt model that then gets filtered by survival

# set parameters
# time series length
n <- 100
S <- 20

# abundance matrix
N <- matrix(data = 0, nrow = S, ncol = n)
N[, 1] <- rpois(S, lambda = 20)

# intra-specific competition
alphas <- runif(S, min = 0.05, max = 0.1)

# competition matrix
A_mat <- comp_matrix2(S, rho = 0, alpha = alphas, num_ngs = 2)

# max intrinsic growth
lambdas <- runif(S, min = 50, max = 100)

# environmental preferences
nu <- runif(S, min = -1, max = 1)

# initiate the environment
env <- vector(mode = "double", length = n)
env[1] <- 0

# data frame to record the simulations
df_sims <- data.frame(
  t = rep(1:n, each = S),
  species = rep(1:S, n),
  abundance = c(N[, 1], rep(0, (n - 1) * S))
)

# survival probabilities
p_surv <- runif(S, min = 0.1, max = 0.2)

# fecundity matrix
fec <- matrix(data = 0, nrow = S, ncol = (n - 1))

# step the process forward
for(t in 2:n){

  # realized growth rates
  lambda_t <- lambdas * exp(-0.5 * (env[t - 1] - nu)^2)

  # fecundity
  fec_t <- rpois(S, (N[,t - 1] * lambda_t / (1 + A_mat %*% N[, t - 1])))
  fec[, t - 1] <- fec_t

  # Abundance in the next year
  N[, t] <- rbinom(S, size = fec_t, prob = p_surv)

  # add to dataframe
  df_sims$abundance[(S * (t - 1) + 1):(S * (t-1) + S)] <- as.double(N[, t])

  # adjust environment for the next year
  env[t] <- 0.5 * env[t - 1] + rnorm(n = 1, sd = 0.5)

}

# plot to check the data
df_sims$species <- as.factor(df_sims$species)
ggplot(data = df_sims, aes(x = t, y = abundance, color = species))+
  geom_line()


# compile data to fit model in JAGS
X_alpha <- cbind(
  rep(1, n),
  env,
  env^2
)
X_beta <- t(N)

datlist <- list(
  N = n - 1,
  P0 = ncol(X_alpha),
  P = ncol(X_beta),
  X_alpha = X_alpha[-n, ],
  X_beta = X_beta[-n, ],
  fec = as.integer(fec[8, ]),
  vm = as.double(N[8, -n]),
  vp = as.double(N[8, -1]),
  slab_scl = 1,
  slab_df = 20,
  tau0 = 0.001
)

# initial values
nchains <- 3

inits <- vector(mode = "list")
for(i in 1:nchains){
  inits[[i]] <- list(
    alpha_std = rnorm(ncol(X_alpha)),
    beta_std = rnorm(ncol(X_beta)),
    local_scl = rexp(ncol(X_beta)),
    tau_std = rexp(1),
    p = runif(1)
  )
}

test <- jags.model(
  here("JAGS/filtered_Pois_BevHolt_FHS.jags"),
  data = datlist,
  inits = inits,
  n.chains = nchains
)

update(test, n.iter=20000)
samples <- coda.samples(test, variable.names=c("alpha", "beta", "p"),
                        n.iter=5000, thin = 5)

samps_all <- rbind(samples[[1]], samples[[2]], samples[[3]])








