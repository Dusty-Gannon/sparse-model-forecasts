
library(rstan)
library(here)
devtools::load_all()

# define sample size
n <- 100

# define sparse coefficient vector
beta <- c(0, 0.5, 1, rep(0, 100))

# define model matrix
X <- cbind(
  rep(1, n),
  matrix(
    rnorm(n * length(beta) - n, sd = 2),
    nrow = n
  )
)

# construct response
y <- as.double(X %*% beta) + rnorm(n)

# compile model
hs_mod <- stan_model(here("Stan/sparse_reg_FHS.stan"))

# ready the data inputs
dat <- list(
  N = n,
  P = ncol(X),
  P0 = 1, # number of coefs getting weakly informative priors
  y = y,
  X = X,
  N_new = 0, # no out-of-sample predictions
  X_new = matrix(data = 1, nrow = 0, ncol = ncol(X)),
  tau0 = tau0(
    y, m0 = 5,
    M = ncol(X) - 1,
    N = n,
    fam = "gaussian"
  ),
  slab_scl = 1, # guess for scale of "large" effects
  slab_df = 6
)

# fit model
mfit <- sampling(hs_mod, data = dat, cores = 4)

# plot results for beta
stan_plot(mfit, pars = "beta")
