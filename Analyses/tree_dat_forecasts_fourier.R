
# Load libraries
library(tidyverse)
library(forecast)
devtools::load_all()

tree_dat <- readRDS(here::here("SparseTS_prismdata/ca719_BM_Seq_1800_2012.rds"))

# average all trees for a given year
tree_dat <- tree_dat %>%
  group_by(year) %>%
  summarise(
    mean_rwi = mean(rwi, na.rm = TRUE)
  )

# set range of frequencies to try
freqs <- seq(80, 200, by = 20)

# create list of fourier term matrices
Xs <- purrr::map(
  periods,
  ~ forecast::fourier(
    ts(tree_dat$mean_rwi, frequency = .x),
    K = .x / 2
  )
)

# create list of training Xs
X_train <- purrr::map(
  Xs,
  ~ .x[1:which(tree_dat$year == 1970), ]
)

# now create list of test Xs
X_test <- purrr::map(
  Xs,
  ~ .x[(which(tree_dat$year == 1970) + 1):nrow(tree_dat), ]
)


# ---- Sparse model fitting ----

# compile stan model
sparse_mod <- rstan::stan_model(
  here::here("Stan/AR-p_err3_FHS_DG.stan")
)

base_mod <- rstan::stan_model(
  here::here("Stan/sparse_reg_FHS.stan")
)

# create list of objects that are shared across all models
# get estimate of tau0
tau_0 <- tau0(
  y = tree_dat$mean_rwi[1:which(tree_dat$year == 1970)],
  m0 = 2,
  M = 10,
  N = nrow(X_train[[1]]),
  fam = "gaussian"
)
shared_obj <- list(
  N = nrow(X_train[[1]]),
  P_0 = 1,
  p = 50,
  tau0_phi = tau_0^(-1),
  slab_scl_phi = 0.5,
  slab_df_phi = 10,
  tau0_beta = tau_0,
  slab_scl_beta = 0.8,
  slab_df_beta = 10,
  y = tree_dat$mean_rwi[1:which(tree_dat$year == 1970)],
  N_new = nrow(X_test[[1]])
)

# fit models
sparse_fits <- purrr::map(
  X_train,
  ~ rstan::sampling(
    sparse_mod,
    data = c(shared_obj, list(P = ncol(P) + 1, X = cbind(1, .x))),
    chains = 4,
    iter = 2000,
    control = list(adapt_delta = 0.95)
  )
)


test <- rstan::sampling(
  sparse_mod,
  data = c(shared_obj, list(
    P = ncol(X_train[[7]][, 1:40]) + 1,
    X = cbind(1, X_train[[7]][, 1:40]),
    X_new = cbind(1, X_test[[7]][, 1:40])
  )
  ),
  chains = 4,
  iter = 2000,
  cores = 4,
  control = list(adapt_delta = 0.97)
)

test2 <- rstan::sampling(
  base_mod,
  data = list(
    N = nrow(X_train[[7]]),
    P0 = 1,
    P = ncol(X_train[[7]]) + 1,
    tau0 = tau_0,
    slab_scl = 0.5,
    slab_df = 10,
    y = tree_dat$mean_rwi[1:which(tree_dat$year == 1970)],
    N_new = nrow(X_test[[7]]),
    P = ncol(X_train[[7]]) + 1,
    X = cbind(1, X_train[[7]]),
    X_new = cbind(1, X_test[[7]])
  ),
  chains = 4,
  iter = 2000,
  cores = 4,
  control = list(adapt_delta = 0.95)
)

# extract posterior predictive draws
y_hat <- rstan::extract(test, pars = "y_rep")$y_rep

plot_df <- data.frame(
  y = tree_dat$mean_rwi,
  y_hat = colMeans(y_hat),
  low = apply(y_hat, 2, quantile, probs = 0.025),
  high = apply(y_hat, 2, quantile, probs = 0.975),
  year = tree_dat$year
)

ggplot(plot_df, aes(x = year)) +
  geom_line(aes(y = y_hat), color = "blue", linewidth = 1) +
  geom_ribbon(aes(ymin = low, ymax = high), fill = "blue", alpha = 0.2) +
  geom_line(aes(y = y), color = "black") +
  theme_minimal()

