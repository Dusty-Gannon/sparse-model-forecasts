# testing arma model with stan

n <- 200

beta <- c(0.5, 0.6, rep(0, 20))

X <- matrix(
  rnorm(n * length(beta)),
  nrow = n
)

y <- as.double(1 + X %*% beta + arima.sim(
  model = list(ar = 0.5),
  n = n,
  sd = 0.5
))

# testing differences in plug-in versus scaling
# global shrinkage parameter, tau. mod defines a model
# where tau0 scales with the error variance
mod <- rstan::stan_model("Stan/AR-p_err3_FHS2_DG.stan")
mod_og <- rstan::stan_model("Stan/AR-p_err3_FHS_DG.stan")

dat <- list(
  N = n - 10,
  y = y[1:(n - 10)],
  P = ncol(X) + 1,
  P_0 = 1,
  p = 15,
  tau0_phi = tau0(y, 2, 15, 200, "gaussian"),
  slab_scl_phi = 0.5,
  slab_df_phi = 6,
  tau0_beta = tau0(y, 4, 21, 200, "gaussian"),
  slab_scl_beta = 1,
  slab_df_beta = 6,
  N_new = 10,
  X = cbind(1, X[1:(n - 10), ]),
  X_new = cbind(1, X[(n - 9):n, ])
)

fit <- rstan::sampling(
  mod,
  dat,
  cores = 4,
  iter = 4000,
  control = list(adapt_delta = 0.95)
)

fit2 <- rstan::sampling(
  mod_og,
  dat,
  cores = 4,
  iter = 3000,
  control = list(adapt_delta = 0.99)
)

# comp <- arima(y, order = c(2, 0, 2))


# create a forecast plot
library(ggplot2)

y_rep <- rstan::extract(fit2)$y_rep

df_plot <- data.frame(
  time = 1:n,
  y = y,
  y_pred = colMeans(y_rep, na.rm = TRUE),
  low = apply(y_rep, 2, quantile, probs = 0.025, na.rm = TRUE),
  high = apply(y_rep, 2, quantile, probs = 0.975, na.rm = TRUE)
)

ggplot(df_plot, aes(x = time)) +
  geom_line(aes(y = y_pred), color = "blue") +
  geom_ribbon(aes(ymin = low, ymax = high), fill = "blue", alpha = 0.2) +
  geom_line(aes(y = y), color = "black") +
  theme_minimal()

# check stationarity
phi_post <- rstan::extract(fit, pars = "phi")$phi
arima.sim(model = list(ar = phi_post[200, ]), n = 20)

# compare to autoarima
test <- forecast::auto.arima(y[1:(n - 10)], xreg = dat$X)

fcs <- forecast::forecast(test, h = 10, xreg = dat$X_new) |>
  as.data.frame()

df_test <- data.frame(
  time = 1:n,
  y = y,
  y_pred = c(fitted(test), fcs$`Point Forecast`),
  low = c(fitted(test), fcs$`Lo 95`),
  high = c(fitted(test), fcs$`Hi 95`)
)

ggplot(df_test, aes(x = time)) +
  geom_line(aes(y = y_pred), color = "blue") +
  geom_ribbon(aes(ymin = low, ymax = high), fill = "blue", alpha = 0.2) +
  geom_line(aes(y = y), color = "black") +
  theme_minimal()



