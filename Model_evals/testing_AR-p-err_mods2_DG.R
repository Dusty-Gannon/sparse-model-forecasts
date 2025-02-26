# testing arma model with stan

n <- 200

y <- 1:n * 0.01 + arima.sim(
  model = list(ar = 0.7),
  n = n,
  sd = 0.5
)

mod <- rstan::stan_model("Stan/AR-p_err_stationary.stan")
mod2 <- rstan::stan_model("Stan/AR-p_err_stationary_FHS.stan")

dat <- list(
  N = n - 10,
  y = y[1:(n - 10)],
  p = 2,
  tau0_phi = 0.0001,
  slab_scl_phi = 0.5,
  slab_df_phi = 6,
  N_new = 10,
  X = cbind(1, 1:(n - 10)),
  X_new = cbind(1, (n - 9):n),
  K = 2
)

fit <- rstan::sampling(
  mod,
  dat,
  cores = 4,
  iter = 4000,
  control = list(adapt_delta = 0.95)
)

fit2 <- rstan::sampling(
  mod2,
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



