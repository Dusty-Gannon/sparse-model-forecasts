
library(rstan)
expose_stan_functions(stanc(file = "Stan/DLM_functions.stan"))

y <- arima.sim(
  model = list(
    ar = c(0.5, 0.1),
    ma = 0.1
  ),
  n = 100
)

dat <- list(
  N = length(y),
  p = 2,
  q = 1,
  m_0 = 0,
  C_0 = 5,
  y = as.double(y)
)

test <- stan_model("Stan/ARMA-p-q-DLM.stan")
dlm2 <- stan_model("Stan/ARMA-p-q-DLM_v2.stan")

fit <- sampling(
  dlm2,
  data = dat,
  chains = 3,
  cores = 3,
  control = list(adapt_delta = 0.9)
)


### Testing reparameterization

arp_v2 <- stan_model("Stan/AR-p_v2.stan")
arp_v3 <- stan_model("Stan/AR-p_v3.stan")

y_ar <- arima.sim(
  model = list(
    ar = c(0.6, -0.1, 0.1)
  ),
  n = 100,
  sd = 3
)

datv2 <- list(
  N = length(y_ar),
  p = 3,
  y = y_ar
)

fitv2 <- sampling(arp_v2, data = datv2)
fitv3 <- sampling(arp_v3, data = datv2)


### Testing AR error model

#### a simple test

n <- 100

beta <- c(1, 0.5)

X <- cbind(
  rep(1, n),
  rnorm(n, sd = 2)
)

mu <- X %*% beta

y <- as.double(mu + arima.sim(list(ar = c(0.5, 0.1)), n = n))

dat_ar_err <- list(
  N = n,
  p = 2,
  P = ncol(X),
  y = y,
  X = X
)

ar_err <- stan_model("Stan/AR-p_err3_DG.stan")

ar_err_fit <- sampling(
  ar_err,
  data = dat_ar_err,
  cores = 4
)


#### extend to seasonal model and compare to flat priors

ar_err_gauss <- stan_model("Stan/AR-p_err3_Gauss_DG.stan")
ar_err_flat <- stan_model("Stan/AR-p_err3_Flat_DG.stan")

n_tot <- 365 * 3

mu <- 0.3 * 1:n_tot/365 + 0.01 * (1:n_tot%%365) +
  0.02 * (1:n_tot%%(29*3)) + 0.05 * (1:n_tot%%12)

y <- ts(mu + arima.sim(list(ar = c(0.5, 0.1)), n = n_tot), frequency = 365)


X_tot <- cbind(
  rep(1, n_tot),
  1:n_tot / 365,
  forecast::fourier(y, K = 180)
)

n <- 365 * 2
X <- X_tot[1:n, ]
X_new <- X_tot[(n + 1):n_tot, ]

dat <- list(
  N = n,
  P = ncol(X),
  p = 10,
  y = as.double(y)[1:n],
  X = X,
  N_new = nrow(X_new),
  X_new = X_new
)

dat_gauss <- c(
  dat,
  list(P_0 = 1)
)

gauss_fit <- sampling(
  ar_err_gauss,
  data = dat_gauss,
  cores = 4
)

flat_fit <- sampling(
  ar_err_flat,
  data = dat,
  cores = 4
)

fits <- list(
  gauss_fit = gauss_fit,
  flat_fit = flat_fit
)


forecast_plot <- function(fit, y_full, fill, freq, ylim){

  y_rep <- rstan::extract(fit, pars = "y_rep")$y_rep

  dat_plot <- as_tibble(
    sapply(
      c(0.025, 0.1, 0.9, 0.975),
      FUN = function(x){
        apply(y_rep, 2, quantile, probs = x)
      }
    )
  )
  names(dat_plot) <- c("low", "mlow", "mhigh", "high")

  dat_plot <- dat_plot %>% mutate(
    y = as.double(y),
    t = 1:length(y) / freq,
  )

  ggplot(data = dat_plot, aes(x = t)) +
    geom_ribbon(aes(ymin = low, ymax = high), fill = fill, alpha = 0.3) +
    geom_ribbon(aes(ymin = mlow, ymax = mhigh), fill = fill, alpha = 0.6) +
    geom_line(aes(y = y), linewidth = 0.2) +
    theme_classic() +
    ylim(ylim)

}

plots <- map(
  fits,
  ~ forecast_plot(
    .x,
    y_full = y,
    fill = "steelblue",
    freq = 365,
    ylim = c(-5,10)
  )
)

plots[[1]] / plots[[2]]

rmse_gauss <- RMSE_bayes(
  as.double(y)[(n + 1):n_tot],
  ppreds = rstan::extract(gauss_fit, "y_rep")$y_rep[,(n + 1):n_tot]
)
rmse_flat <- RMSE_bayes(
  as.double(y)[(n + 1):n_tot],
  ppreds = rstan::extract(flat_fit, "y_rep")$y_rep[,(n + 1):n_tot]
)

par(mfrow = c(1,2))
hist(rmse_flat, breaks = 100, xlim = c(0, 3.5))
abline(v = mean(rmse_flat), col = "blue")
hist(rmse_gauss, breaks = 100, xlim = c(0, 3.5))
abline(v = mean(rmse_gauss), col = "blue")


