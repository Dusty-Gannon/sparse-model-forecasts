
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

