
library(rstan)
library(tidyverse)
library(patchwork)

#### testing seasonal model and compare priors

ar_err_gauss <- stan_model("Stan/AR-p_err3_Gauss_DG.stan")
ar_err_flat <- stan_model("Stan/AR-p_err3_Flat_DG.stan")
ar_err_hs <- stan_model("Stan/AR-p_err3_FHS_DG.stan")

# create time series of 2.5 years of daily data
n_tot <- 365 * 2 + 180
set.seed(8968)

# create seasonally fluctuating mean
mu2 <-  1.5 * 1:n_tot/365 + 2 * cos(pi * 1:n_tot / 365) + sin(pi * 1:n_tot / 30) +
  0.5 * sin(pi * 1:n_tot / 10)

# add the AR-2 errors
y <- ts(mu2 + arima.sim(list(ar = c(0.5, 0.1)), n = n_tot), frequency = 365)

# create the model matrix using K = 100 fourier components
X_tot <- cbind(
  rep(1, n_tot),                      # intercept
  1:n_tot / 365,                      # trend
  forecast::fourier(y, K = 100)       # sparse seasonality terms
)

# fit the model to the first two years
n <- 365 * 2
X <- X_tot[1:n, ]
X_new <- X_tot[(n + 1):n_tot, ]

# compile data
dat_flat <- list(
  N = n,
  P = ncol(X),
  p = 15,
  y = as.double(y)[1:n],
  X = X,
  N_new = nrow(X_new),
  X_new = X_new
)

# don't regularize the trend for the
# gaussian fit
dat_gauss <- c(
  dat_flat,
  list(P_0 = 2)
)

# add in HS-specific inputs
dat_hs <- c(
  dat_gauss,
  list(
    tau0_beta = 0.01,
    slab_scl_beta = 0.5,
    slab_df_beta = 4,
    tau0_phi = 0.001,
    slab_scl_phi = 0.5,
    slab_df_phi = 4
  )
)

# fit the models
gauss_fit <- sampling(
  ar_err_gauss,
  data = dat_gauss,
  cores = 4
)

flat_fit <- sampling(
  ar_err_flat,
  data = dat_flat,
  cores = 4
)

hs_fit <- sampling(
  ar_err_hs,
  dat = dat_hs,
  cores = 4,
  control = list(adapt_delta = 0.9)
)

# put fits into a list
fits <- list(
  flat_fit = flat_fit,
  gauss_fit = gauss_fit,
  hs_fit = hs_fit
)

# create plotting function
forecast_plot <- function(fit, y_full, fill, freq, ylim, title, rmse = NULL, xlim_rmse = NULL, ylim_rmse = NULL){

  library(patchwork)
  library(dplyr)
  library(ggplot2)

  y_rep <- rstan::extract(fit, pars = "y_rep")$y_rep

  dat_plot <- dplyr::as_tibble(
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

  ts_plot <- ggplot(data = dat_plot, aes(x = t)) +
    geom_ribbon(aes(ymin = low, ymax = high), fill = fill, alpha = 0.3) +
    geom_ribbon(aes(ymin = mlow, ymax = mhigh), fill = fill, alpha = 0.6) +
    geom_line(aes(y = y), linewidth = 0.2) +
    theme_classic() +
    ylim(ylim) +
    ggtitle(title) +
    xlab("") + ylab("")

  if(is.null(rmse)){
    return(ts_plot)
  } else {
    rmse_df <- data.frame(
      x = rmse
    )
    rmse_plot <- ggplot(data = rmse_df, aes(x = x)) +
      geom_density(fill = fill) +
      theme_classic() +
      xlab("") + ylab("") +
      xlim(xlim_rmse) +
      ylim(ylim_rmse)
    return(
      ts_plot + rmse_plot +
        plot_layout(widths = c(2, 1))
    )
  }

}

# compute forecast rmse for each fit
rmse <- list(
  rmse_flat = RMSE_bayes(
    as.double(y)[(n + 1):n_tot],
    ppreds = rstan::extract(fits$flat_fit, "y_rep")$y_rep[,(n + 1):n_tot]
  ),
  rmse_gauss = RMSE_bayes(
    as.double(y)[(n + 1):n_tot],
    ppreds = rstan::extract(fits$gauss_fit, "y_rep")$y_rep[,(n + 1):n_tot]
  ),
  rmse_hs = RMSE_bayes(
    as.double(y)[(n + 1):n_tot],
    ppreds = rstan::extract(fits$hs_fit, "y_rep")$y_rep[, (n + 1):n_tot]
  )
)

# construct the plots
plots <- mapply(
  FUN = forecast_plot,
  fit = fits,
  title = c("Flat", "Gaussian", "Horseshoe"),
  rmse = rmse,
  MoreArgs = list(
    ylim = c(-5, 10),
    xlim_rmse = c(0, 8.5),
    ylim_rmse = c(0, 1.5),
    y_full = y,
    fill = "brown",
    freq = 365
  ),
  SIMPLIFY = F
)

# organize the plots
xlab_plot1 <- ggplot(data.frame(l = "Time (years)", x = 1, y = 1)) +
  geom_text(aes(x, y, label = l)) +
  theme_void() +
  coord_cartesian(clip = "off")
xlab_plot2 <- ggplot(data.frame(l = "Forecast RMSE", x = 1, y = 1)) +
  geom_text(aes(x, y, label = l)) +
  theme_void() +
  coord_cartesian(clip = "off")

p4 <- xlab_plot1 + xlab_plot2 +
  plot_layout(widths = c(2, 1))

plots[[1]] / plots[[2]] / plots[[3]] / p4 +
  plot_layout(heights = c(2, 2, 2, 1))




