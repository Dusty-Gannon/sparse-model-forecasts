###############################################################
# Conceptual figure: bias-variance tradeoff
#
# True model: Gaussian curve f(x) = exp(-0.5 * (x/1.5)^2)
# Three noisy realizations are fit with:
#   a) a full (degree n-1) interpolating polynomial
#   b) the same polynomial basis regularized with a horseshoe prior
###############################################################

library(rstan)
library(here)
library(tidyverse)
library(patchwork)
devtools::load_all(here())

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

set.seed(260731)

# ---- Setup -------------------------------------------------------

true_f  <- function(x) {
  exp(-0.5 * (x / 1.5)^2)
}

n_obs   <- 20
degree  <- 15
sigma_e <- 0.15

x_obs  <- seq(-3, 3, length.out = n_obs)
x_pred <- seq(-3, 3, length.out = 400)

# three dataset colors
pal <- c("#c03728", "#028090", "#e6a000")

# ---- Generate 3 realizations ------------------------------------

datasets <- lapply(1:3, function(i) {
  data.frame(
    x   = x_obs,
    y   = true_f(x_obs) + rnorm(n_obs, 0, sigma_e),
    rep = factor(i)
  )
})

obs_df  <- bind_rows(datasets)
true_df <- data.frame(x = x_pred, y = true_f(x_pred))

# ---- Panel a: full interpolating polynomial ----

poly_df <- bind_rows(lapply(1:3, function(i) {
  d   <- datasets[[i]]
  fit <- lm(y ~ poly(x, degree), data = d)
  data.frame(
    x   = x_pred,
    y   = predict(fit, newdata = data.frame(x = x_pred)),
    rep = factor(i)
  )
}))

# ---- Panel b: horseshoe-regularized polynomial ------------------

hs_mod <- stan_model(here("Stan/sparse_reg_FHS.stan"))

hs_df <- bind_rows(lapply(1:3, function(i) {
  d  <- datasets[[i]]

  # orthogonal polynomial basis anchored to training x
  pb <- poly(d$x, degree)
  X_train <- cbind(1, as.matrix(pb))
  X_new   <- cbind(1, predict(pb, newdata = x_pred))

  datlist <- list(
    N        = n_obs,
    P0       = 1,                   # intercept unshrunk
    P        = ncol(X_train),
    y        = d$y,
    X        = X_train,
    tau0     = tau0(d$y, m0 = 2, M = degree, N = n_obs, fam = "gaussian"),
    slab_scl = 2,
    slab_df  = 4,
    N_new    = length(x_pred),
    X_new    = X_new
  )

  fit   <- sampling(
    hs_mod, data = datlist,
    chains = 4, cores = 4,
    iter = 2000, warmup = 1000,
    refresh = 0
  )

  y_rep <- rstan::extract(fit, pars = "y_rep")$y_rep
  y_new <- y_rep[, (n_obs + 1):(n_obs + length(x_pred))]

  data.frame(
    x   = x_pred,
    y   = colMeans(y_new),
    lo  = apply(y_new, 2, quantile, 0.1),
    hi  = apply(y_new, 2, quantile, 0.9),
    rep = factor(i)
  )
}))

# ---- Shared theme -----

fig_theme <- theme_classic() +
  theme(
    legend.position = "none",
    axis.title      = element_text(size = 8),
    axis.text       = element_text(size = 7),
    plot.title      = element_text(size = 9),
    plot.margin     = margin(8, 12, 8, 8)
  )

y_lim <- c(-0.5, 1.5)

# ---- Figure ----

## ---- Panel a ----

p1 <- ggplot() +
  geom_line(
    data = true_df,
    aes(x = x, y = y),
    color = "grey40", linewidth = 0.8, linetype = "dashed"
  ) +
  geom_point(
    data = obs_df,
    aes(x = x, y = y, color = rep),
    size = 1, alpha = 0.7
  ) +
  geom_line(
    data  = poly_df,
    aes(x = x, y = y, color = rep),
    linewidth = 0.4
  ) +
  scale_color_manual(values = pal) +
  coord_cartesian(ylim = y_lim) +
  labs(x = "x", y = "y",
       title = "a) Low bias,\nhigh variance") +
  fig_theme

## ---- Panel b ----

p2 <- ggplot() +
  geom_line(data  = true_df,
            aes(x = x, y = y),
            color = "grey40", linewidth = 0.8, linetype = "dashed") +
  geom_point(
    data = obs_df,
    aes(x, y, color = rep),
    size = 1,
    alpha = 0.7
  ) +
  geom_smooth(
    data = obs_df,
    aes(x, y, color = rep),
    linewidth = 0.4,
    se = F, method = "lm"
  ) +
  scale_color_manual(values = pal) +
  coord_cartesian(ylim = y_lim) +
  labs(x = "x", y = "y",
       title = "b) High bias,\nlow variance") +
  fig_theme

## ---- Panel c ----

p3 <- ggplot() +
  geom_line(
    data = true_df,
    aes(x = x, y = y),
    color = "grey40",
    linewidth = 0.8,
    linetype = "dashed"
  ) +
  geom_point(
    data = obs_df,
    aes(x = x, y = y, color = rep),
    size = 1,
    alpha = 0.7
  ) +
  geom_smooth(
    data = hs_df,
    aes(x = x, y = y, color = rep),
    se = F,
    linewidth = 0.4
  ) +
  scale_color_manual(values = pal) +
  coord_cartesian(ylim = y_lim) +
  labs(x = "x", y = "y",
       title = "c) Balancing bias\nand variance") +
  fig_theme


# ---- Combine and save ----

fig <- p1 + p2 + p3

ggsave(
  here("Figures/bias_variance_tradeoff.pdf"),
  plot = fig,
  width = 6.5,
  height = 2,
  units = "in",
  dpi = 400,
  device = "pdf"
)





# ---- Priors intuition ----
d  <- datasets[[1]]

# orthogonal polynomial basis anchored to training x
pb <- poly(d$x, degree)
X_train <- cbind(1, as.matrix(pb))
X_new   <- cbind(1, predict(pb, newdata = x_pred))

datlist <- list(
  N        = n_obs,
  P0       = 1,                   # intercept unshrunk
  P        = ncol(X_train),
  y        = d$y,
  X        = X_train,
  tau0     = tau0(d$y, m0 = 2, M = degree, N = n_obs, fam = "gaussian"),
  slab_scl = 2,
  slab_df  = 4,
  N_new    = length(x_pred),
  X_new    = X_new
)

fit   <- sampling(
  hs_mod, data = datlist,
  chains = 4, cores = 4,
  iter = 2000, warmup = 1000,
  refresh = 0
)


post_draws <- rstan::extract(fit, pars = c("tau", "sigma", "local_scale_tilde"))

kappa_draws <- apply(
  post_draws$local_scale_tilde,
  MARGIN = 2,
  FUN = function(x){
    mapply(
      function(x, s, t, n = nrow(X_train)) {
        1 / (1 + n * s^(-2) * t^2 * x^2)
      },
      x = x,
      s = post_draws$sigma,
      t = post_draws$tau
    )
  }
)

## ---- Density plots for shrinkage coefficients ----

p_kappa2 <- ggplot(data.frame(value = kappa_draws[, 2]), aes(x = value)) +
  geom_density(fill = "#c03728", color = "#c03728", alpha = 0.4, linewidth = 0.6) +
  scale_x_continuous(limits = c(0, 1), expand = c(0.01, 0)) +
  labs(x = expression(kappa[2]), y = "Density", title = "Shrinkage coefficients") +
  fig_theme

p_kappa5 <- ggplot(data.frame(value = kappa_draws[, 5]), aes(x = value)) +
  geom_density(fill = "#c03728", color = "#c03728", alpha = 0.4, linewidth = 0.6) +
  scale_x_continuous(limits = c(0, 1), expand = c(0.01, 0)) +
  labs(x = expression(kappa[5]), y = "Density") +
  fig_theme

p_kappa <- p_kappa2 + p_kappa5




## ---- Density plots for regression coefficients ----

beta_draws <- rstan::extract(fit, pars = "beta")$beta
beta_mle   <- coef(lm(d$y ~ 0 + X_train))

p_beta2 <- ggplot(data.frame(value = beta_draws[, 3]), aes(x = value)) +
  geom_density(fill = "#c03728", color = "#c03728", alpha = 0.4, linewidth = 0.6) +
  geom_vline(xintercept = beta_mle[3], linewidth = 0.7, linetype = "dashed") +
  annotate(
    "text",
    x     = beta_mle[3],
    y     = Inf,
    label = "hat(beta[MLE])",
    parse = TRUE,
    hjust = -0.3,
    vjust = 1.1,
    size  = 2.5
  ) +
  labs(x = expression(beta[2]), y = "Density", title = "Regression coefficients") +
  fig_theme

p_beta5 <- ggplot(data.frame(value = beta_draws[, 5]), aes(x = value)) +
  geom_density(fill = "#c03728", color = "#c03728", alpha = 0.4, linewidth = 0.6) +
  geom_vline(xintercept = beta_mle[5], linewidth = 0.7, linetype = "dashed") +
  annotate(
    "text",
    x     = beta_mle[5],
    y     = Inf,
    label = "hat(beta[MLE])",
    parse = TRUE,
    hjust = -0.2,
    vjust = 1.3,
    size  = 2.5
  ) +
  labs(x = expression(beta[5]), y = "Density") +
  fig_theme

p_beta <- p_beta2 + p_beta5

## ---- Prior hierarchy + shrinkage + coefficient layout ----

library(magick)
library(grid)

system2("pdflatex", args = c(
  "-interaction=nonstopmode",
  "-output-directory", here("Figures"),
  here("Figures/prior_hierarchy.tex")
))
system2("pdfcrop", args = c(
  here("Figures/prior_hierarchy.pdf"),
  here("Figures/prior_hierarchy_crop.pdf")
))

hier_img  <- image_read_pdf(here("Figures/prior_hierarchy_crop.pdf"), density = 300)
hier_grob <- rasterGrob(as.raster(hier_img), interpolate = TRUE)
p_hier    <- wrap_elements(hier_grob) +
  plot_annotation(
    title = "Prior hierarchy",
    theme = theme(plot.title = element_text(size = 9, margin = margin(b = 4)))
  )

p_hier + (p_kappa2 / p_kappa5) + (p_beta2 / p_beta5) +
  plot_layout(widths = c(1.5, 1, 1))

ggsave(
  here::here("Figures/priors_and_shrinkage.pdf"),
  device = "pdf",
  dpi = 600,
  width = 6.5,
  height = 4,
  units = "in"
)
