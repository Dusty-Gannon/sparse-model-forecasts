################################################################
# This script is used to create a conceptual figure illustrating
# covariate shift with a simple, 2 covariate example.
################################################################

# load libraries
devtools::load_all()
library(tidyverse)
library(patchwork)
set.seed(3778)

# series length
n <- 120
t <- 1:n

#### functions to create time series ####

construct_main <- function(m, h, tmin, tmax, n, sd){
  t <- seq(tmin, tmax, length.out = n)
  y <- m * t / (h + t) + arima.sim(list(ar = 0.6), n = n, sd = sd)
  return(data.frame(t = t, y = y))
}

construct_grad_2 <- function(t0, m, h, tmin, tmax, n, sd){
  t <- seq(tmin, tmax, length.out = n)
  a <- m * t0 / (h + t0)
  b <- h * m / (h + t0)^2
  y <- a + b * (t - t0) +
    arima.sim(list(ar = 0.6), n = n, sd = sd)
  return(data.frame(t = t, y = y))
}

construct_abrupt_2 <- function(t0, m, h, tmin, tmax, n, sd){
  t <- seq(tmin, tmax, length.out = n)
  a <- m * t0 / (h + t0)
  b <- h * m / (h + t0)^2
  b2 <- -h * m / (h + t0)^3
  y <- a + b * (t - t0) + b2 * (t - t0)^2 +
    arima.sim(list(ar = 0.6), n = n, sd = sd)
  return(data.frame(t = t, y = y))
}

m <- 3
h <- 0.4 * n |> round()
tmin = 1; tmax = n
t0 <- 35



# Now construct 1000 simulations from each model
x1 <- purrr::map(
  1:1000,
  ~ construct_main(m, h, tmin, tmax, n, sd = 0.1)
)

x2_grad <- purrr::map(
  1:1000,
  ~ construct_grad_2(t0, m, h, tmin, tmax, n, sd = 0.12)
)

x2_abrupt <- purrr::map(
  1:1000,
  ~ construct_abrupt_2(t0, m, h, tmin, tmax, n, sd = 0.12)
)

# construct some example time series data
ts_data <- data.frame(
  x_var = rep(c("x1", "x2", "x3"), each = n),
  x = c(x1[[1]]$y, x2_grad[[1]]$y, x2_abrupt[[1]]$y),
  time = 1:n %>% rep(3)
)


# plot the time series
(ts_plot <- ggplot(ts_data, aes(time, x, color = x_var)) +
    geom_rect(aes(
      xmin = 0, xmax = round(n * 0.75),
      ymin = 0, ymax = 3.5
    ), color = "darkgreen", fill = NA,
    linewidth = 0.8, linetype = "dashed") +
    geom_rect(aes(
      xmin = round(n * 0.75) + 1, xmax = n,
      ymin = 0, ymax = 3.5
    ), color = "darkorange4", fill = NA,
    linewidth = 1, linetype = "dashed") +
    geom_line() +
    theme_classic() +
    scale_color_manual(
      values = c("black", PNWColors::pnw_palette("Bay", 5)[c(1,5)]),
      labels = expression(x[1], x[2], x[3])
    ) +
    theme(
      legend.title = element_blank(),
      text = element_text(family = "Arial")
    ) +
    guides(color = guide_legend(override.aes = list(linewidth = 1)))
)

# extract the first half of each for a density plot
dens_pre <- data.frame(
  x = purrr::map(
    x1,
    ~ .x$y[1:round(n * 0.75)]
  ) %>% unlist(),
  y1 = purrr::map(
    x2_grad,
    ~ .x$y[1:round(n * 0.75)]
  ) %>% unlist(),
  y2 = purrr::map(
    x2_abrupt,
    ~ .x$y[1:round(n * 0.75)]
  ) %>% unlist()
)

# extract the second half of each for a density plot
dens_post <- data.frame(
  x = purrr::map(
    x1,
    ~ .x$y[(round(n * 0.75) + 1):n]
  ) %>% unlist(),
  y1 = purrr::map(
    x2_grad,
    ~ .x$y[(round(n * 0.75) + 1):n]
  ) %>% unlist(),
  y2 = purrr::map(
    x2_abrupt,
    ~ .x$y[(round(n * 0.75) + 1):n]
  ) %>% unlist()
)




(dplot_grad_pre <- ggplot(dens_pre, aes(x, y1)) +
    geom_density2d(color = "darkgreen") +
    theme_classic() +
    theme(axis.text = element_blank()) +
    xlab(expression(x[1])) +
    ylab(expression(x[2])) +
    theme(
      axis.title.y = element_text(
        color = PNWColors::pnw_palette("Bay", 5)[1]
      )
    )
)

(dplot_grad_post <- ggplot(dens_post, aes(x, y1)) +
    geom_density2d(color = "darkorange4") +
    theme_classic() +
    theme(axis.text = element_blank()) +
    xlab(expression(x[1])) +
    ylab(expression(x[2])) +
    theme(
      axis.title.y = element_text(
        color = PNWColors::pnw_palette("Bay", 5)[1]
      )
    )
)

(dplot_abrupt_pre <- ggplot(dens_pre, aes(x, y2)) +
    geom_density2d(color = "darkgreen") +
    theme_classic() +
    theme(axis.text = element_blank()) +
    xlab(expression(x[1])) +
    ylab(expression(x[3])) +
    theme(
      axis.title.y = element_text(
        color = PNWColors::pnw_palette("Bay", 5)[5]
      )
    )
)

(dplot_abrupt_post <- ggplot(dens_post, aes(x, y2)) +
    geom_density2d(color = "darkorange4") +
    theme_classic() +
    theme(axis.text = element_blank()) +
    xlab(expression(x[1])) +
    ylab(expression(x[3])) +
    theme(
      axis.title.y = element_text(
        color = PNWColors::pnw_palette("Bay", 5)[5]
      )
    )
)

# lab1 <- cowplot::get_plot_component(
#   ggplot() +
#     labs(x = "Joint densities prior to\nforecast horizon"),
#   "xlab-b"
# )
#
# lab2 <- cowplot::get_plot_component(
#   ggplot() +
#     labs(x = "Joint densities after\nforecast horixon"),
#   "xlab-b"
# )


### begin constructing plot as a whole

lo <- "
AAAAA
BC#DE
"

(ts_plot + dplot_grad_pre + dplot_abrupt_pre +
    dplot_grad_post + dplot_abrupt_post +
    plot_layout(design = lo))

ggsave(
  here::here("Figures/covariate_shift_concept.png"),
  device = "png",
  width = 6, height = 3.5,
  units = "in"
)

