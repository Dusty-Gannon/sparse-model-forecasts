################################################################
# This script is used to create a conceptual figure illustrating
# covariate shift with a simple, 2 covariate example.
#
# AUTHOR: Dusty Gannon
# LAST UPDATED: 12 Sept 2025
################################################################

# Load required libraries
devtools::load_all()
library(tidyverse)
library(patchwork)
set.seed(3778)

# Set series length
n <- 120
t <- 1:n

################################################################
# Define functions to construct different types of time series
################################################################

# Function to construct a baseline main time series
construct_main <- function(m, h, tmin, tmax, n, sd) {
  t <- seq(tmin, tmax, length.out = n)
  y <- m * t / (h + t) + arima.sim(list(ar = 0.6), n = n, sd = sd)
  return(data.frame(t = t, y = y))
}

# Function to construct a time series with a gradual change
construct_grad_2 <- function(t0, m, h, tmin, tmax, n, sd) {
  t <- seq(tmin, tmax, length.out = n)
  a <- m * t0 / (h + t0)
  b <- h * m / (h + t0)^2
  y <- a + b * (t - t0) + arima.sim(list(ar = 0.6), n = n, sd = sd)
  return(data.frame(t = t, y = y))
}

# Function to construct a time series with an abrupt change
construct_abrupt_2 <- function(t0, m, h, tmin, tmax, n, sd) {
  t <- seq(tmin, tmax, length.out = n)
  a <- m * t0 / (h + t0)
  b <- h * m / (h + t0)^2
  b2 <- -h * m / (h + t0)^3
  y <- a + b * (t - t0) + b2 * (t - t0)^2 + arima.sim(list(ar = 0.6), n = n, sd = sd)
  return(data.frame(t = t, y = y))
}

# Set model parameters
m <- 3
h <- 0.4 * n |> round()
tmin <- 1
tmax <- n
t0 <- 35

################################################################
# Simulate 1000 realizations from each time series model
################################################################

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

################################################################
# Create example time series data for plotting
################################################################

ts_data <- data.frame(
  x_var = rep(c("x1", "x2", "x3"), each = n),
  x = c(x1[[1]]$y, x2_grad[[1]]$y, x2_abrupt[[1]]$y),
  time = 1:n %>% rep(3)
)

# Plot the time series with plot annotations
(ts_plot <- ggplot(ts_data, aes(time, x, color = x_var)) +
    annotate(
      geom = "rect",
      xmin = 0, xmax = round(n * 0.75), ymin = min(ts_data$x), ymax = 3.5,
      fill = "darkgreen", alpha = 0.2, color = NA
    ) +
    annotate(
      geom = "rect",
      xmin = round(n * 0.75) + 1, xmax = n, ymin = min(ts_data$x), ymax = 3.5,
      fill = "darkorange4", color = NA, alpha = 0.2
    ) +
    geom_line() +
    theme_classic() +
    scale_color_manual(
      values = c("black", PNWColors::pnw_palette("Bay", 5)[c(1, 5)]),
      labels = expression(x[1], x[2], x[3])
    ) +
    theme(
      text = element_text(family = "Arial"),
      legend.title = element_text(size = 10)
    ) +
    guides(color = guide_legend(override.aes = list(linewidth = 1))) +
    labs(color = "Covariate")
)

################################################################
# Extract pre- and post-forecast data for density plots
################################################################

# Pre-forecast horizon
dens_pre <- data.frame(
  x = purrr::map(x1, ~ .x$y[1:round(n * 0.75)]) %>%
    unlist(),
  y1 = purrr::map(x2_grad, ~ .x$y[1:round(n * 0.75)]) %>%
    unlist(),
  y2 = purrr::map(x2_abrupt, ~ .x$y[1:round(n * 0.75)]) %>%
    unlist()
)

# Post-forecast horizon
dens_post <- data.frame(
  x = purrr::map(x1, ~ .x$y[(round(n * 0.75) + 1):n]) %>%
    unlist(),
  y1 = purrr::map(x2_grad, ~ .x$y[(round(n * 0.75) + 1):n]) %>%
    unlist(),
  y2 = purrr::map(x2_abrupt, ~ .x$y[(round(n * 0.75) + 1):n]) %>%
    unlist()
)

################################################################
# Create density plots for each condition
################################################################

(dplot_grad_pre <- ggplot(dens_pre, aes(x, y1)) +
   geom_density2d(color = "darkgreen", bins = 7) +
   theme_classic() +
   theme(axis.text = element_blank()) +
   xlab(expression(x[1])) +
   ylab(expression(x[2])) +
   theme(
     axis.title.y = element_text(color = PNWColors::pnw_palette("Bay", 5)[1])
   )
)

(dplot_grad_post <- ggplot(dens_post, aes(x, y1)) +
    geom_density2d(color = "darkorange4", bins = 7) +
    theme_classic() +
    theme(axis.text = element_blank()) +
    xlab(expression(x[1])) +
    ylab(expression(x[2])) +
    theme(
      axis.title.y = element_text(color = PNWColors::pnw_palette("Bay", 5)[1])
    )
)

(dplot_abrupt_pre <- ggplot(dens_pre, aes(x, y2)) +
    geom_density2d(color = "darkgreen", bins = 7) +
    theme_classic() +
    theme(axis.text = element_blank()) +
    xlab(expression(x[1])) +
    ylab(expression(x[3])) +
    theme(
      axis.title.y = element_text(color = PNWColors::pnw_palette("Bay", 5)[5])
    )
)

(dplot_abrupt_post <- ggplot(dens_post, aes(x, y2)) +
    geom_density2d(color = "darkorange4", bins = 7) +
    theme_classic() +
    theme(axis.text = element_blank()) +
    xlab(expression(x[1])) +
    ylab(expression(x[3])) +
    theme(
      axis.title.y = element_text(color = PNWColors::pnw_palette("Bay", 5)[5])
    )
)

################################################################
# Combine plots and save the figure
################################################################

l1 <- grid::textGrob(
  label = "Joint densities during\ntraining period",
  x = unit(0.55, "npc"), y = unit(0.9, "npc"),
  gp = grid::gpar(fontsize = 10),
  vjust = 1
)

l2 <- grid::textGrob(
  label = "Joint densities during\ntesting period",
  x = unit(0.4, "npc"), y = unit(0.9, "npc"),
  gp = grid::gpar(fontsize = 10, fontfamily = "Arial"),
  vjust = 1
)

# Define plot layout
lo <- "
AAAAAAAAA
AAAAAAAAA
AAAAAAAAA
BBCC#DDEE
BBCC#DDEE
FFFF#GGGG
"

(ts_plot + dplot_grad_pre + dplot_abrupt_pre +
    dplot_grad_post + dplot_abrupt_post + l1 + l2 +
    plot_layout(design = lo))

# Save the plot as a PNG file
ggsave(
  here::here("Figures/covariate_shift_concept.jpg"),
  device = "jpeg",
  width = 6, height = 4.5,
  units = "in", dpi = "retina"
)
