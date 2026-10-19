#

library(ggplot2)
library(dplyr)
library(tidyr)
library(here)
library(patchwork)

# ---- load data ----
arp_err_sims <- readr::read_csv(
    here("Data/fourier_sims/ARp_err_sims_20260928_condensed.csv")
)

# ---- RMSE plot ----

colors <- PNWColors::pnw_palette("Cascades", n = 4)
fills <- adjustcolor(colors, alpha.f = 0.6)

labs <- c(
    "0" = "0 environmental features",
    "2" = "2 environmental features",
    "4" = "4 environmental features"
)

p <- ggplot(arp_err_sims, aes(x = rmse_forecast, y = model)) +
    facet_wrap(vars(betas), ncol = 1, labeller = as_labeller(labs)) +
    geom_violin(aes(fill = model, color = model), adjust = 2) +
    theme_bw() +
    stat_summary(
        color = "grey3",
        geom = "errorbar",
        fun.min = \(x){quantile(x, probs = 0.025)},
        fun.max = \(x){quantile(x, probs = 0.975)},
        width = 0,
        linewidth = 0.5
    ) +
    stat_summary(
        color = "grey3",
        geom = "errorbar",
        fun = mean,
        fun.min = \(x){quantile(x, probs = 0.1)},
        fun.max = \(x){quantile(x, probs = 0.9)},
        linewidth = 1,
        width = 0
    ) +
    stat_summary(
        geom = "point",
        aes(color = model),
        fun = mean,
        size = 1.5
    ) +
    stat_summary(
        geom = "point",
        fun = mean,
        color = "white",
        size = 0.7
    ) +
    geom_vline(xintercept = 1, linetype = "dashed") +
    scale_fill_manual(values = fills[c(2,4)]) +
    scale_color_manual(values = colors[c(2,4)]) +
    scale_y_discrete(labels = c("stepAIC +\nauto-ARIMA", "RHS")) +
    theme(
        legend.position = "none"
    ) +
    xlab("Forecast RMSE") +
    ylab("")

ggsave(
    p,
    file = here("Figures/pred_rmse_fourier.pdf"),
    device = "pdf",
    width = 4,
    height = 6,
    units = "in",
    dpi = 300
)

# ---- Pred Interval coverage plot ----

pi_cov <- ggplot(arp_err_sims, aes(x = PIcoverage, y = model)) +
    facet_wrap(vars(betas), ncol = 1, labeller = as_labeller(labs)) +
    stat_summary(
        geom = "pointrange",
        aes(color = model),
        fun = mean,
        fun.min = \(x) { quantile(x, probs = 0.025) },
        fun.max = \(x) { quantile(x, probs = 0.975) }
    ) +
    theme_bw() +
    geom_vline(xintercept = 0.95, linetype = "dashed") +
    scale_fill_manual(values = fills[c(2,4)]) +
    scale_color_manual(values = colors[c(2,4)]) +
    scale_y_discrete(labels = c("stepAIC +\nauto-ARIMA", "RHS")) +
    theme(
        legend.position = "none",
        plot.title = element_text(hjust = -0.1)
    ) +
    ylab("") +
    xlab("Coverage rate") +
    ggtitle("a) Prediction interval coverage")

pi_width <- ggplot(arp_err_sims, aes(x = PIwidth, y = model)) +
    facet_wrap(vars(betas), ncol = 1, labeller = as_labeller(labs)) +
    stat_summary(
        geom = "pointrange",
        aes(color = model),
        fun = mean,
        fun.min = \(x) { quantile(x, probs = 0.025) },
        fun.max = \(x) { quantile(x, probs = 0.975) }
    ) +
    theme_bw() +
    scale_fill_manual(values = fills[c(2,4)]) +
    scale_color_manual(values = colors[c(2,4)]) +
    scale_y_discrete(labels = c("stepAIC +\nauto-ARIMA", "RHS")) +
    theme(
        legend.position = "none",
        plot.title = element_text(hjust = -0.1)
    ) +
    ylab("") +
    xlab("Width") +
    ggtitle("b) Prediction interval width")

pi_plot <- pi_cov + pi_width + plot_layout(axes = 'collect_y')

ggsave(
    pi_plot,
    file = here("Figures/pred_int_res_fourier.pdf"),
    width = 7, height = 7,
    units = "in",
    device = "pdf",
    dpi = 300
)

# ---- In-text figures ----

arp_err_sims %>% 
    group_by(model, betas) %>%
    summarise(
        n = n(),
        mean = mean(rmse_forecast),
        min = min(rmse_forecast),
        q95 = quantile(rmse_forecast, probs = 0.95),
        pi_cov = mean(PIcoverage),
        width = mean(PIwidth)
    )


