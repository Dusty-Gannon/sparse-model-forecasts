# Description: Combine the forecast + RMSE panels from the env-covariates
#              and fourier-basis tree growth analyses into one 4-panel
#              figure (a-b: env covariates; c-d: fourier basis), without
#              rerunning either script's Stan sampling.
#
#              Run Analyses/tree_dat_forecasts_env_covariates.R and
#              Analyses/tree_dat_forecasts_fourier.R first to produce the
#              two panel .rds files this script reads.

library(patchwork)
library(ggplot2)

panel_env_covs <- readRDS(
  here::here("Figures/tree_growth_forecasts_env_covs_panel.rds")
)
panel_fourier <- readRDS(
  here::here("Figures/tree_forecasts_fourier_panel.rds")
)

combined_panel <- panel_env_covs / panel_fourier

combined_panel

ggsave(
  here::here("Figures/tree_growth_forecasts_combined.pdf"),
  plot = combined_panel,
  width = 7,
  height = 7,
  device = "pdf",
  units = "in",
  dpi = 300
)
