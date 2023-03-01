##########################################################
# This script takes output from fitting the Ricker models,
# with different ratios of focal demographic stochasticity
# relative to the variance of heterospecific abundances
##########################################################

library(tidyverse)
library(here)
library(gridExtra)

# load data sims and results from fits
  fits <- readRDS(here("Data/terrestrial_sim_data/lnorm_ricker/vr_tests.rds"))
  sims <- readRDS(here("Data/terrestrial_sim_data/lnorm_ricker/lnorm_ricker_sims_200steps_rho0_S5_s55.rds"))
  mult <- c(0.5, 1, 2, 4, 10)

# sum over the confusion matrices for each category
  # initialize
  conf_mats_sum <- purrr::map(
    1:length(fits[[1]]),
    ~ matrix(data = 0, nrow = 2, ncol = 2)
  )
  # add to it
  for(i in 1:length(fits)){
    for(j in 1:length(fits[[i]])){
      conf_mats_sum[[j]] <- conf_mats_sum[[j]] +
        fits[[i]][[j]]$conf_mat
    }
  }

  # lumping the simulations into categories for the multiplier
  conf_df_lumped <- data.frame(
    multiplier = c("full", mult),
    TPR = purrr::map_dbl(
      conf_mats_sum,
      ~ {.x[1, 1] / sum(.x[1, ])}
    ),
    TPR_tot = purrr::map_dbl(
      conf_mats_sum,
      ~ {sum(.x[1, ])}
    ),
    TNR = purrr::map_dbl(
      conf_mats_sum,
      ~ {.x[2, 2] / sum(.x[2, ])}
    ),
    TNR_tot = purrr::map_dbl(
      conf_mats_sum,
      ~ {sum(.x[2, ])}
    ),
    Precision = purrr::map_dbl(
      conf_mats_sum,
      ~ {.x[1, 1] / sum(.x[, 1])}
    ),
    prec_tot = purrr::map_dbl(
      conf_mats_sum,
      ~ {sum(.x[, 1])}
    ),
    Accuracy = purrr::map_dbl(
      conf_mats_sum,
      ~ {sum(diag(.x)) / sum(.x)}
    ),
    acc_tot = purrr::map_dbl(
      conf_mats_sum,
      ~ {sum(.x)}
    )
  )

# create dataframe with unique vrs from simulations and metrics of interest
  N_vrtests_list <- purrr::map(
    sims,
    ~ {.x$N_vrtests}
  ) %>% flatten()

  confmats_manip <- purrr::map(
    fits,
    ~ {.x[2:6]}
  ) %>% flatten() %>% purrr::map(., ~ {.x$conf_mat})

# create vector of focal scl of dem. stoch.
  sigmas <- rep(
    purrr::map_dbl(
      sims,
      ~ {.x$sim_params$sigmas[1]}
    ),
    each = 5
  ) / rep(sqrt(mult), 200)


# fill in a dataframe with unique vrs
  conf_df_unique <- data.frame(
    ratio = purrr::map2_dbl(
      N_vrtests_list,
      sigmas,
      ~ {mean(apply(.x, 1, sd)) / .y}
    ),
    TPR = purrr::map_dbl(
      confmats_manip,
      ~ {.x[1, 1] / sum(.x[1, ])}
    ),
    TNR = purrr::map_dbl(
      confmats_manip,
      ~ {.x[2, 2] / sum(.x[2, ])}
    ),
    Precision = purrr::map_dbl(
      confmats_manip,
      ~ {.x[1, 1] / sum(.x[, 1])}
    ),
    Accuracy = purrr::map_dbl(
      confmats_manip,
      ~ {sum(diag(.x)) / sum(.x)}
    )
  )

# pivot longer for ggplot
  sumdf_long <- cbind(
    conf_df_lumped %>% select(
      !contains("tot")
    ) %>% pivot_longer(
      .,
      cols = TPR:Accuracy,
      names_to = "metric",
      values_to = "rate"
    ),
    conf_df_lumped %>% select(
      contains("tot")
    ) %>% pivot_longer(
      .,
      cols = everything(),
      names_to = "type",
      values_to = "total"
    )
  ) %>% select(., !type)

# add columns for confidence intervals
  sumdf_long <- sumdf_long %>% mutate(
    low = purrr::map2_dbl(
      rate, total,
      ~ binom.test(.x * .y, n = .y)$conf.int[1]
    ),
    high = purrr::map2_dbl(
      rate, total,
      ~ binom.test(.x * .y, n = .y)$conf.int[2]
    )
  )

# split for full and maniptulated data
  sumdf_manip <- filter(sumdf_long, multiplier != "full")
  sumdf_fullcomm <- filter(sumdf_long, multiplier == "full") %>%
    select(., !multiplier)

# long data for unique estimates
  uniqdf_long <- conf_df_unique %>% pivot_longer(
    cols = TPR:Accuracy,
    names_to = "metric",
    values_to = "rate"
  )

# plot and save
  p1 <- ggplot(data = sumdf_manip, aes(x = as.numeric(multiplier), y = rate)) +
    facet_grid(rows = vars(metric)) +
    geom_line(color = "brown") +
    geom_errorbar(aes(ymin = low, ymax = high), color = "brown", width = 0.1) +
    geom_point(size = 2) +
    geom_point(color = "brown") +
    theme_bw() +
    xlab("multiplier")

  p2 <- ggplot(data = uniqdf_long, aes(x = ratio, y = rate)) +
    facet_grid(rows = vars(metric)) +
    geom_point(color = "brown", alpha = 0.3) +
    geom_smooth(se = F, color = "black", span = 20) +
    theme_bw()

  png(
    here("Figures/Ricker_confusion_metrics_vr_tests.png"),
    width = 1500, height = 1800, units = "px",
    res = 300
  )
  grid.arrange(
    arrangeGrob(p1, top = "Multiplier"),
    arrangeGrob(p2, top = "Exact ratio"),
    nrow = 1
  )
  dev.off()









