##########################################################
# This script takes output from fitting the Ricker models,
# with different ratios of focal demographic stochasticity
# relative to the variance of heterospecific abundances
##########################################################

library(tidyverse)
library(here)
library(gridExtra)
library(scales)

# load data sims and results from fits
  fits <- readRDS(here("Data/terrestrial_sim_data/lnorm_ricker/vr_tests_1param.rds"))

# constant to change the ratios
  het_vr <- c(0.5, 1, 2, 5, 20)
  grps <- length(het_vr)

  # pull out confusion matrices
  conf_mats <- list(
    full = purrr::map(fits, ~ .x$N_full$conf_mat),
    cor = purrr::map(fits, ~ .x$N_cor$conf_mat)
  )

  # summarize the full community confusion matrices
  confmat_full <- Reduce("+", conf_mats[[1]])

  # confusion matrices by group for the manipulated communities
  conf_mats_sum <-  vector(mode = "list", length = grps)

    for(j in 1:grps){
      sub_ids <- ((j - 1) * 100 + 1):(j * 100)
      conf_mats_sum[[j]] <- Reduce("+", conf_mats[[2]][sub_ids])
    }

  confusion_metrics <- function(mat){

    return(
      data.frame(
        TPR = mat[1, 1] / sum(mat[1, ]),
        TPR_low = binom.test(mat[1, 1], n = sum(mat[1, ]))$conf.int[1],
        TPR_high = binom.test(mat[1, 1], n = sum(mat[1, ]))$conf.int[2],
        TNR = mat[2, 2] / sum(mat[2, ]),
        TNR_low = binom.test(mat[2, 2], n = sum(mat[2, ]))$conf.int[1],
        TNR_high = binom.test(mat[2, 2], n = sum(mat[2, ]))$conf.int[2],
        Precision = mat[1, 1] / sum(mat[, 1]),
        prec_low = binom.test(mat[1, 1], n = sum(mat[, 1]))$conf.int[1],
        prec_high = binom.test(mat[1, 1], n = sum(mat[, 1]))$conf.int[2],
        Accuracy = sum(diag(mat)) / sum(mat),
        acc_low = binom.test(sum(diag(mat)), n = sum(mat))$conf.int[1],
        acc_high = binom.test(sum(diag(mat)), n = sum(mat))$conf.int[2]
      )
    )

  }

  metrics_df_full <- confusion_metrics(confmat_full)

  metrics_df_manip <- data.frame(NULL)
  for(j in 1:grps){
    metrics_df_manip <- rbind(
      metrics_df_manip, confusion_metrics(conf_mats_sum[[j]])
    )
  }

  metrics_df_manip <- metrics_df_manip %>% mutate(
    constant = het_vr
  )

  # convert to long format for plotting
  conf_metrics_long <- cbind(
    metrics_df_manip %>% select(
      TPR, Precision, Accuracy, constant
    ) %>% pivot_longer(
      cols = TPR:Accuracy,
      values_to = "rate",
      names_to = "metric"
    ),
    metrics_df_manip %>% select(
      TPR_low, prec_low, prec_high
    ) %>% pivot_longer(
      cols = everything(),
      values_to = "low",
      names_to = "name"
    ),
    metrics_df_manip %>% select(
      TPR_high, prec_high, acc_high
    ) %>% pivot_longer(
      cols = everything(),
      values_to = "high",
      names_to = "name2"
    )
  )

  # remove the names_to columns
  conf_metrics_long <- conf_metrics_long %>%
    select(!starts_with("name"))

  # create long_format df for the full sims
  full_df_long <- cbind(
    metrics_df_full %>% select(
      TPR, Precision, Accuracy
    ) %>% pivot_longer(
      everything(),
      values_to = "rate",
      names_to = "metric"
    ),

    metrics_df_full %>% select(
      TPR_low, prec_low, acc_low
    ) %>% pivot_longer(
      everything(),
      values_to = "low",
      names_to = "name1"
    ),

    metrics_df_full %>% select(
      TPR_high, prec_high, acc_high
    ) %>% pivot_longer(
      everything(),
      values_to = "high",
      names_to = "name2"
    )
  )

  # remove the names_to columns
  full_df_long <- full_df_long %>%
    select(!starts_with("name"))


  # create the plots
  confmet_plot <- ggplot(conf_metrics_long, aes(x = constant, y = rate)) +
    facet_grid(rows = vars(metric)) +
    geom_hline(data = full_df_long, aes(yintercept = rate), col = "grey") +
    geom_hline(data = full_df_long, aes(yintercept = low), color = "grey", linetype = "dashed") +
    geom_hline(data = full_df_long, aes(yintercept = high), color = "grey", linetype = "dashed") +
    geom_line(color = "brown") +
    geom_point(color = "brown") +
    geom_errorbar(aes(ymin = low, ymax = high), width = 0.1, color = "brown") +
    theme_bw() +
    theme(
      legend.title = element_blank(),
      panel.grid = element_blank()
    ) +
    scale_y_continuous(labels = label_number(accuracy = 0.1)) +
    xlab("amplification factor") +
    ylab("rate")


  # save the plot
  ggsave(
    filename = here("Figures/confusion_metrics_1param.png"),
    width = 3, height = 5,
    units = "in"
  )







