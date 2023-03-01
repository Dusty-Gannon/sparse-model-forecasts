##########################################################
# This script takes output from fitting the Ricker models,
# with different thinning intervals.
##########################################################

library(tidyverse)
library(here)
library(gridExtra)
library(scales)

# load data sims and results from fits
  fits_ordered <- readRDS(here("Data/terrestrial_sim_data/lnorm_ricker/ordered_thin_tests.rds"))
  sims_ordered <- readRDS(here("Data/terrestrial_sim_data/lnorm_ricker/lnorm_ricker_thin_sims_ordered_S5_s55.rds"))

# pull out confusion matrices
  conf_mats <- purrr::map(
    fits_ordered,
    ~ .x$conf_summaries$conf_mat
  )

# confusion matrices by group for the manipulated communities
  thin_freq <- purrr::map_dbl(
    sims_ordered,
    ~ .x$sim_params$thin_freq
  )
  grps <- sort(unique(thin_freq))
  conf_mats_sum <-  vector(mode = "list", length = length(grps))

  for(j in 1:length(grps)){
    sub_ids <- which(thin_freq == grps[j])
    conf_mats_sum[[j]] <- Reduce("+", conf_mats[sub_ids])
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

  metrics_df <- data.frame(NULL)
  for(j in 1:length(grps)){
    metrics_df <- rbind(
      metrics_df, confusion_metrics(conf_mats_sum[[j]])
    )
  }

  metrics_df <- metrics_df %>% mutate(
    thin_freq = grps
  )

  # convert to long format for plotting
  conf_metrics_long <- cbind(
    metrics_df %>% select(
      TPR, TNR, Precision, Accuracy, thin_freq
    ) %>% pivot_longer(
      cols = TPR:Accuracy,
      values_to = "rate",
      names_to = "metric"
    ),
    metrics_df %>% select(
      TPR_low, TNR_low, prec_low, prec_high
    ) %>% pivot_longer(
      cols = everything(),
      values_to = "low",
      names_to = "name"
    ),
    metrics_df %>% select(
      TPR_high, TNR_high, prec_high, acc_high
    ) %>% pivot_longer(
      cols = everything(),
      values_to = "high",
      names_to = "name2"
    )
  )

  # remove the names_to columns
  conf_metrics_long <- conf_metrics_long %>%
    select(!starts_with("name"))

  # create the plots
  confmet_plot <- ggplot(conf_metrics_long, aes(x = thin_freq, y = rate)) +
    facet_grid(rows = vars(metric)) +
    geom_line(color = "brown") +
    geom_point(color = "brown") +
    geom_errorbar(aes(ymin = low, ymax = high), width = 0.1, color = "brown") +
    theme_bw() +
    # theme(
    #   legend.title = element_blank(),
    #   panel.grid = element_blank()
    # ) +
    # scale_y_continuous(labels = label_number(accuracy = 0.1)) +
    xlab("thinning frequency") +
    ylab("rate")


  # save the plot
  ggsave(
    filename = here("Figures/confusion_metrics_ordered_thin.png"),
    width = 3, height = 5,
    units = "in"
  )







