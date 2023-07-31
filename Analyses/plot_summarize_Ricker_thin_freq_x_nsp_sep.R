##########################################################
# This script takes output from fitting the Ricker models,
# with different thinning intervals and proportion of the
# community thinned, where results from model fits are
# stored in separate .rds files, and plots confusion
# metrics in a heatmap.
##########################################################

library(tidyverse)
library(here)
library(patchwork)
devtools::load_all()

# load data sims and results from fits
# arguments are:
#  1. relative filepath to directory with model fit results
#  2. relative filepath for figure output
#  3. optional search pattern for listing files
  args <- commandArgs(trailingOnly = T)
  args

# create a vector of the file names
  if(length(args) == 3){
    file_list <- list.files(
      here(args[1]),
      full.names = T,
      pattern = args[3]
    )
  } else {
    file_list <- list.files(
      here(args[1]),
      full.names = T
    )
  }

# load results into a list
  results_list <- lapply(
    file_list,
    FUN = function(f){
      x <- readRDS(f)
      list(
        conf_mat = pluck(x, "conf_summaries", "conf_mat"),
        thin_freq = pluck(x, "sim_params", "thin_freq"),
        prop_cthin = pluck(x, "sim_params", "prop_cthin")
      )
    }
  )

# convert to a tibble
  conf_df <- tibble(
    conf_mats = map(
      results_list,
      ~ pluck(.x, "conf_mat")
    ),
    thin_interval = map_dbl(
      results_list,
      ~ pluck(.x, "thin_freq")
    ),
    prop_community = map_dbl(
      results_list,
      ~ pluck(.x, "prop_cthin")
    )
  )

# create a summarized version
  conf_df_sum <- conf_df %>% group_by(
    thin_interval, prop_community
  ) %>% summarise(
    conf_mat_sum = list(Reduce("+", conf_mats))
  ) %>% ungroup()

# function to calculate confusion metrics
  confusion_metrics <- function(mat){

    return(
      data.frame(
        metric = c("TPR", "TNR", "Precision", "Accuracy"),
        value = c(
          if(sum(mat[1, ]) > 0){mat[1, 1] / sum(mat[1, ])} else {NA},
          if(sum(mat[2, ]) > 0){mat[2, 2] / sum(mat[2, ])} else {NA},
          if(sum(mat[, 1]) > 0){mat[1, 1] / sum(mat[, 1])} else {NA},
          sum(diag(mat)) / sum(mat)
        ),
        low = c(
          if(sum(mat[1, ]) > 0){binom.test(mat[1, 1], n = sum(mat[1, ]))$conf.int[1]} else {NA},
          if(sum(mat[2, ]) > 0){binom.test(mat[2, 2], n = sum(mat[2, ]))$conf.int[1]} else {NA},
          if(sum(mat[, 1]) > 0){binom.test(mat[1, 1], n = sum(mat[, 1]))$conf.int[1]} else {NA},
          binom.test(sum(diag(mat)), n = sum(mat))$conf.int[1]
        ),
        high = c(
          if(sum(mat[1, ]) > 0){binom.test(mat[1, 1], n = sum(mat[1, ]))$conf.int[2]} else {NA},
          if(sum(mat[1, ]) > 0){binom.test(mat[2, 2], n = sum(mat[2, ]))$conf.int[2]} else {NA},
          if(sum(mat[, 1]) > 0){binom.test(mat[1, 1], n = sum(mat[, 1]))$conf.int[2]} else {NA},
          binom.test(sum(diag(mat)), n = sum(mat))$conf.int[2]
        )
      )
    )

  }

# calculate metrics for each row and convert to long format
  metrics_list <- lapply(conf_df_sum$conf_mat_sum, confusion_metrics)
  metrics_df <- Reduce(rbind, metrics_list)

# add columns for the input parameters
  metrics_df <- cbind(
    metrics_df,
    select(
      conf_df_sum,
      thin_interval:prop_community
    )
  )


# add rows for the missing combos to make the heatmap square
  empties <- data.frame(
    metric = rep(unique(metrics_df$metric), 20),
    value = NA,
    low = NA,
    high = NA,
    thin_interval = rep(c(rep(0, 10), unique(conf_df_sum$thin_interval)[-1]), each = 4),
    prop_community = rep(c(unique(conf_df_sum$prop_community)[-1], rep(0, 10)), each = 4)
  )
  metrics_df <- rbind(
    empties,
    metrics_df
  )

  # create the plots
  plot_hm <- function(df, title){
    library(ggplot2)
    ggplot(data = df, aes(x = thin_interval, y = prop_community)) +
      geom_tile(aes(fill = value)) +
      theme(
        panel.background = element_blank(),
        axis.line = element_line(color = "darkgrey", size = 0.1),
        legend.title = element_blank()
      ) +
      scale_fill_gradient(low = "grey", high = "brown", na.value = "white") +
      ggtitle(title) +
      xlab("") +
      ylab("")
  }

  dat_nest <- metrics_df %>%
    group_by(metric) %>%
    nest()

  plots <- map2(
    dat_nest$metric,
    dat_nest$data,
    ~ plot_hm(.y, .x)
  )

  xlabel <- ggplot(data.frame(l = "Thinning interval", x = 1, y = 1)) +
    geom_text(aes(x, y, label = l), size = 5) +
    theme_void()
  ylabel <- ggplot(data.frame(l = "Proportion of community thinned", x = 1, y = 1)) +
    geom_text(aes(x, y, label = l), angle = 90, size = 5) +
    theme_void()

  # layout
  lo <- '
    ABC
    ADE
    #FF
  '

  # save the plots
  ggsave(
    ylabel + plots[[1]] + plots[[2]] + plots[[3]] + plots[[4]] + xlabel +
      plot_layout(design = lo, widths = unit(c(0.75, 4, 4), "cm"), heights = unit(c(4, 4, 0.75), "cm")),
    filename = here(
      paste0("Figures/confusion_metrics_freq_x_nsp_thin_", args[2], ".png")
    ),
    height = 7, width = 7,
    units = "in"
  )







