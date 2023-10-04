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

# create a tibble to store results
  results_list <- lapply(
    file_list,
    FUN = function(f){
      x <- readRDS(f)
      list(
        conf_mat = pluck(x, "conf_summaries", "conf_mat"),
        dist_prob = pluck(x, "sim_params", "dist_prob"),
        dist_int = pluck(x, "sim_params", "dist_int"),
        prop_cdist = pluck(x, "sim_params", "prop_cdist")
      )
    }
  )

# convert to a tibble
  conf_df <- tibble(
    conf_mats = map(
      results_list,
      ~ pluck(.x, "conf_mat")
    ),
    dist_prob = map_dbl(
      results_list,
      ~ pluck(.x, "dist_prob")
    ),
    dist_int = map_dbl(
      results_list,
      ~ pluck(.x, "dist_int")
    ),
    prop_cdist = map_dbl(
      results_list,
      ~ pluck(.x, "prop_cdist")
    )
  )

# create a summarized version
  conf_df_sum <- conf_df %>% group_by(
    dist_int, prop_cdist, dist_prob
  ) %>% summarise(
    conf_mat_sum = list(Reduce("+", conf_mats))
  )

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
    conf_df_sum[rep(1:nrow(conf_df_sum), each = 4), 1:3]
  )

  df_subs <- vector(mode = "list", length = 2)
  df_subs[[1]] <- rbind(
    metrics_df[1:4, ],
    dplyr::filter(metrics_df, dist_int == 0.5)
  )
  df_subs[[2]] <- rbind(
    metrics_df[1:4, ],
    dplyr::filter(metrics_df, dist_int == 0.8)
  )


# add rows for the missing combos to make the heatmap square
  empties <- data.frame(
    dist_prob = rep(c(rep(0, 10), unique(conf_df_sum$dist_prob)[-1]), each = 4),
    prop_cdist = rep(c(unique(conf_df_sum$prop_cdist)[-1], rep(0, 10)), each = 4),
    dist_int = NA,
    metric = rep(unique(metrics_df$metric), 20),
    value = NA,
    low = NA,
    high = NA
  )

# insert into the other dfs
  for(i in 1:length(df_subs)){
    df_subs[[i]] <- rbind(
      df_subs[[i]][1:4, ],
      empties,
      df_subs[[i]][5:nrow(df_subs[[i]]), ]
    )
  }


  # create the plots
  plot_hm <- function(df, title){
    library(ggplot2)
    ggplot(data = df, aes(x = dist_prob, y = prop_cdist)) +
      geom_tile(aes(fill = value)) +
      theme(
        panel.background = element_blank(),
        axis.line = element_line(color = "darkgrey", size = 0.1),
        legend.title = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)
      ) +
      scale_fill_gradient(low = "grey", high = "brown", na.value = "white") +
      scale_x_continuous(labels = scales::label_number(accuracy = 0.001)) +
      ggtitle(title) +
      xlab("") +
      ylab("")
  }

  # create nested datasets
  dfs_nest <- purrr::map(
    df_subs,
    ~ {.x %>%
    group_by(metric) %>%
    nest()}
  )

  # create plots
  plots1 <- map2(
    dfs_nest[[1]]$metric,
    dfs_nest[[1]]$data,
    ~ plot_hm(.y, .x)
  )

  plots2 <- purrr::map2(
    dfs_nest[[2]]$metric,
    dfs_nest[[2]]$data,
    ~ plot_hm(.y, .x)
  )

  xlabel <- ggplot(data.frame(l = "Disturbance probability", x = 1, y = 1)) +
    geom_text(aes(x, y, label = l), size = 5) +
    theme_void()
  ylabel <- ggplot(data.frame(l = "Proportion of community disturbed", x = 1, y = 1)) +
    geom_text(aes(x, y, label = l), angle = 90, size = 5) +
    theme_void()
  p1label <- ggplot(data.frame(l = "Disturbance intensity = 0.5", x = 1, y = 1)) +
    geom_text(aes(x, y, label = l), size = 5) +
    theme_void()
  p2label <- ggplot(data.frame(l = "Disturbance intensity = 0.8", x = 1, y = 1)) +
    geom_text(aes(x, y, label = l), size = 5) +
    theme_void()

  # layout
  lo <- '
    #AA
    ###
    #BC
    DEF
    D##
    DGG
    D##
    DHJ
    #KL
    #MM
  '

  p <- p1label + plots1[[1]] + plots1[[2]] + ylabel + plots1[[3]] + plots1[[4]] +
    p2label + plots2[[1]] + plots2[[2]] + plots2[[3]] + plots2[[4]] + xlabel +
   plot_layout(
     design = lo,
     widths = unit(c(0.75, 4, 4), "cm"),
     heights = unit(c(0.5, 0.2, 4, 4, 0.2, 0.5, 0.2, 4, 4, 0.5), "cm")
   )

  # save the plots
  ggsave(
    p,
    filename = here(
      paste0("Figures/confusion_metrics_dist_freq_x_int_", args[2], ".png")
    ),
    height = 14, width = 7,
    units = "in"
  )







