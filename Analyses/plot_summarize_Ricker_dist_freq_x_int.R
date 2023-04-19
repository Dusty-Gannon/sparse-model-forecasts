##########################################################
# This script takes output from fitting the Ricker models,
# with different thinning intervals and proportion of the
# community thinned, and plots confusion metrics in a
# heatmap.
##########################################################

library(tidyverse)
library(here)
library(patchwork)
devtools::load_all()

# load data sims and results from fits
  args <- commandArgs(trailingOnly = T)
  args
  sims <- readRDS(here(
    paste0("Data/terrestrial_sim_data/lnorm_ricker/", args[1])
  ))

  file_list <- list.files(
    here(
      paste0("Data/terrestrial_sim_data/lnorm_ricker/disturb_results/", args[2])
    ),
    full.names = T,
    pattern = "disturb_tests"
  )
  # put the file list in order
  sort_pat <- as.numeric(
    stringr::str_extract(file_list, "tests_[:digit:]+") %>%
      stringr::str_extract(., "[:digit:]+")
  )
  file_list <- file_list[order(sort_pat)]

  # loop through and add the components we want to a full list
  conf_mats <- vector(mode = "list")
  # rmse <- vector(mode = "list")
  for(i in 1:length(file_list)){
    l_i <- readRDS(file_list[i])
    conf_mats_i <- purrr::map(
      l_i,
      ~ .x$conf_summaries$conf_mat
    )
    conf_mats <- c(conf_mats, conf_mats_i)
    # rmse_i <- purrr::map(
    #   l_i,
    #   ~ .x$conf_summaries$rmse
    # )
    # rmse <- c(rmse, rmse_i)
  }

# create dataframe with crossed treatments
  conf_df_full <- data.frame(
    dist_prob = purrr::map_dbl(
      sims,
      ~ .x$sim_params$dist_prob
    ),
    prop_cdist = purrr::map_dbl(
      sims,
      ~ .x$sim_params$prop_cdist
    ),
    dist_int = purrr::map_dbl(
      sims,
      ~ .x$sim_params$dist_int
    )
  )

# create a summarized version
  conf_df_sum <- unique(conf_df_full)

# confusion matrices by group
  grps <- nrow(conf_df_sum)
  reps <- length(sims) / grps

  conf_mats_sum <-  vector(mode = "list", length = grps)
  for(j in 1:grps){
    sub_ids <- (reps * (j - 1) + 1):(reps * j)
    conf_mats_sum[[j]] <- Reduce("+", conf_mats[sub_ids])
  }

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

  metrics_df <- data.frame(NULL)
  for(j in 1:grps){
    metrics_df <- rbind(
      metrics_df, confusion_metrics(conf_mats_sum[[j]])
    )
  }

  conf_df_sum <- cbind(
    conf_df_sum[rep(1:grps, each = 4), ],
    metrics_df
  )

  df_subs <- vector(mode = "list", length = 2)
  df_subs[[1]] <- rbind(
    conf_df_sum[1:4, ],
    dplyr::filter(conf_df_sum, dist_int == 0.5)
  )
  df_subs[[2]] <- rbind(
    conf_df_sum[1:4, ],
    dplyr::filter(conf_df_sum, dist_int == 0.8)
  )


# add rows for the missing combos to make the heatmap square
  empties <- data.frame(
    dist_prob = rep(c(rep(0, 10), unique(conf_df_sum$dist_prob)[-1]), each = 4),
    prop_cdist = rep(c(unique(conf_df_sum$prop_cdist)[-1], rep(0, 10)), each = 4),
    dist_int = NA,
    metric = rep(unique(conf_df_sum$metric), 20),
    value = NA,
    low = NA,
    high = NA
  )

# combine with the other dfs
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
        axis.text.x = element_text(angle = 45, hjust = 0)
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
     heights = unit(c(0.5, 0.1, 4, 4, 0.1, 0.5, 0.1, 4, 4, 0.75), "cm")
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







