##########################################################
# This script takes output from fitting the Ricker models,
# with different thinning intervals and proportion of the
# community thinned.
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
      paste0("Data/terrestrial_sim_data/lnorm_ricker/freq_x_nsp_results/", args[2])
    ),
    full.names = T
  )
  # put the file list in order
  sort_pat <- as.numeric(
    stringr::str_extract(file_list, "thintests_[:digit:]+") %>%
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

# create dataframe with thinning interval and proportion of community thinned
  conf_df_full <- data.frame(
    thin_interval = purrr::map_dbl(
      sims,
      ~ .x$sim_params$thin_freq
    ),
    prop_community = purrr::map_dbl(
      sims,
      ~ .x$sim_params$prop_cthin
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

# add rows for the missing combos to make the heatmap square
  empties <- data.frame(
    thin_interval = rep(c(rep(0, 10), unique(conf_df_sum$thin_interval)[-1]), each = 4),
    prop_community = rep(c(unique(conf_df_sum$prop_community)[-1], rep(0, 10)), each = 4),
    metric = rep(unique(conf_df_sum$metric), 20),
    value = NA,
    low = NA,
    high = NA
  )
  conf_df_sum <- rbind(
    conf_df_sum[1:4, ],
    empties,
    conf_df_sum[5:nrow(conf_df_sum), ]
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

  dat_nest <- conf_df_sum %>%
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

  # save the plot
  ggsave(
    ylabel + plots[[1]] + plots[[2]] + plots[[3]] + plots[[4]] + xlabel +
      plot_layout(design = lo, widths = unit(c(0.75, 4, 4), "cm"), heights = unit(c(4, 4, 0.75), "cm")),
    filename = here(
      paste0("Figures/confusion_metrics_freq_x_nsp_thin_", args[2], ".png")
    ),
    height = 7, width = 7,
    units = "in"
  )







