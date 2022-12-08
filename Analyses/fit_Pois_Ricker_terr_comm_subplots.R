###################################################
# Fitting Ricker models to simulated community data
###################################################


# libraries
  library(rstan)
  library(here)
  # library(parallel)
  devtools::load_all()

# load user-defined functions (teton way, not devtools way)
  src_files <- list.files(here("R/"), pattern = "*.R", full.names = T)
  sapply(src_files, source, .GlobalEnv)

# pull in arguments from the command line
  args <- commandArgs(trailingOnly = T)
  n_obs <- as.numeric(args[1])
  start <- as.numeric(args[2])

# function to fit the model
  fit_ricker <- function(
    X, dat_all, n_obs,
    stan_mod,
    start = 200, ann = as.character(1:5),
    pip = 0.7,
    S_init = 50,
    sub_cells = 7,
    diagnostic_plots_dir = "Data/diagnostic_plots/"
  ){

    # store some useful variables
    dat <- dat_all[[X]]
    dat_id <- X
    n <- start + n_obs

    # convert cover data to wide format
    sp <- unique(as.character(dat$subplot$species))
    cover_wide <- matrix(
      data = unique(dat$subplot$t),
      ncol = 1
    )
    for(s in sp){
      cover_wide <- cbind(
        cover_wide,
        subset(dat$subplot, species == s)$cover
      )
    }
    pcover_wide <- cbind(
      cover_wide[, 1],
      cover_wide[, -1] * sub_cells *sub_cells
    )
    colnames(pcover_wide) <- c("t", sp)
    cover_df <- as.data.frame(pcover_wide)

    # add the environmental covariate
    cover_df$env <- dat$env

    # add methods for handling errors
    e <- simpleError("No species in the dataset fits the bill...")
    myTryCatch <- function(expr) {
      warn <- err <- NULL
      value <- withCallingHandlers(
        tryCatch(expr, error=function(e) {
          err <<- e
          NULL
        }), warning=function(w) {
          warn <<- w
          invokeRestart("muffleWarning")
        })
      list(value = value, warning = warn, error = err)
    }

    # find annual with greatest abundance at the end that also has at least
    #  one strong competitor
    col_ids_a <- which(names(cover_df) %in% ann)
    foc_a <- myTryCatch(foc_a <- choose_focal(
      df = cover_df,
      col_ids = col_ids_a,
      t = n, num_ngs = 2,
      exclude_names = c("t", "env")
    ))

    # find a perennial with high abundance and at least one
    #  strong competitor
    col_ids_p <- which(!(names(cover_df) %in% c(ann, "t", "env")))
    foc_p <- myTryCatch(choose_focal(
      df = cover_df,
      col_ids = col_ids_p,
      t = n, num_ngs = 2,
      exclude_names = c("t", "env")
    ))

    # define time steps
    tsteps_p1 <- (start + 1):n
    tsteps_m1 <- start:(n - 1)

    ### Fitting the model for the annual focal species ###
    if(!is.null(foc_a$value)){
      tau_0_a <- tau0(
        y = cover_df[tsteps_p1, foc_a$value],
        m0 = 2,
        M = ncol(cover_df) - 2,
        N = n_obs,
        fam = "poisson"
      )

      # define the model matrix
      X_beta_a <- as.matrix(
        cover_df[tsteps_m1, -which(colnames(cover_df) %in% c("t", "env"))]
      )

      # some species may be extinct towards the end of the time series
      extincts_a <- which(
        apply(X_beta_a, 2, function(x){
          sd(x)
        }) == 0
      )
      if(length(extincts_a) > 0){
        X_beta_a <- X_beta_a[, -extincts_a]
      }

      # create model matrix for non-shrinking effects
      env <- cover_df$env
      X_alpha <- cbind(
        rep(1, n_obs),
        env[tsteps_m1],
        (env[tsteps_m1])^2
      )

      # define augmented data
      y_star <- as.integer(cover_df[tsteps_p1, foc_a$value])
      y_star[y_star == 0] <- 0.5

      # compile data list for fitting the model
      datlist_a <- list(
        N = n_obs,
        P0 = ncol(X_alpha),
        P = ncol(X_beta_a),
        y = as.integer(cover_df[tsteps_p1, foc_a$value]),
        y_star = y_star,
        X_alpha = X_alpha,
        X_beta = X_beta_a,
        tau0 = tau_0_a,
        slab_scl = 2,
        slab_df = 10
      )

      # fit the model
      ricker_fit_a <- sampling(
        stan_mod,
        data = datlist_a,
        cores = 3,
        chains = 3,
        control = list(adapt_delta = 0.99, max_treedepth = 15)
      )

      # check the residuals
      y_rep_a <- rstan::extract(ricker_fit_a, pars = "y_rep")$y_rep
      resids_a <- bayes_qresids(datlist_a$y, y_rep = y_rep_a)

      df_resid_a <- data.frame(
        t = (start + 1):n,
        resids_a
      )

      resids_plot_a <- ggplot(df_resid_a, aes(x = t, y = resids_a))+
        geom_point()+
        theme_bw()

      # save the residuals plot
      if(!dir.exists(here(diagnostic_plots_dir))){
        dir.create(here(diagnostic_plots_dir), recursive = T)
      }
      fname_a <- paste0(
        here(diagnostic_plots_dir),
        "resids_", dat_id, "a", ".png"
      )
      png(
        filename = fname_a,
        width = 4, height = 3,
        units = "in", res = 300
      )
        resids_plot_a
      dev.off()

      # define the true positive betas
      comp_sp_a <- as.integer(foc_a$value) + c(0,1,2)
      comp_sp_id_a <- which(colnames(X_beta_a) %in% comp_sp_a)

      # define the true negative betas
      noncomp_sp_id_a <- which(!(colnames(X_beta_a) %in% comp_sp_a))

      # compute confusion matrix
      beta_post_a <- rstan::extract(ricker_fit_a, pars = "beta")$beta
      nzs <- apply(
        beta_post_a, MARGIN = 2,
        FUN = function(x, p){
          as.integer(
            {mean(x < 0) > p | mean(x < 0) < 1 - p}
          )
        },
        p = pip
      )
      confusion_mat_a <- matrix(
        data = c(
          sum(which(nzs == 1) %in% comp_sp_id_a),
          sum(which(nzs == 1) %in% noncomp_sp_id_a),
          sum(which(nzs == 0) %in% comp_sp_id_a),
          sum(which(nzs == 0) %in% noncomp_sp_id_a)
        ),
        ncol = 2, nrow = 2
      )
      # define the annual component of the return list
      ann_ret_list <- list(
        confusion_mat = confusion_mat_a,
        perc_diverged = mean(get_divergent_iterations(ricker_fit_a))
      )
    } else{
      ann_ret_list = NULL
    }

    ### Fitting the model for the focal perennial ###
    if(!is.null(foc_p$value)){
      tau_0_p <- tau0(
        y = cover_df[tsteps_p1, foc_p$value],
        m0 = 2,
        M = ncol(cover_df) - 2,
        N = n_obs,
        fam = "poisson"
      )

      # define the model matrix
      X_beta_p <- as.matrix(
        cover_df[tsteps_m1, -which(colnames(cover_df) %in% c("t", "env"))]
      )

      # some species may be extinct towards the end of the time series
      extincts_p <- which(
        apply(X_beta_p, 2, function(x){
          sd(x)
        } == 0
      ))
      if(length(extincts_p) > 0){
        X_beta_p <- X_beta_p[, -extincts_p]
      }

      # create model matrix for non-shrinking effects
      env <- as.double(cover_df$env)
      X_alpha <- cbind(
        rep(1, n_obs),
        env[tsteps_m1],
        (env[tsteps_m1])^2
      )

      # compile data list for fitting the model
      y_star <- as.integer(cover_df[tsteps_p1, foc_p$value])
      y_star[y_star == 0] <- 0.5

      datlist_p <- list(
        N = n_obs,
        P0 = ncol(X_alpha),
        P = ncol(X_beta_p),
        y = as.integer(cover_df[tsteps_p1, foc_p$value]),
        y_star = y_star,
        X_alpha = X_alpha,
        X_beta = X_beta_p,
        tau0 = tau_0_p,
        slab_scl = 2,
        slab_df = 10
      )

      # fit the model
      ricker_fit_p <- sampling(
        stan_mod,
        data = datlist_p,
        cores = 3,
        chains = 3,
        control = list(adapt_delta = 0.99, max_treedepth = 15)
      )

      # check the residuals
      y_rep_p <- rstan::extract(ricker_fit_p, pars = "y_rep")$y_rep
      resids_p <- bayes_qresids(obs = datlist_p$y, y_rep = y_rep_p, z = T)

      df_resid_p <- data.frame(
        t = (start + 1):n,
        resids = resids_p
      )

      resids_plot_p <- ggplot(df_resid_p, aes(x = t, y = resids))+
        geom_point()+
        theme_bw()

      # save the residuals plot
      if(!dir.exists(here(diagnostic_plots_dir))){
        dir.create(here(diagnostic_plots_dir), recursive = T)
      }
      fname_p <- paste0(
        here(diagnostic_plots_dir),
        "resids_", dat_id, "p", ".png"
      )
      ggsave(
        filename = fname_p, plot = resids_plot_p,
        device = "png", width = 4, height = 3,
        units = "in"
      )

      # define the true positive betas
      if(sum(as.integer(foc_p$value) + c(0, 1, 2) > S_init) == 0){
        comp_sp_p <- as.integer(foc_p$value) + c(0, 1, 2)
      } else if(sum(as.integer(foc_p$value) + c(0, 1, 2) > S_init) == 1){
        comp_sp_p <- as.integer(foc_p$value) + c(-1, 0, 1)
      } else{
        comp_sp_p <- as.integer(foc_p$value) - c(0, 1, 2)
      }
      comp_sp_id_p <- which(colnames(X_beta_p) %in% comp_sp_p)

      # define the true negative betas
      noncomp_sp_id_p <- which(!(colnames(X_beta_p) %in% comp_sp_p))

      # compute confusion matrix
      beta_post_p <- rstan::extract(ricker_fit_p, pars = "beta")$beta
      nzs <- apply(
        beta_post_p, MARGIN = 2,
        FUN = function(x, p){
          as.integer(
            {mean(x < 0) > p | mean(x < 0) < 1 - p}
          )
        },
        p = pip
      )
      confusion_mat_p <- matrix(
        data = c(
          sum(which(nzs == 1) %in% comp_sp_id_p),
          sum(which(nzs == 1) %in% noncomp_sp_id_p),
          sum(which(nzs == 0) %in% comp_sp_id_p),
          sum(which(nzs == 0) %in% noncomp_sp_id_p)
        ),
        ncol = 2, nrow = 2
      )
      # define the perennial component of the return list
      per_ret_list <- list(
        confusion_mat = confusion_mat_p,
        perc_diverged = mean(get_divergent_iterations(ricker_fit_p))
      )
    } else{
      per_ret_list = NULL
    }

    # return each list of results
    return(list(
      annual = ann_ret_list,
      perennial = per_ret_list
    ))

  }




### fitting the models ###

  # load data on director node
  dat_list <- readRDS(
    here("Data/terrestrial_sim_data/simdat_500reps_500steps_S2s48_5ann_env0.05.rds")
  )

  # compile stan model
  ricker_mod <- stan_model(here("Stan/Pois_ricker_FHS.stan"))


  # # use parallel package to fit the models
  # cl <- makeCluster(12)
  #
  # # load libraries and functions on each node
  # clusterEvalQ(cl, {library(rstan); library(here)})
  # clusterEvalQ(
  #   cl,
  #   {
  #     src_files <- list.files(here("R/"), pattern = "*.R", full.names = T);
  #     sapply(src_files, source, .GlobalEnv)
  #     }
  # )

  # create file path for diagnostic plots
  fp <- paste0("Data/terrestrial_sim_data/Ricker_", n_obs, "_env0.05", "/diagnostic_plots/")

  results <- pblapply(
    X = 1:length(dat_list),
    FUN = fit_ricker,
    dat_all = dat_list,
    n_obs = n_obs,
    stan_mod = ricker_mod,
    pip = 0.65,
    diagnostic_plots_dir = fp,
    start = 200
  )

  # results <- vector(mode = "list", length = 100)
  # for(i in 1:100){
  #
  #   results[[i]] <- fit_ARMA_pq(
  #     dat = dat_list[[i]], n_obs = 100, p = 1,
  #     q = 1, stan_mod = arma_pq_FHS,
  #     dat_id = i,
  #     pip = 0.7, diagnostic_plots_dir = fp,
  #     start = 300
  #   )
  #
  # }

#  stopCluster(cl)

  fp_fits <- paste0("Data/terrestrial_sim_data/Ricker_", n_obs, "_env0.05", "/model_fits.rds")

  # save the model fits
  saveRDS(results, file = here(fp_fits))

  confmat_a <- matrix(data = 0, nrow = 2, ncol = 2)
  confmat_p <- confmat_a

  for(i in 1:length(results)){

    if(!is.null(results[[i]]$annual)){
      confmat_a <- confmat_a + results[[i]]$annual$confusion_mat
    }

    if(!is.null(results[[i]]$perennial)){
      confmat_p <- confmat_p + results[[i]]$perennial$confusion_mat
    }

  }







