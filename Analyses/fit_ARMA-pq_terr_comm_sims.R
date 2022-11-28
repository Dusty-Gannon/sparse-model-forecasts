#########################################################
# Fitting ARMA(p,q) models to simulated community data
#########################################################


# libraries
  library(rstan)
  library(here)
  library(parallel)
 # devtools::load_all()

# load user-defined functions (teton way, not devtools way)
  src_files <- list.files(here("R/"), pattern = "*.R", full.names = T)
  sapply(src_files, source, .GlobalEnv)

# pull in arguments from the command line
  args <- commandArgs(trailingOnly = T)
  n_obs <- as.numeric(args[1])
  start <- as.numeric(args[2])
  p <- as.numeric(args[3])
  q <- as.numeric(args[4])

# function to fit the model
  fit_ARMA_pq <- function(
    X, dat_all, n_obs, p, q,
    start = 300, ann = as.character(1:5),
    pip = 0.7,
    S_init = 50,
    diagnostic_plots_dir = "Data/diagnostic_plots"
  ){

    # store some useful variables
    dat <- dat_all[[X]]
    dat_id <- X
    n <- start + n_obs
    m <- max(p, q)

    # convert cover data to wide format
    sp <- unique(as.character(dat$cover$species))
    cover_wide <- matrix(
      data = unique(dat$cover$t),
      ncol = 1
    )
    for(s in sp){
      cover_wide <- cbind(
        cover_wide,
        subset(dat$cover, species == s)$cover
      )
    }
    pcover_wide <- cbind(
      cover_wide[, 1],
      cover_wide[, -1] * 100
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
    col_ids_p <- which(!(names(cover_df) %in% c(ann, "t")))
    foc_p <- myTryCatch(choose_focal(
      df = cover_df,
      col_ids = col_ids_p,
      t = n, num_ngs = 2
    ))

    # define time steps
    tsteps_p1 <- (start + 1):n
    tsteps_m1 <- start:(n - 1)

    # proceed to model fitting as long as focal species exist
    if(!is.null(foc_a$value) | !is.null(foc_p$value)){
      # compile stan model
      arma_pq_FHS <- stan_model(here("Stan/ARMA-p-q_FHS.stan"))
    }

    ### Fitting the model for the annual focal species ###
    if(!is.null(foc_a$value)){
      tau_0_a <- tau0(
        y = cover_df[tsteps_p1, foc_a$value],
        m0 = 2,
        M = ncol(cover_df) - 2,
        N = n_obs,
        fam = "gaussian"
      )

      # define the model matrix
      X_beta_a <- as.matrix(scale(
        cover_df[tsteps_m1, -which(colnames(cover_df) %in% c("t", "env", foc_a$value))]
      ))

      # some species may be extinct towards the end of the time series
      extincts_a <- which(
        apply(X_beta_a, 2, function(x){
          sum(is.nan(x))
        }) > 0
      )
      X_beta_a <- X_beta_a[, -extincts_a]

      # create model matrix for non-shrinking effects
      env_std <- as.double(scale(cover_df$env))
      X_alpha <- cbind(
        rep(1, n_obs),
        env_std[tsteps_m1],
        (env_std[tsteps_m1])^2
      )

      # compile data list for fitting the model
      datlist_a <- list(
        N = n_obs,
        P0 = ncol(X_alpha),
        P = ncol(X_beta_a),
        p = p,
        q = q,
        y = as.double(scale(cover_df[tsteps_p1, foc_a$value])),
        X_alpha = X_alpha,
        X_beta = X_beta_a,
        tau0 = tau_0_a,
        slab_scl = 2,
        slab_df = 10
      )

      # fit the model
      armafit_a <- sampling(
        arma_pq_FHS,
        data = datlist_a,
        control = list(adapt_delta = 0.99, max_treedepth = 15)
      )

      # check the residuals
      y_rep_a <- rstan::extract(armafit_a, pars = "y_rep")$y_rep
      resids_a <- t(apply(
        y_rep_a, MARGIN = 1,
        FUN = function(x, obs){
          obs - x
        },
        obs = datlist_a$y
      ))

      df_resid_a <- data.frame(
        t = (m + 1):n_obs,
        mean = apply(resids_a[, (m + 1):n_obs], 2, mean),
        low = apply(resids_a[, (m + 1):n_obs], 2, quantile, probs = 0.05),
        high = apply(resids_a[, (m + 1):n_obs], 2, quantile, probs = 0.95)
      )

      resids_plot_a <- ggplot(df_resid_a, aes(x = t, y = mean))+
        geom_errorbar(aes(ymin = low, ymax = high), width = 0)+
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
      ggsave(
        filename = fname_a, plot = resids_plot_a,
        device = "png", width = 4, height = 3,
        units = "in"
      )

      # define the true positive betas
      comp_sp_a <- as.integer(foc_a$value) + c(1,2)
      comp_sp_id_a <- which(colnames(X_beta_a) %in% comp_sp_a)

      # define the true negative betas
      noncomp_sp_id_a <- which(!(colnames(X_beta_a) %in% comp_sp_a))

      # compute confusion matrix
      beta_post_a <- rstan::extract(armafit_a, pars = "beta")$beta
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
        fit = armafit_a,
        confusion_mat = confusion_mat_a,
        perc_diverged = mean(get_divergent_iterations(armafit_a))
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
        fam = "gaussian"
      )

      # define the model matrix
      X_beta_p <- as.matrix(scale(
        cover_df[tsteps_m1, -which(colnames(cover_df) %in% c("t", "env", foc_p$value))]
      ))

      # some species may be extinct towards the end of the time series
      extincts_p <- which(
        apply(X_beta_p, 2, function(x){
          sum(is.nan(x))
        }) > 0
      )
      X_beta_p <- X_beta_p[, -extincts_p]

      # create model matrix for non-shrinking effects
      env_std <- as.double(scale(cover_df$env))
      X_alpha <- cbind(
        rep(1, n_obs),
        env_std[tsteps_m1],
        (env_std[tsteps_m1])^2
      )

      # compile data list for fitting the model
      datlist_p <- list(
        N = n_obs,
        P0 = ncol(X_alpha),
        P = ncol(X_beta_p),
        p = p,
        q = q,
        y = as.double(scale(cover_df[tsteps_p1, foc_p$value])),
        X_alpha = X_alpha,
        X_beta = X_beta_p,
        tau0 = tau_0_p,
        slab_scl = 2,
        slab_df = 10
      )

      # fit the model
      armafit_p <- sampling(
        arma_pq_FHS,
        data = datlist_p,
        control = list(adapt_delta = 0.99, max_treedepth = 15)
      )

      # check the residuals
      y_rep_p <- rstan::extract(armafit_p, pars = "y_rep")$y_rep
      resids_p <- t(apply(
        y_rep_p, MARGIN = 1,
        FUN = function(x, obs){
          obs - x
        },
        obs = datlist_p$y
      ))

      df_resid_p <- data.frame(
        t = m:n_obs,
        mean = apply(resids_p, 2, mean),
        low = apply(resids_p, 2, quantile, probs = 0.05),
        high = apply(resids_p, 2, quantile, probs = 0.95)
      )

      resids_plot_p <- ggplot(df_resid_p, aes(x = t, y = mean))+
        geom_errorbar(aes(ymin = low, ymax = high), width = 0)+
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
      if(sum(as.integer(foc_p$value) + c(1, 2) > S_init) == 0){
        comp_sp_p <- as.integer(foc_p$value) + c(1, 2)
      } else if(sum(as.integer(foc_p$value) + c(1, 2) > S_init) == 1){
        comp_sp_p <- as.integer(foc_p$value) + c(-1, 1)
      } else{
        comp_sp_p <- as.integer(foc_p$value) - c(1, 2)
      }
      comp_sp_id_p <- which(colnames(X_beta_p) %in% comp_sp_p)

      # define the true negative betas
      noncomp_sp_id_p <- which(!(colnames(X_beta_p) %in% comp_sp_p))

      # compute confusion matrix
      beta_post_p <- rstan::extract(armafit_p, pars = "beta")$beta
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
      # define the annual component of the return list
      per_ret_list <- list(
        fit = armafit_p,
        confusion_mat = confusion_mat_p,
        perc_diverged = mean(get_divergent_iterations(armafit_p))
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


# load data
  dat_list <- readRDS(
    here("Data/terrestrial_sim_data/simdat_500reps_500steps_S2s48_5ann.rds")
  )

  # use parallel package to fit the models

  cl <- makeCluster(20)

  clusterEvalQ(cl, {library(rstan); library(here)})
  clusterEvalQ(
    cl,
    {
      src_files <- list.files(here("R/"), pattern = "*.R", full.names = T);
      sapply(src_files, source, .GlobalEnv)
      }
  )

  # create file path for diagnostic plots
  fp <- paste0("Data/terrestrial_sim_data/ARMA_", p, q, "_n", n_obs, "/diagnostic_plots/")

  results <- parLapply(
    cl = cl,
    X = 1:length(dat_list),
    fun = fit_ARMA_pq,
    n_obs = n_obs, p = p, q = q,
    diagnostic_plots_dir = fp,
    start = start
  )

  stopCluster(cl)

  fp_fits <- paste0("Data/terrestrial_sim_data/ARMA_", p, q, "_n", n_obs, "/model_fits.rds")

  # save the model fits
  saveRDS(results, file = here(fp_fits))







