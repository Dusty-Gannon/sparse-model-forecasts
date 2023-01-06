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
    start = 100,
    pip = 0.7,
    S_init = 30,
    lambda_knowledge = c("full", "covariates", "naive"),
    ...
  ){

    xtra_args <- list(...)
    # store some useful variables
    dat <- dat_all[[X]]
    dat_id <- X
    n <- start + n_obs
    nsp <- length(unique(dat$abundance_long$species))

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

    # find species with at least one strong competitor, preferably
    #  a competitor that is also dynamic
    foc_sp <- myTryCatch(choose_focal2(
      df = dat$abundance_wide[,-c(1:3)],
      num_ngs = 2, dyn_sp = dat$dynamic_sp,
      tw = c(start, n)
    ))

    # define time steps
    tsteps <- (start + 1):n

    ### Fitting the model for the focal species ###
    if(!is.null(foc_sp$value)){

      foc <- paste0("s", foc_sp$value)
      df <- dat$abundance_wide[tsteps, ]

      if(lambda_knowledge == "full"){
        # construct known lambda
        lambda <- gauss_env_effect(
          env = df$v1,
          lambda_max = dat$lambda_max[foc_sp$value],
          optims = dat$sp_optims[foc_sp$value]
        )

        # # define the model matrices
        # X_lambda <- model.matrix(~ v1 + I(v1^2), data = df)

        X_beta0 <- as.matrix(
          df[, -which(names(df) %in% c("t", "v1", "v2", foc))]
        )
        X_beta0 <- X_beta0[, -which(apply(X_beta0, 2, mean) == 0)] /
          mean(df[, foc])

        X_beta <- cbind(
          X_beta0,
          apply(X_beta0, 2, function(x, v){x * v}, v = df$v2)
        )

        # define tau0
        tau_0 <- tau0(
          y = df[, foc],
          m0 = 2,
          M = ncol(X_beta),
          N = n_obs,
          fam = "poisson"
        )

        # compile data list for fitting the model
        datlist <- list(
          N = n_obs,
          P_h = ncol(X_beta0),
          P = ncol(X_beta),
          scl_X = mean(df[, foc]),
          y = as.integer(df[, foc]),
          lambda = lambda,
          X_beta0 = X_beta0,
          X_beta = X_beta,
          tau0 = tau_0,
          slab_scl = 0.2,
          slab_df = 10
        )

        # fit the model
        ricker_fit <- rstan::sampling(
          stan_mod,
          data = datlist,
          cores = 3,
          chains = 3,
          control = list(adapt_delta = 0.99, max_treedepth = 15)
        )
      }

      # model fit when information on lambda is partial
      if(lambda_knowledge == "covariates"){

        # define the model matrices
        X_lambda <- model.matrix(~ v1 + I(v1^2), data = df)

        X_beta0 <- as.matrix(
          df[, -which(names(df) %in% c("t", "v1", "v2", foc))]
        )
        X_beta0 <- X_beta0[, -which(apply(X_beta0, 2, mean) == 0)] /
          mean(df[, foc])

        X_beta <- cbind(
          X_beta0,
          apply(X_beta0, 2, function(x, v){x * v}, v = df$v2)
        )

        tau_0 <- tau0(
          y = df[, foc],
          m0 = 2,
          M = ncol(X_beta),
          N = n_obs,
          fam = "poisson"
        )

        # compile data list for fitting the model
        datlist <- list(
          N = n_obs,
          P_lambda = ncol(X_lambda),
          P_h = ncol(X_beta0),
          P = ncol(X_beta),
          scl_X = mean(df[, foc]),
          y = as.integer(df[, foc]),
          X_lambda = X_lambda,
          X_beta0 = X_beta0,
          X_beta = X_beta,
          tau0 = tau_0,
          slab_scl = 0.5,
          slab_df = 10,
          lambda_star = dat$lambda_max[foc_sp$value]
        )

        # fit the model
        ricker_fit <- rstan::sampling(
          stan_mod,
          data = datlist,
          cores = 3,
          chains = 3,
          control = list(adapt_delta = 0.99, max_treedepth = 15)
        )

      }

      # fit the model with a prior on lambda
      if(lambda_knowledge == "naive"){

        hyper_pars <- c(xtra_args$a, xtra_args$b)

        # construct the model matrices
        X_beta0 <- as.matrix(
          df[, -which(names(df) %in% c("t", "v1", "v2", foc))]
        )
        X_beta0 <- X_beta0[, -which(apply(X_beta0, 2, mean) == 0)] /
          mean(df[, foc])

        X_beta <- cbind(
          X_beta0,
          apply(X_beta0, 2, function(x, v){x * v}, v = df$v2)
        )
        # rename columns to have unique names
        colnames(X_beta) <- c(
          colnames(X_beta0),
          paste0("v", colnames(X_beta0))
        )

        tau_0 <- tau0(
          y = df[, foc],
          m0 = 2,
          M = ncol(X_beta),
          N = n_obs,
          fam = "poisson"
        )

        # compile data list for fitting the model
        datlist <- list(
          N = n_obs,
          P_h = ncol(X_beta0),
          P = ncol(X_beta),
          scl_X = mean(df[, foc]),
          a_lambda = hyper_pars[1], b_lambda = hyper_pars[2],
          y = as.integer(df[, foc]),
          X_beta0 = X_beta0,
          X_beta = X_beta,
          tau0 = tau_0,
          slab_scl = 0.2,
          slab_df = 10
        )

        # fit the model
        ricker_fit <- rstan::sampling(
          stan_mod,
          data = datlist,
          cores = 3,
          chains = 3,
          control = list(adapt_delta = 0.99, max_treedepth = 15)
        )

      }

      # define the true positive betas
      ng_species <- as.integer(foc_sp$value) + c(1:dat$sim_params$S)
      ng_coefs <- paste0("s", ng_species)

      # define whether the ng_competitors have dynamic effects
      for(i in 1:length(ng_species)){
        if(ng_species[i] %in% dat$dynamic_sp){
          ng_coefs <- c(
            ng_coefs,
            paste0("vs", ng_species[i])
          )
        }
      }
      ng_betas <- which(colnames(X_beta) %in% ng_coefs)

      # define the true zero betas
      g_betas <- which(!(colnames(X_beta) %in% ng_coefs))

      # compute confusion matrix
      beta_post <- rstan::extract(ricker_fit, pars = "beta")$beta
      nzs <- apply(
        beta_post, MARGIN = 2,
        FUN = function(x, p){
          as.integer(
            {mean(x < 0) > p | mean(x < 0) < 1 - p}
          )
        },
        p = pip
      )
      confusion_mat <- matrix(
        data = c(
          sum(which(nzs == 1) %in% ng_betas),
          sum(which(nzs == 1) %in% g_betas),
          sum(which(nzs == 0) %in% ng_betas),
          sum(which(nzs == 0) %in% g_betas)
        ),
        ncol = 2, nrow = 2
      )
      # return the relevant objects
      if(lambda_knowledge == "full"){
        lambda_post <- lambda
      }
      if(lambda_knowledge == "naive"){
        lambda_post <- rstan::extract(ricker_fit, pars = "lambda")$lambda
      }
      if(lambda_knowledge == "covariates"){
        lambda_post <- rstan::extract(ricker_fit, pars = "gamma")$gamma
      }
      return(
        list(
          data = dat,
          fit = list(
            confusion_mat = confusion_mat,
            perc_diverged = mean(get_divergent_iterations(ricker_fit)),
            beta_post = beta_post,
            beta0_post = rstan::extract(ricker_fit, pars = "beta0")$beta0,
            alpha_post = rstan::extract(ricker_fit, pars = "alpha")$alpha,
            lambda_post = lambda_post
          )
        )
      )
    } else{
      return(NULL)
    }

  }




### fitting the models ###

  # load data on director node
  dat_list <- readRDS(
    here("Data/terrestrial_sim_data/simdat_100reps_500steps_S2s28_10dyn_2env.rds")
  )

  # compile stan model
  ricker1 <- stan_model(here("Stan/Pois_ricker_known_lambda_FHS.stan"))


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







