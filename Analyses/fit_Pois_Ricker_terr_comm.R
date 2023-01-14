###################################################
# Fitting Ricker models to simulated community data
###################################################


# libraries
  library(rstan)
  library(here)
  library(parallel)
  devtools::load_all()

# # load user-defined functions (teton way, not devtools way)
#   src_files <- list.files(here("R/"), pattern = "*.R", full.names = T)
#   sapply(src_files, source, .GlobalEnv)


# function to fit the model
  fit_ricker <- function(
    X, n_obs,
    stan_mod,
    start = 100,
    pip = 0.7,
    S_init = 50, num_ngs = 3,
    lambda_knowledge = c("full", "covariates", "naive"),
    vary_comp = FALSE,
    ...
  ){

    xtra_args <- list(...)
    if(is.null(xtra_args$beta0_scl)){xtra_args$beta0_scl <- 0}
    if(is.null(xtra_args$a)){xtra_args$a <- 2}
    if(is.null(xtra_args$b)){xtra_args$b <- 2}
    # store some useful variables
    dat <- X
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
    if(isTRUE(vary_comp)){
      foc_sp <- myTryCatch(choose_focal2(
        df = dat$abundance_wide[,-c(1:3)],
        num_ngs = num_ngs, dyn_sp = dat$dynamic_sp,
        tw = c(start, n)
      ))
    } else{
      foc_sp <- myTryCatch(choose_focal(
        df = dat$abundance_wide[, -c(1, 2)],
        num_ngs = num_ngs, start = start,
        tw = n_obs
      ))
    }

    ### Fitting the model for the focal species ###
    if(!is.null(foc_sp$value)){

      foc <- paste0("s", foc_sp$value$foc)
      tsteps <- foc_sp$value$tw
      df <- dat$abundance_wide[tsteps, ]

      if(lambda_knowledge == "full"){
        # construct known lambda
        lambda <- gauss_env_effect(
          env = df$v1,
          lambda_max = dat$lambda_max[foc_sp$value$foc],
          optims = dat$sp_optims[foc_sp$value$foc]
        )

        # columns not to be included
        xclude <- which(names(df) %in% c("t", "v1", "v2", foc))
        X_beta0 <- as.matrix(df[, -xclude])
        X_beta0 <- X_beta0[, -which(apply(X_beta0, 2, mean) == 0)] /
          mean(df[, foc])

        # define design matrix for beta depending on whether the
        # competition varies through time
        if(vary_comp){
          X_beta <- cbind(
            X_beta0,
            apply(X_beta0, 2, function(x, v){x * v}, v = df$v2)
          )
          # rename columns to have unique names
          colnames(X_beta) <- c(
            colnames(X_beta0),
            paste0("v", colnames(X_beta0))
          )
        } else{
          X_beta <- X_beta0
        }

        # define tau0
        tau_0 <- tau0(
          y = df[, foc],
          m0 = 5,
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
          beta0_scl = xtra_args$beta0_scl,
          tau0 = tau_0,
          slab_scl = 1,
          slab_df = 8
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

        if(vary_comp){
          X_beta <- cbind(
            X_beta0,
            apply(X_beta0, 2, function(x, v){x * v}, v = df$v2)
          )
          # rename columns to have unique names
          colnames(X_beta) <- c(
            colnames(X_beta0),
            paste0("v", colnames(X_beta0))
          )
        } else{
          X_beta <- X_beta0
        }

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
          beta0_scl = xtra_args$beta0_scl,
          tau0 = tau_0,
          slab_scl = 1,
          slab_df = 8,
          lambda_star = dat$lambda_max[foc_sp$value$foc]
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

        if(vary_comp){
          X_beta <- cbind(
            X_beta0,
            apply(X_beta0, 2, function(x, v){x * v}, v = df$v2)
          )
          # rename columns to have unique names
          colnames(X_beta) <- c(
            colnames(X_beta0),
            paste0("v", colnames(X_beta0))
          )
        } else{
          X_beta <- X_beta0
        }

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
          beta0_scl = xtra_args$beta0_scl,
          tau0 = tau_0,
          slab_scl = 1,
          slab_df = 8
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
      ng_species <- which(
        (dat$B_mat[, 2, foc_sp$value$foc] != 0)
      )
      ng_species <- ng_species[-which(
        ng_species == foc_sp$value$foc
      )]
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
          foc_sp = foc_sp$value$foc,
          nonzero_betas = ng_betas,
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




  ####################
  # fitting the models
  ####################

  # pull in arguments from the command line
  args <- commandArgs(trailingOnly = T)
  datfile <- args[1]
  n_obs <- as.numeric(args[2])
  start <- as.numeric(args[3])

  # load data on director node
  dat_list <- readRDS(
    here(paste0("Data/terrestrial_sim_data/", datfile))
  )[1]

  # compile stan model
  modfile <- paste0("Stan/Pois_ricker_", args[4], "_lambda_FHS.stan")
  stan_mod <- rstan::stan_model(here(modfile))

  # set knowledge parameter
  if(args[4] == "known"){
    lambda_knowledge <- "full"
  }
  if(args[4] == "partial"){
    lambda_knowledge <- "covariates"
  }
  if(args[4] == "fixed"){
    lambda_knowledge <- "naive"
  }

  # use parallel package to fit the models
  cl <- makeCluster(4)

  # load libraries and functions on each node
  clusterEvalQ(cl, {library(rstan); library(here); devtools::load_all()})
  # clusterEvalQ(
  #   cl,
  #   {
  #     src_files <- list.files(here("R/"), pattern = "*.R", full.names = T);
  #     sapply(src_files, source, .GlobalEnv)
  #   }
  # )

  results <- parLapply(
    cl = cl,
    X = dat_list,
    fun = fit_ricker,
    n_obs = n_obs,
    stan_mod = stan_mod,
    pip = 0.65,
    start = start,
    lambda_knowledge = lambda_knowledge,
    a = 2, b = 2
  )

  stopCluster(cl)

  # save fits into unique directory and filename
  dir_fits <- paste0("Data/terrestrial_sim_data/Ricker_lambda_", args[3], "/")
  dir.create(here(dir_fits))

  fname <- paste0("fits_n", n_obs, "_lambda_", args[3], ".rds")

  fits_file <- paste0(dir_fits, fname)

  # save the model fits
  saveRDS(results, file = here(fits_file))







