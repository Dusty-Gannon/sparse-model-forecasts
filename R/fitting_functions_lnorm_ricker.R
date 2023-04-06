###########################################
# Functions to fit and summarize the models
###########################################


##### Function to fit growth models #####

#' Fitting log-normal growth models
#'
#' @param N Matrix of community abundances, one species per row, with time over columns
#' @param stan_mod compiled stan model to fit
#' @param tsteps Vector of time steps to use
#' @param dist Logical indicating whether a covariate for disturbance times should be
#' included in the model.
#' @param ... Extra arguments to the control list for Stan
#'
#' @return Matrix of posterior draws for the competition coefficients
#'
fit_growth_models <- function(N, stan_mod, tsteps, dist_vec = NULL, ...){

  # get list of control arguments
  cntrl_args <- list(...)

  # compile data to feed into Stan
  N_het <- t(N[-1, tsteps])
  N_het_std <- scale(N_het)
  N_foc <- as.double(N[1, tsteps])
  n <- length(tsteps)
  y <- log(N_foc[2:n] / N_foc[1:(n - 1)])

  # determine if the focal went extinct at any point and
  # create an indicator vector for when present in the
  # community
  z <- 1 - as.numeric(
    is.infinite(y) | is.nan(y)
  )

  # replace the ys that will not get used
  if(sum(z) < length(z)){
    y[which(z == 0)] <- 1
  }

  # determine global shrinkage prior
  tau_0 <- tau0(
    y = y,
    m0 = min(5, ncol(N_het) - 1),
    M = ncol(N_het),
    N = n - 1,
    fam = "gaussian"
  )

  # compile non-shrinking variables
  if(is.null(dist_vec)){
    X_alpha <- matrix(
      data = as.double(scale(N_foc[1:(n - 1)])),
      ncol = 1
    )
  } else{
    X_alpha <- cbind(
      as.double(scale(N_foc[1:(n - 1)])),
      dist_vec[tsteps][1:(n - 1)]
    )
  }

  datlist <- list(
    N = n - 1,
    P0 = ncol(X_alpha),
    P = ncol(N_het),
    y = y,
    z = z,
    dens_foc = N_foc[1:(n - 1)],
    X_alpha = X_alpha,
    X_beta = N_het_std[1:(n - 1), ],
    error_scl = 0.5,
    tau0 = tau_0,
    slab_scl = 0.25,
    slab_df = 6
  )

  if(is.null(cntrl_args)){
    mfit <- rstan::sampling(
      stan_mod,
      data = datlist,
      cores = 3, chains = 3,
      control = list(adapt_delta = 0.99, max_treedepth = 15)
    )
  } else {
    mfit <- rstan::sampling(
      stan_mod,
      data = datlist,
      cores = 3, chains = 3,
      control = cntrl_args
    )
  }


  beta_post <- rstan::extract(mfit, pars = "beta")$beta

  return(beta_post %*% solve(diag(apply(N_het, 2, sd))))

}


##### Function to summarize results from posterior betas #####

#' Summarize a matrix of posterior draws for the competition coefficients
#'
#' @param beta_post Matrix of posterior draws
#' @param A_mat Original matrix of competition coefficients used to simulate the data
#' @param pip Posterior inclusion probability
#'
#' @return List with a confusion matrix and posterior RMSE for parameter draws
#'
conf_mat_summaries <- function(beta_post, A_mat, pip = 0.9){

  nz_estims <- which(apply(beta_post, 2, function(x){
    mean(x < 0) >= pip | mean(x > 0) >= pip
  }))
  z_estims <- which(!(1:ncol(beta_post) %in% nz_estims))
  nz_true <- which(A_mat[1, -1] != 0)
  z_true <-  which(!(1:ncol(beta_post) %in% nz_true))

  conf_mat <- matrix(
    c(
      sum(nz_estims %in% nz_true),
      sum(nz_estims %in% z_true),
      sum(z_estims %in% nz_true),
      sum(z_estims %in% z_true)
    ),
    nrow = 2, ncol = 2
  )

  rmse <- RMSE_bayes(
    as.double(A_mat[1, -1]),
    ppreds = beta_post
  )

  return(
    list(
      conf_mat = conf_mat,
      rmse = rmse
    )
  )

}

##### Wrapper function to fit and summarize #####


#' Wrapper to fit and summarize log-normal growth models
#'
#' @param X List with matrix of abundances called N in it
#' @param stan_mod Compiled stan model to fit
#' @param tsteps Vector of time steps to use as data
#' @param pip Posterior inclusion probability
#' @param dist Logical indicating whether disturbances that affect the
#' focal species should be accounted for in the model
#'
#' @return List with summaries and posterior draws
#'
fit_n_summarize <- function(X, stan_mod, tsteps = 51:100, pip = 0.9, dist = F){

  if(dist){
    if(X$sim_params$dist_prob > 0){
      beta_post <- fit_growth_models(
        N = X$N,
        stan_mod = stan_mod,
        tsteps = tsteps,
        dist_vec = X$dist_foc[2:length(X$dist_foc)]
      )
    }
  } else{
    beta_post <- fit_growth_models(
      N = X$N,
      stan_mod = stan_mod,
      tsteps = tsteps
    )
  }

  summaries <- list(
    beta_post = beta_post,
    conf_summaries = conf_mat_summaries(
      beta_post = beta_post,
      A_mat = X$sim_params$A_mat,
      pip = pip
    )
  )

  return(summaries)

}
