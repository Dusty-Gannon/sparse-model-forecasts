# libraries
library(here)
library(rstan)
library(parallel)
library(dplyr)
library(tidyr)

# # load functions on personal machine
# devtools::load_all()

# load user-defined functions (teton way, not devtools way)
src_files <- list.files(here("R/"), pattern = "*.R", full.names = T)
sapply(src_files, source, .GlobalEnv)

# function to compute posterior mode

# compile the stan model
  pois_garma11 <- stan_model(here("Stan/Pois_GLARMA-1-1_simple.stan"))

# function to simulate from a GLARMA(1,1) process and fit the model with
# 10 different fixed values of c and one random for augmenting the raw data

  # set parameter values
  params = list(
    n = 100,
    beta = c(-0.5, 1, 0.5),
    phi = 0.6,
    theta = 0.5
  )

threshold_eval <- function(X, stan_mod, params){

  # handy function for posterior mode

  # simulating from a GLARMA(1,1) process

  # set parameters
    # time series length
    n <- params$n

    # reg. coefficients
    beta <- params$beta

    # AR and MA parameters
    phi <- params$phi
    theta <- params$theta

    # model matrix
    X <- matrix(
      rnorm(n = n * (length(beta) - 1)),
      nrow = n
    )
    X <- cbind(
      rep(1, n),
      X
    )

    # threshold parameters
    thresh <- seq(0.05, 0.95, length.out = 10)

  # initiate the process
    y <- vector(mode = "double", length = n)
    y_star_true <- vector(mode = "double", length = n)
    mu <- vector(mode = "double", length = n)
    z <- vector(mode = "double", length = n)
    y[1] <- rpois(1, lambda = exp(X[1, ] %*% beta))
    y_star_true[1] <- max(y[1], runif(1))
    mu[1] <- y_star_true[1]
    z[1] <- 0

  # continue the process
    for(t in 2:n){

      # arma term
      z_t <- phi * (log(y_star_true[t-1]) - X[t-1, ] %*% beta) + theta * log(y_star_true[t-1]/mu[t-1])
      mu[t] <- exp(X[t, ] %*% beta + z_t)
      y[t] <- rpois(1, lambda = mu[t])
      y_star_true[t] <- max(y[t], runif(1))

    }

  # now define y_star according to different thresholds
    y_star <- vector(mode = "list", length = length(thresh) + 1)
    for(j in 1:length(thresh)){
      y_star[[j]] <- sapply(
        y,
        function(x, x2){
          max(x, x2)
        },
        x2 = thresh[j]
      )
    }
    # randomized threshold
    y_star[[length(y_star)]] <- sapply(
      1:n,
      function(i, x, x2){
        max(x[i], x2[i])
      },
      x = y,
      x2 = runif(n)
    )

    # compile data EXCEPT for y_star
    dat_remdat <- list(
      N = n,
      P = length(beta),
      y = y,
      X = X
    )

    # function to fit the model
    fit <- function(x, remain_dat, mod){

      # compile data to feed into stan
      datlist <- list(
        N = remain_dat$N,
        P = remain_dat$P,
        y = remain_dat$y,
        X = remain_dat$X,
        y_star = x
      )

      # fit mod
      mfit <- rstan::sampling(
        mod,
        data = datlist,
        cores = 3,
        chains = 3
      )

      return(mfit)

    }


    # apply the fit to each y_star
    fits_all <- lapply(
      y_star,
      FUN = fit,
      remain_dat = dat_remdat,
      mod = stan_mod
    )

    # extract the posterior modes for each parameter of interest
    res <- vector(mode = "list")
    res$phi <- sapply(
      fits_all,
      FUN = function(x){
        phi_post <- rstan::extract(x, pars = "phi")$phi
        return(post_mode(phi_post))
      }
    )
    res$theta <- sapply(
      fits_all,
      FUN = function(x){
        theta_post <- rstan::extract(x, pars = "theta")$theta
        return(post_mode(theta_post))
      }
    )
    res$beta <- sapply(
      fits_all,
      FUN = function(x){
        beta_post <- rstan::extract(x, pars = "beta")$beta
        return(
          apply(beta_post, 2, post_mode)
        )
      }
    )

    return(res)
}


# use parallel package to repeat the simulations 500 times
  reps <- 500

  cl <- makeCluster(20)

  clusterEvalQ(cl, {library(rstan)})
  clusterEvalQ(cl, {
    post_mode <- function(x){
     id <- which.max(density(x)$y)
     return(density(x)$x[id])
  }
  })

  results <- parLapply(
    cl = cl,
    X = 1:reps,
    fun = threshold_eval,
    stan_mod = pois_garma11,
    params = params
  )

  stopCluster(cl)


 saveRDS(results, file = here("Data/Pois_GLARMA_thresh_eval/pois_GLARMA_thresh_evals.rds"))




#### Figure to visualize results ####
  df_results <- data.frame(
    c = rep(
      c(seq(0.05, 0.95, length.out = 10), "random"),
      length(results)
    ),
    phi = NA,
    theta = NA,
    beta0 = NA,
    beta1 = NA,
    beta2 = NA
  )

  # now loop through and add the appropriate values
  for(i in 1:length(results)){

    row_ids <- (11 * (i - 1) + 1):(11 * (i - 1) + 11)

    df_results$phi[row_ids] <- results[[i]]$phi
    df_results$theta[row_ids] <- results[[i]]$theta
    df_results$beta0[row_ids] <- results[[i]]$beta[1, ]
    df_results$beta1[row_ids] <- results[[i]]$beta[2, ]
    df_results$beta2[row_ids] <- results[[i]]$beta[3, ]

  }

  # long format
  df_results_long <- pivot_longer(
    df_results,
    cols = phi:beta2,
    names_to = "param",
    values_to = "estimate"
  )

  # summarize results for each threshold value
  df_plot <- group_by(df_results_long, param, c) %>%
    summarise(
      estim = mean(estimate),
      se = sd(estimate)
    )

  # create another dataframe with true values for each parameter
  df_true <- data.frame(
    param = unique(df_plot$param),
    true = c(params$beta, params$phi, params$theta)
  )

  # create facet labels
  param_labs <- c(
    beta0 = "beta[0]",
    beta1 = "beta[1]",
    beta2 = "beta[2]",
    phi = "phi",
    theta = "theta"
  )


  # create the plot
  p <- ggplot(data = df_plot, aes(x = c))+
    geom_errorbar(aes(ymin = estim - se, ymax = estim + se))+
    geom_point(aes(y = estim))+
    geom_hline(
      data = df_true, aes(yintercept = true),
      color = "brown2", linetype = "dashed"
    )+
    theme_bw()+
    theme(
      axis.text.x = element_text(angle = 45, vjust = 0.5)
    )+
    xlab("Threshold value")+
    ylab("Posterior mode")+
    facet_wrap(
      ~ param,
      labeller = labeller(param = as_labeller(param_labs, label_parsed))
    )

  ggsave(
    filename = here("Figures/Pois_GLARMA_threshold_evals.png"),
    plot = p, width = 6, height = 4, units = "in",
    device = "png"
  )








