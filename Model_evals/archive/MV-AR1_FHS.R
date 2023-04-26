### Testing the MVN time series model ###

# libraries
  library(rstan)
  library(abind)
  library(viridis)

# set simulation parameters

  # length of time series
  n <- 200

  # number of species
  S <- 4

  # means for each species
  mus <- runif(S, min = 2, max = 6)

  # standard deviations of innovations
  sigmas <- runif(S, min = 0.2, max = 0.8)

  # competition coeffs
  alphas <- runif(S, min = 0.3, max = 0.5)
  A_mat <- comp_matrix2(S, rho = 0, alpha = alphas, num_ngs = 2)

  # growth rates
  phi <- runif(S, min = 0.5)

  # initialize response matrix
  y <- matrix(nrow = S, ncol = n)
  y[, 1] <- rnorm(S, mean = mus, sd = sigmas)

  for(t in 2:n){

    y[, t] <- mus + phi * y[, (t - 1)] - A_mat %*% y[, (t - 1)] + rnorm(S, sd = sigmas)

  }

  df_cover <- as.data.frame(t(y))
  names(df_cover) <- c("s1", "s2", "s3", "s4")
  df_cover$t <- 1:n

  df_long <- tidyr::pivot_longer(
    df_cover, cols = s1:s4,
    names_to = "species",
    values_to = "cover"
  )

  ggplot(data = df_long, aes(x = t, y = cover, color = species))+
    geom_line()+
    theme_classic()


# can we fit the model?
  n_obs <- 100
  datlist <- list(
    N = n_obs,
    S = S,
    P = 1,
    y = y[, (n - n_obs + 1):n],
    X_alpha = array(data = 1, dim = c(n_obs, S, 1)),
    tau0 = 0.001,
    slab_scl = 1,
    slab_df = 10
  )

  mvar1 <- stan_model(file = here("Stan/MV-AR1_FHS.stan"))


  mfit <- sampling(
    mvar1,
    data = datlist,
    chains = 3,
    cores = 3,
    control = list(adapt_delta = 0.9)
  )

# posterior predictions
  y_hat <- matrix(nrow = S, ncol = n_obs)
  for(j in 1:S){
    for(t in 1:n_obs){
      y_hat[j, t] <- mean(y_rep[, t, j])
    }
  }

  y_low <- matrix(nrow = S, ncol = n_obs)
  for(j in 1:S){
    for(t in 1:n_obs){
      y_low[j, t] <- quantile(y_rep[, t, j], probs = 0.05)
    }
  }

  y_high <- matrix(nrow = S, ncol = n_obs)
  for(j in 1:S){
    for(t in 1:n_obs){
      y_high[j, t] <- quantile(y_rep[, t, j], probs = 0.95)
    }
  }

  y_hat_t <- as.data.frame(cbind(
    t(y_hat),
    (n - n_obs + 1):n
  ))
  names(y_hat_t) <- names(df_cover)

  y_low_t <- as.data.frame(cbind(
    t(y_low),
    (n - n_obs + 1):n
  ))
  names(y_low_t) <- names(df_cover)

  y_high_t <- as.data.frame(cbind(
    t(y_high),
    (n - n_obs + 1):n
  ))
  names(y_high_t) <- names(df_cover)

  y_hat_long <- tidyr::pivot_longer(
    y_hat_t,
    cols = s1:s4,
    names_to = "species",
    values_to = "cover_estim"
  )
  y_low_long <- tidyr::pivot_longer(
    y_low_t,
    cols = s1:s4,
    names_to = "species",
    values_to = "cover_low"
  )
  y_high_long <- tidyr::pivot_longer(
    y_high_t,
    cols = s1:s4,
    names_to = "species",
    values_to = "cover_high"
  )

  df_cover_wpred <- left_join(
    df_long[df_long$t %in% (n - n_obs + 1):n, ],
    y_hat_long,
    by = c("t", "species")
  ) %>% left_join(
    .,
    y_low_long,
    by = c("t", "species")
  ) %>% left_join(
    .,
    y_high_long,
    by = c("t", "species")
  )

  ggplot(data = df_cover_wpred)+
    geom_ribbon(aes(x = t, ymin = cover_low, ymax = cover_high, fill = species), alpha = 0.5)+
    geom_line(aes(x = t, y = cover_estim, color = species), linetype = "dashed")+
    geom_line(aes(x = t, y = cover, color = species), size = 0.5)+
    geom_point(aes(x = t, y = cover, color = species), size = 0.5)+
    theme_classic()


# now let's see about the matrix
  df_Amat <- data.frame(
    x = rep(1:S, each = S),
    y = rep(1:S, S),
    z = as.vector(A_mat)
  )

# extract posteriors
  beta_post <- as.data.frame(rstan::extract(mfit, pars = "beta"))
  df_Amat$pred <- apply(beta_post, 2, mean)
  df_Amat$include <- apply(
    beta_post, 2,
    FUN = function(x){
      as.numeric(mean(x < 0) > 0.75 |
        mean(x < 0) < 0.25)
    }
  )

  true_A <- ggplot(data = df_Amat, aes(x = x, y = y))+
    geom_tile(aes(fill = z))

  estim_a <- ggplot(data = df_Amat, aes(x = x, y = y))+
    geom_tile(aes(fill = include))

  gridExtra::grid.arrange(
    true_A, estim_a
  )









