#################################################################
# Testing Ricker models with log-normal demographic stochasticity
#################################################################


# libraries
  library(rstan)
  library(here)
  library(parallel)
  devtools::load_all()


##### Generate some parameter values that will be consistent through the following tests #####

  set.seed(19873)
  steps <- 100                                      # number of time steps
  nsp <- 40                                         # number of species in the community
  lambdas <- runif(nsp, min = 1.2, max = 1.8)       # intrinsic growth rates
  alphas <- runif(nsp, min = 0.005, max = 0.01)     # intraspecific competition
  sigma <- runif(nsp, 0.1, 0.5)                     # demographic stochasticity
  N0 <- rpois(nsp, lambda = 20)                     # initial abundances

  # competition matrix
  A_mat <- comp_matrix2(
    n_sp = nsp,
    rho = 0,
    num_ngs = 2,
    alpha = alphas,
    ng_range = c(0.1, 0.5)
  )

  # model
  growth_mod <- stan_model(here("Stan/pop_growth_rate_FHS.stan"))


##### Simulate full community #####

  N <- matrix(nrow = nsp, ncol = steps)
  N[, 1] <- N0

  for(t in 2:steps){

    N[, t] <- N[, t - 1] * lambdas * exp(- A_mat %*% N[, t - 1] + rnorm(nsp) * sigma / sqrt(N[, t - 1]))

  }

  # # plot the community
  # df_sim <- data.frame(
  #   t = rep(1:50, nsp),
  #   species = as.factor(rep(1:nsp, each = steps)),
  #   abundance = as.vector(t(N))
  # )
  # com_plot <- ggplot(data = df_sim, aes(x = t, y = abundance, color = species)) +
  #   geom_line() +
  #   theme_bw()

  N_het_full <- t(N[-1, ])
  N_foc_full <- as.double(N[1, ])
  tau_0_full <- tau0(
    y = log(N_foc_full[2:steps] / N_foc_full[1:(steps - 1)]),
    m0 = 2,
    M = nsp - 1,
    N = steps - 1,
    fam = "gaussian"
  )

  datlist_full <- list(
    N = steps,
    P = nsp - 1,
    y = N_foc_full,
    X_beta = scale(N_het_full),
    error_scl = 0.5,
    tau0 = tau_0_full,
    slab_scl = 0.5,
    slab_df = 20
  )

  mfit_full <- sampling(
    growth_mod,
    data = datlist_full,
    cores = 3, chains = 3,
    control = list(adapt_delta = 0.99, max_treedepth = 15)
  )

  beta_post_full <- rstan::extract(mfit_full, pars = "beta")$beta %*%
    solve(diag(apply(N_het_full, 2, sd)))
  df_full <- data.frame(
    param = paste0("alpha_1", 2:nsp),
    truth = as.double(A_mat[1, -1]),
    estim = apply(beta_post_full, 2, mean),
    low = apply(beta_post_full, 2, quantile, probs = 0.025),
    high = apply(beta_post_full, 2, quantile, probs = 0.975)
  )

  p_full <- ggplot(data = df_full, aes(x = param))+
    geom_point(aes(y = estim)) +
    geom_errorbar(aes(ymin = low, ymax = high), width = 0) +
    geom_point(aes(y = -truth), color = "brown", shape = 2, size = 4) +
    theme_classic() +
    theme(axis.text.x = element_text(hjust = 1, angle = 45)) +
    ggtitle("Full community simulations") +
    xlab("Parameter") +
    ylab("Value")


##### Simulate with randomly selected heterospecific abundances #####

  # initialize response
  N_foc <- vector(mode = "double", length = steps)
  N_foc[1] <- N0[1]

  # randomly generate heterospecific abundances for each time step
  N_het_ind <- mvtnorm::rmvnorm(
    n = steps,
    mean = apply(N_het_full[20:steps, ], 2, mean),
    sigma = diag(apply(N_het_full[20:steps, ], 2, var))
  )

  # demographic stochasticity of the focal species
  sigma_foc <- sigma[1]

  # create response vector
  for(t in 2:steps){
    N_foc[t] <- N_foc[t - 1] * lambdas[1] * exp(-alphas[1] * N_foc[t - 1] - N_het_ind[t - 1, ] %*% A_mat[1, -1] +
      rnorm(1) * sigma_foc / sqrt(N_foc[t - 1]))
  }

  # compile data list
  tau_0 <- tau0(
    y = log(N_foc[2:steps]/N_foc[1:(steps - 1)]),
    m0 = 2,
    M = ncol(N_het_ind),
    N = steps - 1,
    fam = "gaussian"
  )

  datlist_ind <- list(
    N = steps,
    P = ncol(N_het_ind),
    y = N_foc,
    X_beta = scale(N_het_ind),
    error_scl = 0.5,
    tau0 = tau_0,
    slab_scl = 0.5,
    slab_df = 6
  )

  # fit the model
  mfit_ind <- sampling(
    growth_mod,
    data = datlist_ind,
    cores = 3, chains = 3,
    control = list(adapt_delta = 0.99, max_treedepth = 15)
  )

# plots
  beta_post_ind <- rstan::extract(mfit_ind, pars = "beta")$beta %*%
    solve(diag(apply(N_het_ind, 2, sd)))
  df_ind <- data.frame(
    param = paste0("alpha_1", 2:nsp),
    truth = as.double(A_mat[1, -1]),
    estim = apply(beta_post_ind, 2, mean),
    low = apply(beta_post_ind, 2, quantile, probs = 0.025),
    high = apply(beta_post_ind, 2, quantile, probs = 0.975)
  )

  p_ind <- ggplot(data = df_ind, aes(x = param))+
    geom_point(aes(y = estim)) +
    geom_errorbar(aes(ymin = low, ymax = high), width = 0) +
    geom_point(aes(y = -truth), color = "brown", shape = 2, size = 4) +
    theme_classic() +
    theme(axis.text.x = element_text(hjust = 1, angle = 45)) +
    ggtitle("Independent heterospecific abundances") +
    annotate("point", color = "brown", shape = 2, x = 30, y = -0.001) +
    annotate("text", hjust = 0, label = "True value", x = 31, y = -0.001) +
    annotate("point", x = 30, y = -0.0015) +
    annotate("text", hjust = 0, label = "Estimate", x = 31, y = -0.0015) +
    xlab("") +
    ylab("Value")




##### Is this due simply to collinearity? #####

  # initialize another round of responses
  N_foc_cor <- vector(mode = "double", length = steps)
  N_foc_cor[1] <- N0[1]

  # draw values from multivariate normal distribution with empirical covariance matrix
  N_het_cor <- mvtnorm::rmvnorm(
    n = steps,
    mean = apply(N_het_full[20:steps, ], 2, mean),
    sigma = cov(N_het_full[20:steps, ])
  )

  # create response vector
  for(t in 2:steps){
    N_foc_cor[t] <- N_foc_cor[t - 1] * lambdas[1] * exp(-alphas[1] * N_foc_cor[t - 1] - N_het_cor[t - 1, ] %*% A_mat[1, -1]) +
      rnorm(1) * sigma_foc / sqrt(N_foc_cor[t - 1])
  }

  # compile data list
  tau_0_cor <- tau0(
    y = log(N_foc_cor[2:steps]/N_foc_cor[1:(steps - 1)]),
    m0 = 2,
    M = ncol(N_het_cor),
    N = steps - 1,
    fam = "gaussian"
  )

  datlist_cor <- list(
    N = steps,
    P = ncol(N_het_cor),
    y = N_foc_cor,
    X_beta = scale(N_het_cor),
    error_scl = 0.5,
    tau0 = tau_0_cor,
    slab_scl = 0.5,
    slab_df = 6
  )

  # fit the model
  mfit_cor <- sampling(
    growth_mod,
    data = datlist_cor,
    cores = 3, chains = 3,
    control = list(adapt_delta = 0.99, max_treedepth = 15)
  )


  # get the plot
  beta_post_cor <- rstan::extract(mfit_cor, pars = "beta")$beta %*%
    solve(diag(apply(N_het_cor, 2, sd)))
  df_cor <- data.frame(
    param = paste0("alpha_1", 2:nsp),
    truth = as.double(A_mat[1, -1]),
    estim = apply(beta_post_cor, 2, mean),
    low = apply(beta_post_cor, 2, quantile, probs = 0.025),
    high = apply(beta_post_cor, 2, quantile, probs = 0.975)
  )

  p_cor <- ggplot(data = df_cor, aes(x = param))+
    geom_point(aes(y = estim)) +
    geom_errorbar(aes(ymin = low, ymax = high), width = 0) +
    geom_point(aes(y = -truth), color = "brown", shape = 2, size = 4) +
    theme_classic() +
    theme(axis.text.x = element_text(hjust = 1, angle = 45)) +
    ggtitle("Correlated abundances") +
    xlab("") +
    ylab("Value")


# save the plots
  png(
    filename = here("Figures/param_checks_terr_comm_eef.png"),
    height = 3600, width = 1500,
    units = "px", res = 300
  )
    gridExtra::grid.arrange(p_ind, p_cor, p_full, ncol = 1)
  dev.off()







