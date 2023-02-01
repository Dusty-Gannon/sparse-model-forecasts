# document setup
library(here)
library(tidyverse)
devtools::load_all()


##############################################################
# This script simulates competitive communities using
# a Ricker population model for each competing species.
# Demographic stochasticity enters in the exponent,
# resulting in log-normally distributed errors. For each
# set of parameters, we simulate data from three scenarios:
# 1) Full community dynamics
# 2) Same parameters, but random heterospecific abundances
#    with the same means and variances as the full community
#    but zero covariances.
# 2) Same parameters, but random and correlated heterospecific
#    abundances.
##############################################################


###### Function to simulate competitive communities #####
simulate_communities <- function(sim_params){

  # generate full communities
  N_full <- with(sim_params, {

    # initialize tracking matrix
    N <- matrix(0, nrow = nsp, ncol = steps)
    N[, 1] <- N0

    for(t in 2:steps){

      # tracking extinct species
      extinct <- which(round(as.double(N[, t - 1])) == 0)
      N[extinct, t - 1] <- 0
      nesp <- (1:nsp)[-extinct]
      if(length(extinct) > 0){
        N[nesp, t] <- N[nesp, t - 1] * lambdas[nesp] *
          exp(- A_mat[nesp, nesp] %*% N[nesp, t - 1] + rnorm(length(nesp)) * sigmas[nesp] / sqrt(N[nesp, t - 1]))
      } else{
        N[, t] <- N[, t - 1] * lambdas * exp(- A_mat %*% N[, t - 1] + rnorm(nsp) * sigmas / sqrt(N[, t - 1]))
      }

    }

    N

  })

  # simulate another round but with random draws for the
  #  heterospecifics that amplify the ratio of variation in
  #  heterospecifics relative to demographic stochasticity

  # get species that are not extinct at the end
  nesp <- which(!(N_full[, sim_params$steps] == 0))
  ext <- which(N_full[, sim_params$steps] == 0)
  N_full_final <- N_full[nesp, ]
  foc <- nesp[1]

  # compute means and variancs on abundance scale
  mus <- apply(N_full_final[-1, 50:steps], 1, mean)
  s_N <- apply(N_full_final[-1, 50:steps], 1, sd)
  D <- diag(sqrt(sim_params$het_vr) * s_N)
  S <- D %*% cor(t(N_full_final[-1, 50:steps])) %*% D

  N_cor <- with(c(sim_params, list(foc = foc, ext = ext, mus = mus, S = S)), {

    # draw from distribution, then exponentiate
    N_het <- t(mvtnorm::rmvnorm(
      steps,
      mean = mus,
      sigma = S
    ))
    N_het[N_het < 0] <- 0

    # initialize focal species vector
    N_foc <- vector(mode = "double", length = steps)
    foc <- nesp[1]
    N_foc[1] <- N0[foc]

    # simulate ricker model
    for(t in 2:steps){
      N_foc[t] <- N_foc[t - 1] * lambdas[foc] *
        exp(-A_mat[foc, foc] * N_foc[t - 1] - A_mat[foc, -c(foc, ext)] %*% N_het[, t - 1] + rnorm(1) * (sigmas[foc] / sqrt(het_vr)) / sqrt(N_foc[t - 1]))
    }

    rbind(N_foc, N_het)

  })

  # simulate another round, but with equal variances across heterospecifics

  N_eqv <- with(c(sim_params, list(foc = foc, ext = ext, mus = mus, S = S)), {

    # create heterospecific densities through time
    m <- rep(quantile(mus, probs = 0.5), length(mus))
    s <- quantile(diag(S), probs = 0.5)
    N_het <- t(mvtnorm::rmvnorm(
      n = steps,
      mean = m,
      sigma = diag(s, nrow = length(mus), ncol = length(mus))
    ))
    N_het[N_het < 0] <- 0

    # initialize focal species vector
    N_foc <- vector(mode = "double", length = steps)
    N_foc[1] <- N0[foc]

    # simulate ricker model
    for(t in 2:steps){
      N_foc[t] <- N_foc[t - 1] * lambdas[foc] * exp(-A_mat[foc, foc] * N_foc[t - 1] - A_mat[foc, -c(foc, ext)] %*% N_het[, t - 1] + rnorm(1) * (sigmas[foc] / sqrt(het_vr)) / sqrt(N_foc[t - 1]))
    }

    rbind(N_foc, N_het)

  })

  return(list(
    N_full = N_full_final,
    N_cor = N_cor,
    N_eqv = N_eqv,
    sim_params = list(
      lambdas = sim_params$lambdas[nesp],
      A_mat = sim_params$A_mat[nesp, nesp],
      sigmas = sim_params$sigmas[nesp]
    )
  ))

}


##### Function to generate the simulation parameters #####
generate_sim_params <- function(
  x, nsp = 40, steps = 200, num_ngs = 3, sigma_rng = c(0.1, 0.5),
  alpha_rng = c(0.005, 0.01), lambda_rng = c(1.2, 1.8),
  ng_range = c(0.2, 0.4), rho = 0, mean_init_abund = 20,
  comp_matrix_type = 2, het_vr = 1
){

  # generate competition matrix
  alpha <- runif(nsp, min = alpha_rng[1], max = alpha_rng[2])
  if(comp_matrix_type == 1){
    A_mat <- sponges::comp_matrix(
      n_sp = nsp, rho = rho, alpha = alpha,
      num_ngs = num_ngs, ng_range = ng_range
    )
  }
  if(comp_matrix_type == 2){
    A_mat <- sponges::comp_matrix2(
      n_sp = nsp, rho = rho, alpha = alpha,
      num_ngs = num_ngs, ng_range = ng_range
    )
  }

  # compile list of return objects
  return(
    list(
      nsp = nsp,
      steps = steps,
      sigmas = runif(nsp, min = sigma_rng[1], max = sigma_rng[2]),
      A_mat = A_mat,
      lambdas = runif(nsp, min = lambda_rng[1], max = lambda_rng[2]),
      N0 = rpois(nsp, lambda = mean_init_abund),
      het_vr = het_vr
    )
  )

}


##### Get an idea of variation from empirical data to use in sims #####

  jr_data <- read.csv(here("Data/JR_Data/JR_count.csv"), header = T)
  jr_p1 <- subset(jr_data, Site == "83A") %>%
    group_by(Year, Species) %>% summarise(count = sum(Data))

  dem_var_estims <- function(sp, df){

    y <- subset(df, Species == sp)$count
    n <- length(y)
    # if(sum(y <= 0) > 0){
    #   return(NA)
    # }

    y2 <- y + 1

    r <- log(y2[2:n] / y2[1:(n - 1)])
    ricker_mod <- lm(r ~ 1 + y2[1:(n - 1)])

    return(
      sum(ricker_mod$residuals^2) / (n - 2)
    )

  }


  jr_demstoch <- sapply(unique(jr_p1$Species), dem_var_estims, df = jr_p1)
  min(jr_demstoch); max(jr_demstoch)


##### Defining ranges for parameter generation #####

  nsp <- 60; steps <- 200; num_ngs <- 3;
  sigma_rng <- c(0.1, sqrt(max(jr_demstoch)) / 2);
  alpha_rng <- c(0.01, 0.05); lambda_rng <- c(1.2, 1.8);
  ng_range <- c(0.1, 0.5); rho <- 0; mean_init_abund <- 50

##### Running the simulations #####

  # number of iterations
  iter = 500

  # generate list of parameters
  params <- lapply(
    1:iter,
    FUN = generate_sim_params,
    nsp = nsp, steps = steps,
    num_ngs = num_ngs, sigma_rng = sigma_rng,
    alpha_rng = alpha_rng, lambda_rng = lambda_rng,
    ng_range = ng_range, rho = rho,
    mean_init_abund = mean_init_abund,
    comp_matrix_type = 2
  )

  # change the variance factor across the sims
  het_vr <- rep(c(0.5, 1, 2, 5, 10), each = 100)
  for(i in 1:iter){
    params[[i]]$het_vr <- het_vr[i]
  }

  sims <- lapply(
    params,
    FUN = simulate_communities
  )

# ##### Quick checks #####
#  sapply(sims, function(x){sum(x$N_full[, steps] > 0)})
#
#   plot_comm <- function(N){
#     library(ggplot2)
#     nsp <- nrow(N)
#     steps <- ncol(N)
#
#     df <- data.frame(
#       t = rep(1:steps, each = nsp),
#       sp = as.factor(rep(1:nsp, steps)),
#       N = as.vector(N)
#     )
#
#     return(
#       ggplot(data = df, aes(x = t, y = N, color = sp)) +
#         geom_line() +
#         theme_classic()
#     )
#   }
#
#   plot_comm(sims[[401]]$N_cor)

##### Save the simulated datasets #####

  fname <- paste0(
    "lnorm_ricker_sims_",
     steps, "steps_",
    "rho", rho, "_S", num_ngs,
    "_s", nsp - num_ngs, ".rds"
  )

  fp <- paste0(
    "Data/terrestrial_sim_data/lnorm_ricker/",
    fname
  )

  saveRDS(sims, file = here(fp))

