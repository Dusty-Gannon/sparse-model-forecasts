# document setup
library(here)
library(tidyverse)
devtools::load_all()


##############################################################
# This script simulates competitive communities using
# a Ricker population model for each competing species.
# Demographic stochasticity enters in the exponent,
# resulting in log-normally distributed errors. For a single
# set of parameters, we simulate data from five scenarios:
# 1) Full community dynamics
# 2) Same parameters, but random heterospecific abundances
#    with the same means and variances as the full community
#    but zero covariances.
# 2) Same parameters, but random and correlated heterospecific
#    abundances.
##############################################################


###### Function to simulate competitive communities #####
simulate_communities <- function(vr, sim_params){

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

  # compute means and variances on abundance scale
  mus <- apply(N_full_final[-1, 50:steps], 1, mean)
  s_N <- apply(N_full_final[-1, 50:steps], 1, sd)
  D <- diag(sqrt(vr) * s_N)
  S <- D %*% cor(t(N_full_final[-1, 50:steps])) %*% D

  N_cor <- with(c(sim_params, list(foc = foc, ext = ext, mus = mus, S = S, vr = vr)), {

    # draw from mvn distribution
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
        exp(-A_mat[foc, foc] * N_foc[t - 1] - A_mat[foc, -c(foc, ext)] %*% N_het[, t - 1] + rnorm(1) * (sigmas[foc] / sqrt(vr)) / sqrt(N_foc[t - 1]))
    }

    rbind(N_foc, N_het)

  })

  return(list(
    N_full = N_full_final,
    N_cor = N_cor,
    sim_params = list(
      lambdas = sim_params$lambdas[nesp],
      A_mat = sim_params$A_mat[nesp, nesp],
      sigmas = sim_params$sigmas[nesp]
    )
  ))

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


##### Defining parameter values #####
  set.seed(789)

  # some useful global variables
  nsp = 60; num_ngs = 5; steps = 200; iter = 1000

  params <- list(
    nsp = nsp, steps = steps, num_ngs = num_ngs,
    sigmas = rep(sqrt(quantile(jr_demstoch, probs = 0.5)) / 2, 60),
    A_mat = comp_matrix2(
      n_sp = nsp, rho = 0, alpha = runif(nsp, 0.01, 0.05),
      num_ngs = num_ngs, ng_range = c(0.1, 0.5)
    ),
    lambdas = runif(nsp, 1.2, 1.8),
    N0 = rpois(nsp, lambda = 50)
  )

# create variance amplification factors
  het_vr <- rep(c(0.5, 1, 2, 5, 20), each = iter / 5)

##### Running the simulations #####

  sims <- lapply(
    het_vr,
    FUN = simulate_communities,
    sim_params = params
  )

# ##### Quick checks #####
 # sapply(sims, function(x){sum(x$N_full[, steps] > 0)})
 #
 #  plot_comm <- function(N){
 #    library(ggplot2)
 #    nsp <- nrow(N)
 #    steps <- ncol(N)
 #
 #    df <- data.frame(
 #      t = rep(1:steps, each = nsp),
 #      sp = as.factor(rep(1:nsp, steps)),
 #      N = as.vector(N)
 #    )
 #
 #    return(
 #      ggplot(data = df, aes(x = t, y = N, color = sp)) +
 #        geom_line() +
 #        theme_classic() +
 #        ylim(c(0, 50))
 #    )
 #  }
 #
 #  plot_comm(sims_final[[450]]$N_full)

##### Save the simulated datasets #####

  # find and remove any cases where the focal went extinct
  to_keep <- which(
    {
      sapply(sims, function(x){sum(is.na(x$N_full))}) +
      sapply(sims, function(x){sum(is.na(x$N_cor))})
    } == 0
  )

  # sample 100 sims from each het_vr scenario
  ids <- vector(mode = "list", length = length(unique(het_vr)))
  for(i in 1:length(unique(het_vr))){
    ids[[i]] <- ((i - 1) * (iter / 5) + 1):(i * (iter / 5))
  }

  sims_final <- vector(mode = "list")
  for(i in 1:length(unique(het_vr))){
    sims_final <- c(
      sims_final,
      sims[sample(to_keep[to_keep %in% ids[[i]]], 100)]
    )
  }

# save simulated data
  fname <- paste0(
    "lnorm_ricker_sims_",
     steps, "steps_",
    "rho0", "_S", num_ngs,
    "_s", nsp - num_ngs, "_1param.rds"
  )

  fp <- paste0(
    "Data/terrestrial_sim_data/lnorm_ricker/",
    fname
  )

  saveRDS(sims_final, file = here(fp))

