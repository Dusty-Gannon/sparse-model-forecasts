#########################################################
# Fitting Latent AR(1) models to simulated community data
#########################################################


# libraries
  library(rstan)
  library(here)
  library(tidyr)
  library(dplyr)
  devtools::load_all()

# load data
  dat <- readRDS(
    here("Data/terrestrial_sim_data/simdat_500reps_S2s48_5ann.rds")
  )[[1]]

# convert cover data to wide format
  cover_wide <- pivot_wider(
    dat$cover,
    names_from = species,
    values_from = cover
  )

# store some handy variables
  n <- nrow(cover_wide)

  # annuals
  ann <- as.character(1:5)

# find annual with greatest abundance at the end that also has at least
#  one strong competitor
  col_ids_a <- which(names(cover_wide) %in% ann)
  foc_a <- choose_focal(
    df = cover_wide,
    col_ids = col_ids_a,
    t = n, num_ngs = 2
  )

# find a perennial with high abundance and at least one
#  strong competitor
  col_ids_p <- which(!(names(cover_wide) %in% c(ann, "t")))
  foc_p <- choose_focal(
    df = cover_wide,
    col_ids = col_ids_p,
    t = n, num_ngs = 2
  )

  tsteps_p1 <- (n - 100 + 1):n
  tsteps_m1 <- (n - 100):(n - 1)
  tau_0 <- tau0(
    y = pull(cover_wide, foc_a)[tsteps_p1],
    m0 = 4,
    M = ncol(cover_wide) - 2,
    N = 100,
    fam = "gamma"
  )
# with the last 100 time steps, try fitting the Lat AR1 model
  datlist_a100 <- list(
    N = 100,
    P0 = 1,
    P = ncol(cover_wide) - 2,
    y = pull(cover_wide[tsteps_p1, ], foc_a),
    X_alpha = matrix(data = 1, nrow = 100, ncol = 1),
    X_beta = as.matrix(
      cover_wide[tsteps_m1, -which(colnames(cover_wide) %in% c("t", foc_a))]
    ),
    tau0 = tau_0,
    slab_scl = 1,
    slab_df = 10
  )

  gamma_latAR1 <- stan_model(here("Stan/Gamma_LatAR1_FHS.stan"))

  mfit_a100 <- sampling(
    gamma_latAR1,
    data = datlist_a100,
    cores = 4,
    chains = 4,
    control = list(adapt_delta=0.99, max_treedepth = 12)
  )

# save the fit
  saveRDS(
    mfit_a100,
    file = here("Data/terrestrial_sim_data/mfit_sim1_a100.rds")
  )











