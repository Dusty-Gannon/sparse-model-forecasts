#########################################################
# Fitting Latent AR(1) models to simulated community data
#########################################################


# libraries
  library(rstan)
  library(here)
  library(dplyr)
#  devtools::load_all()

# load user-defined functions (teton way, not devtools way)
  src_files <- list.files(here("R/"), pattern = "*.R", full.names = T)
  sapply(src_files, source, .GlobalEnv)

# load data
  dat <- readRDS(
    here("Data/terrestrial_sim_data/simdat_500reps_S2s48_5ann.rds")
  )[[1]]

# number of time points to use as data
  n_obs <- 100

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

# store some handy variables
  n <- nrow(cover_df)

  # annuals
  ann <- as.character(1:5)

# find annual with greatest abundance at the end that also has at least
#  one strong competitor
  col_ids_a <- which(names(cover_df) %in% ann)
  foc_a <- choose_focal(
    df = cover_df,
    col_ids = col_ids_a,
    t = n, num_ngs = 2
  )

# find a perennial with high abundance and at least one
#  strong competitor
  col_ids_p <- which(!(names(cover_df) %in% c(ann, "t")))
  foc_p <- choose_focal(
    df = cover_df,
    col_ids = col_ids_p,
    t = n, num_ngs = 2
  )

  tsteps_p1 <- (n - n_obs + 1):n
  tsteps_m1 <- (n - n_obs):(n - 1)
  tau_0 <- tau0(
    y = cover_df[tsteps_p1, foc_a],
    m0 = 4,
    M = ncol(cover_df) - 2,
    N = n_obs,
    fam = "gaussian"
  )
# with the last 200 time steps, try fitting the Lat AR1 model
# try with all variables centered and scaled
  X_beta <- as.matrix(scale(
    cover_df[tsteps_m1, -which(colnames(cover_df) %in% c("t", foc_a))]
  ))
  datlist_a100 <- list(
    N = n_obs,
    P0 = 1,
    P = ncol(cover_df) - 2,
    p = 2,
    y = as.double(scale(cover_df[tsteps_p1, foc_a])),
    X_alpha = matrix(data = 1, nrow = n_obs, ncol = 1),
    X_beta = X_beta,
    tau0 = tau_0,
    slab_scl = 1,
    slab_df = 10
  )

  # some species may be extinct towards the end of the time series
  extincts <- which(
    apply(datlist_a100$X_beta, 2, function(x){
      sum(is.nan(x))
    }) > 0
  )
  datlist_a100$X_beta <- datlist_a100$X_beta[, -extincts]
  datlist_a100$P <- ncol(datlist_a100$X_beta)

  arp_FHS <- stan_model(here("Stan/AR-p_FHS.stan"))

  mfit_a100 <- sampling(
    arp_FHS,
    data = datlist_a100,
    cores = 3,
    chains = 3,
    control = list(adapt_delta=0.95, max_treedepth = 12)
  )

  plot(datlist_a100$y[2:100], pch = 20)
  lines(datlist_a100$y[2:100])
  y_rep <- rstan::extract(mfit_a100, pars = "y_rep")$y_rep
  lines(apply(y_rep, 2, mean), col = "blue")
  lines(apply(y_rep, 2, quantile, probs = 0.05), col = "blue", lty="dashed")
  lines(apply(y_rep, 2, quantile, probs = 0.95), col = "blue", lty="dashed")


  # check the residuals
  resids <- rstan::extract(mfit_a100, pars = "resid")$resid

  df_resid <- data.frame(
    t = 3:100,
    mean = apply(resids, 2, mean),
    low = apply(resids, 2, quantile, probs = 0.05),
    high = apply(resids, 2, quantile, probs = 0.95)
  )

  ggplot(df_resid, aes(x = t, y = mean))+
    geom_errorbar(aes(ymin = low, ymax = high), width = 0)+
    geom_point()+
    theme_bw()

  plot(qresids)
  plot(seq(-3, 3, length.out = 98), y = sort(qresids))
  abline(a = 0, b = 1)

  # check for a mean-variance relationship
  y <- cover_df[tsteps_p1, foc_a]
  thresh <- seq(min(y) - 0.1, max(y) + 0.1, length.out = 8)
  df_bins <- data.frame(
    y = y,
    bin = 1
  )
  for(i in 1:nrow(df_bins)){
    b <- 1
    while(df_bins$y[i] > thresh[b]){
      b <- b + 1
    }
    df_bins$bin[i] <- b - 1
  }

  df_mv <- group_by(df_bins, bin) %>%
    summarise(mean = mean(y), var = var(y))

# save the fit
  saveRDS(
    mfit_a100,
    file = here("Data/terrestrial_sim_data/mfit_sim1_a100.rds")
  )






