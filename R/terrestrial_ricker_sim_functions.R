

#' Generate an array of parameters that determine competition
#'
#' This function generates an array of parameters that, when combined with
#' a model matrix of environmental variables, determines the matrix of
#' competition coefficients.
#'
#' @param nsp Number of species in the community
#' @param num_ngs Number of non-generic competitors (usually, stronger competitors)
#' @param num_env Number of environmental variables that help to determine
#' the linear effect of the variables on the competitive effects of some select dynamic
#' species whose competitive abilities vary with the environment. The effect is drawn from
#' a uniform distribution on the interval \eqn{(-\alpha_{ii}, \alpha_{ii})},
#' where \eqn{\alpha_{ii}} is the intra-specific effect of the species for which the
#' dynamic species serves as a non-generic competitor.
#' \eqn{\beta_0} is the generic effect. This helps to keep inter-specific effects weaker than
#' intra-specific effects.
#' @param num_dynamic Number of species whose competitive abilities vary with the environment
#' @param generic Effect of adding an individual of a generic competitor to the growth of the
#' focal species
#' @param ng_strength Factor to multiply the generic effect by for the deviation from the
#'  generic strength to the non-generic competitor's effect
#' @param intra_strength Factor to multiply the generic effect by for the deviation from the
#' generic strength to the intraspecific effect
#'
#' @return An array with \eqn{S} slices, each with dimensions \eqn{S \times V + 2}, where \eqn{S} is the
#' number of species in the community and \eqn{V} is the number of environmental variables.
#' @export
#'
#' @examples
#'
#' ricker_comp_matrix(
#'   nsp = 10, num_ngs = 3,
#'   env_sd = c(0.5, 0.1), num_dynamic = 2,
#'   generic = 0.001
#' )
#'
ricker_comp_array <- function(
    nsp, num_ngs, num_env,
    num_dynamic, generic = 0.001,
    ng_strength = c(20, 50), intra_strength = c(80, 100)
){

  # randomly sample the species that get more or less competitive
  #  depending on the environment
  dyn_sp <- sample(1:nsp, size = num_dynamic)

  # initialize array
  B <- array(data = 0, dim = c(nsp, num_env + 2, nsp))

  # fill in the first column of each slice with the generic effect
  B[, 1, ] <- generic

  # proceed to fill in the remaining elements of each matrix
  for(j in 1:nsp){

    # intraspecific effect
    B[j, 2, j] <- generic * runif(1, min = intra_strength[1], max = intra_strength[2])

    # non generic effects
    if(nsp - j >= num_ngs){
      B[j + c(1:num_ngs), 2, j] <- generic * runif(num_ngs, min = ng_strength[1], max = ng_strength[2])
      for(i in j + c(1:num_ngs)){
        if(i %in% dyn_sp){
          B[i, (1:num_env) + 2, j] <- runif(
            num_env,
            min = - B[j, 2, j],
            max = B[j, 2, j]
          )
        }
      }
    }

    # begin wrapping once we start to reach the edge
    if(nsp - j < num_ngs & nsp - j > 0){
      # get ng_sp ids
      ng_ids <- c(
        1:(num_ngs - (nsp - j)),
        j + c(1:(nsp - j))
      )
      B[ng_ids, 2, j] <- generic * runif(num_ngs, min = 20, max = 50)
      for(i in ng_ids){
        if(i %in% dyn_sp){
          B[i, (1:num_env) + 2, j] <- runif(
            num_env,
            min = -B[j, 2, j],
            max = B[j, 2, j]
          )
        }
      }
    }

    if(nsp - j == 0){
      B[1:num_ngs, 2, j] <- generic * runif(num_ngs, min = 20, max = 50)
      for(i in (1:num_ngs)){
        if(i %in% dyn_sp){
          B[i, (1:num_env) + 2, j] <- runif(
            num_env,
            min = -B[j, 2, j],
            max = B[j, 2, j]
          )
        }
      }
    }

  }

  return(list(
    B = B,
    dynamic_spids = dyn_sp
  ))

}






#' Gaussian function to determine population growth given the environment
#'
#' This function uses optimal environmental growing conditions for a given species
#' and the current environment at time \eqn{t} to determine intrinsic growth at time \eqn{t}.
#'
#' @param env_t Environment at time t (Numeric).
#' @param lambda_max Vector of max growth rates for each species, achieve when \code{env_t} is
#' equal to a species' optimum.
#' @param optims Vector of optimal growing conditions for each species.
#' @param tau Scalar or vector of sensitivity parameters that determine how quickly a species
#' intrinsic growth rate deviates from \code{lambda_max} with changes in environmental conditions.
#'
#' @return Vector of intrinsic growth rates given the environment.
#'
#' @examples
#' env <- rnorm(1)
#' lambda_max <- runif(10, min = 1, max = 2)
#' optims <- runif(10, min = -0.5, max = 0.5)
#' gauss_env_effect(env, lambda_max, optims)
#'
#'
gauss_env_effect <- function(env, lambda_max, optims, tau = 1){

  lambda_max * exp(-tau * (env - optims)^2)

}







#' Ricker competition model with Poisson-distributed demographic stochasticity
#'
#' This function computes the abundance of \eqn{S} species at time
#' \eqn{t + 1} given their abundances at time \eqn{t} based on a
#' Ricker population growth model
#'
#' @param N_t Vector of population abundances at time \eqn{t}
#' @param lambdas Vector of intrinsic growth rates
#' @param A_mat Matrix of competition coefficients
#' @param stochastic Logical determining whether there is demographic stochasticity or not
#'
#' @return A vector of species abundances at time \eqn{t + 1}
#' @export
#'
#' @examples
#' N_t <- rpois(5, lambda = 10)
#' lambdas <- runif(5)
#' A_mat <- comp_matrix(
#'   n_sp = 5, rho = c(0.1, 0.1),
#'   alpha = runif(5, min = 0.05, max = 0.1),
#'   num_ngs = 0, num_regs = 0
#' )
#' ricker_step(N_t, lambdas, A_mat)
#'
ricker_step_pois <- function(N_t, lambdas, A_mat, stochastic = T){

  # get number of species in the community
  S <- length(N_t)

  # define abundance vector for time t plus 1
  N_tp1 <- vector(mode = "double", length = S)

  # define the next step
  for(s in 1:S){
    N_tp1[s] <- N_t[s] * lambdas[s] * exp(
      - A_mat[s, ] %*% N_t
    )
  }

  # return either a stochastic or deterministic abundance for the next step
  if(isTRUE(stochastic)){
    return(rpois(S, lambda = N_tp1))
  } else{
    return(N_tp1)
  }

}









#' Generate a list of simulation parameters for lognormal Ricker model
#'
#' @param x Dummy variable so that this function is straight forward to use with apply statements
#' @param nsp Number of species to start with
#' @param steps Number of steps over which to simulate
#' @param num_ngs Number of non-generic species
#' @param sigma_rng Interval over which to draw scales of demographic stochasticity
#' @param alpha_rng Interval from which to draw intraspecific competition coefficients
#' @param lambda_rng Interval from which to draw low-density growth rates
#' @param ng_range Range over which to draw the fraction of intra-specific growth for the non-generics
#' (i.e., 0.5 means that non-generic competition coefficient would be half as large as intra-specific
#' competition)
#' @param rho Fraction of intraspecific competition that is the generic competitive effect
#' @param mean_init_abund Average initial abundance for the random initial states
#' @param comp_matrix_type Non-generic competition can either be random (\code{comp_matrix_type = 1})
#' or structured so that no species gets overloaded with strong competitors (\code{comp_matrix_type = 2}).
#' @param het_vr Multiplier by which to amplify the variance of heterospecific abundances.
#'
#' @return List of parameters for the simulation
#'
generate_sim_params_vrtests <- function(
    x, nsp = 40, steps = 200, num_ngs = 3,
    sigma_rng = c(0.1, 0.5), alpha_rng = c(0.005, 0.01),
    lambda_rng = c(1.2, 1.8), ng_range = c(0.2, 0.4), rho = 0,
    mean_init_abund = 20, comp_matrix_type = 2, het_vr = 1
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









#' Generate simulation parameters for single TS thinning experiments
#'
#' @param x Dummy variable so that this function is straight forward to use with apply statements
#' @param nsp Number of species to start with
#' @param init_steps Number of steps before thinning
#' @param tot_steps Total number of steps over which to simulate
#' @param num_ngs Number of non-generic species
#' @param sigma_rng Interval over which to draw scales of demographic stochasticity
#' @param alpha_rng Interval from which to draw intraspecific competition coefficients
#' @param lambda_rng Interval from which to draw low-density growth rates
#' @param ng_range Range over which to draw the fraction of intra-specific growth for the non-generics
#' (i.e., 0.5 means that non-generic competition coefficient would be half as large as intra-specific
#' competition)
#' @param rho Fraction of intraspecific competition that is the generic competitive effect
#' @param mean_init_abund Average initial abundance for the random initial states
#' @param comp_matrix_type Non-generic competition can either be random (\code{comp_matrix_type = 1})
#' or structured so that no species gets overloaded with strong competitors (\code{comp_matrix_type = 2}).
#' @param thin_freq Frequency of thinning treatments. \code{thin_freq = 0} will proceed without
#' thinning treatments, while a value of \code{thin_freq = c} with c > 0, will perform a
#' thinning treatment every c years.
#' @param prop_cthin Proportion of the community that should be thinned throughout the time series.
#' @param thin_factor Factor by which to thin a species. \code{thin_freq = 0.1} with thin a
#' species to 10 percent of its population density in the previous year.
#'
#' @return List of simulation parameters
#'
generate_sim_params_thin <- function(
    nsp = 40, init_steps = 100, tot_steps = 500, num_ngs = 3,
    sigma_rng = c(0.1, 0.5), alpha_rng = c(0.005, 0.01),
    lambda_rng = c(1.2, 1.8), ng_range = c(0.2, 0.4), rho = 0,
    mean_init_abund = 20, comp_matrix_type = 2, thin_freq = 2,
    prop_cthin = 1, thin_factor = 0.1, target_thin = TRUE
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
      init_steps = init_steps, tot_steps = tot_steps,
      sigmas = runif(nsp, min = sigma_rng[1], max = sigma_rng[2]),
      A_mat = A_mat,
      lambdas = runif(nsp, min = lambda_rng[1], max = lambda_rng[2]),
      N_0 = rpois(nsp, lambda = mean_init_abund),
      thin_factor = thin_factor, thin_freq = thin_freq,
      prop_cthin = prop_cthin, target_thin = target_thin
    )
  )

}






generate_sim_params_dist <- function(
    nsp = 40, steps = 500, num_ngs = 5,
    sigma_rng = c(0.1, 0.5), alpha_rng = c(0.005, 0.01),
    lambda_rng = c(1.2, 1.8), ng_range = c(0.2, 0.4), rho = 0,
    mean_init_abund = 20, comp_matrix_type = 2, dist_prob = 0,
    dist_int = 0, prop_cdist = 0
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
      nsp = nsp, steps = steps,
      sigmas = runif(nsp, min = sigma_rng[1], max = sigma_rng[2]),
      A_mat = A_mat,
      lambdas = runif(nsp, min = lambda_rng[1], max = lambda_rng[2]),
      N_0 = rpois(nsp, lambda = mean_init_abund),
      dist_prob = dist_prob, dist_int = dist_int,
      prop_cdist = prop_cdist
    )
  )

}





#' Simulating a Ricker population model with log-normal demographic stochasticity
#'
#' @param N_0 Initial abundance of the focal species
#' @param lambda Intrinsic growth rate of the focal species
#' @param A_i Vector of competition coefficients with intra-specific competition in the
#' first index
#' @param sigma_i Scale of demographic stochasticity for the focal species
#' @param N_het Matrix of heterospecific abundances through time
#' @param steps Number of steps to simulate
#' @param het_vr Scalar multiplier by which to reduce demographic stochasticity
#'
#' @return Vector of focal species abundances through time
#'
#'
ricker_ts_lnorm_foc <- function(N_0, lambda, A_i, sigma_i, N_het, steps = 300, het_vr = 1){

  N <- vector(mode = "double", length = steps)
  N[1] <- N_0
  for(t in 1:(steps - 1)){
    N[t + 1] <- N[t] * lambda * exp(-A_i[1] * N[t] - N_het[t, ] %*% A_i[-1] + rnorm(1) * (sigma_i / sqrt(het_vr)) / sqrt(N[t]))
  }

  return(N)

}







#' Full community Ricker simulations
#'
#' @param N_0 Vector of inital population abundances
#' @param lambdas Vector of low-density growth rates
#' @param A_mat Competition matrix
#' @param sigmas Vector of scales of demographic stochasticity
#' @param steps Number of time steps
#' @param dist_prob Probability the community is disturbed in a given year.
#' @param dist_int Disturbance intensity, on the unit interval. An intensity of 0.5 translates
#' to a disturbance that, on average, reduces 50% of the species to 50% of their density.
#' @param thin_freq Frequency of thinning treatments. \code{thin_freq = 0} will proceed without
#' thinning treatment, while a value of \code{thin_freq = c} with c > 0, will perform a
#' thinning treatment every c years.
#' @param thin_factor Factor by which to thin a species. \code{thin_freq = 0.1} with thin a
#' species to 10 percent of its population density in the previous year.
#' @param thin_order Order in which to thin species. If \code{NULL}, the order will be
#' generated at random
#'
#' @return Matrix with as many rows as there are species in the simulation and as many columns
#' as steps. If the communities had natural disturbance, a list of which species were disturbed
#' and when is also returned.
#'
#'
ricker_ts_lnorm <- function(
    N_0, lambdas, A_mat, sigmas, steps,
    dist_prob = 0, dist_int = 0, prop_cdist = 0,
    dist_min_thresh = 1, thin_freq = 0, thin_factor = 0.1,
    thin_order = NULL, thin_levels = NULL
){

  nsp <- length(N_0)
  # initialize tracking matrix
  N <- matrix(0, nrow = nsp, ncol = steps)
  N[, 1] <- N_0

  ### No thinning or disturbance ###
  if(thin_freq == 0 & dist_prob == 0){
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
  }

  ### Thinning but no disturbance ###
  if(thin_freq > 0 & dist_prob == 0){

    # create vector of thinning factors
    thin <- rep(1, steps)
    if(is.null(thin_levels)){
      thin[thin_freq * c(1:floor(steps/thin_freq))] <- thin_factor
    } else{
      thin[thin_freq * c(1:floor(steps/thin_freq))] <- 0
    }


    # create vector that cycles through the thinning order
    if(is.null(thin_order)){thin_order <- sample(1:nsp)}
    sp2thin <- rep(
      rep(thin_order, each = thin_freq),
      ceiling(steps / (thin_freq * length(thin_order)))
    )

    # now proceed with simulations
    for(t in 2:steps){

      # tracking extinct species
      extinct <- which(round(as.double(N[, t - 1])) == 0)
      N[extinct, t - 1] <- 0

      # thin species if it was a thinning year
      if(is.null(thin_levels)){
        N[sp2thin[t - 1], t - 1] <- N[sp2thin[t - 1], t - 1] * thin[t - 1]
      } else{
        thinsp_id <- which(thin_order == sp2thin[t - 1])
        if(thin[t - 1] == 0){
          N[sp2thin[t - 1], t - 1] <- min(N[sp2thin[t - 1], t - 1], thin_levels[thinsp_id])
        }
      }

      # step the process forward
      nesp <- (1:nsp)[-extinct]
      if(length(extinct) > 0){
        N[nesp, t] <- N[nesp, t - 1] * lambdas[nesp] *
          exp(- A_mat[nesp, nesp] %*% N[nesp, t - 1] + rnorm(length(nesp)) * sigmas[nesp] / sqrt(N[nesp, t - 1]))
      } else{
        N[, t] <- N[, t - 1] * lambdas * exp(- A_mat %*% N[, t - 1] + rnorm(nsp) * sigmas / sqrt(N[, t - 1]))
      }

    }

  }


  ### Disturbance but no thinning ###
  if(thin_freq == 0 & dist_prob > 0){

    # create vector of disturbance events
    dist <- as.numeric(purrr::rbernoulli(steps, p = dist_prob))

    # create vectors for which species get disturbed
    dist_vecs <- purrr::map(
      1:steps,
      ~ as.numeric(purrr::rbernoulli(nsp, p = 1 - prop_cdist))
    )

    # change the disturbed values to a multiplier
    dist_vecs <- lapply(
      dist_vecs,
      FUN = function(x, p){
        x[x == 0] <- p
        x
      },
      p = (1 - dist_int)
    )

    # now proceed with simulations
    for(t in 2:steps){

      # tracking extinct species
      extinct <- which(round(as.double(N[, t - 1])) == 0)
      N[extinct, t - 1] <- 0

      # disturb if it was a disturbance year
      N_tm1 <- N[, t - 1] * (1 - dist[t - 1]) + N[, t - 1] * dist[t - 1] * dist_vecs[[t - 1]]


      # step the process forward
      nesp <- (1:nsp)[-extinct]
      if(length(extinct) > 0){
        N[nesp, t] <- N[nesp, t - 1] * lambdas[nesp] *
          exp(- A_mat[nesp, nesp] %*% N[nesp, t - 1] + rnorm(length(nesp)) * sigmas[nesp] / sqrt(N[nesp, t - 1]))
      } else{
        N[, t] <- N[, t - 1] * lambdas * exp(- A_mat %*% N[, t - 1] + rnorm(nsp) * sigmas / sqrt(N[, t - 1]))
      }

    }

  }

  if(dist_prob > 0){
    return(list(
      N = N,
      disturbances = lapply(
        1:steps,
        FUN = function(x, dvecs, indic){
          (1 - dvecs[[x]]) * indic[x]
        },
        dvecs = dist_vecs,
        indic = dist
      )
    ))
  } else{
    return(N)
  }

}









#' Plot a simulated community
#'
#' @param N The matrix of species densities in rows through time over columns.
#'
#' @return ggplot
#'
plot_comm <- function(N){
  library(ggplot2)
  nsp <- nrow(N)
  steps <- ncol(N)

  df <- data.frame(
    t = rep(1:steps, each = nsp),
    sp = as.factor(rep(1:nsp, steps)),
    N = as.vector(N)
  )

  return(
    ggplot(data = df, aes(x = t, y = N, color = sp)) +
      geom_line() +
      theme_classic()
  )
}





