

#' Generate an array of parameters that determine competition
#'
#' This function generates an array of parameters that, when combined with
#' a model matrix of environmental variables, determines the matrix of
#' competition coefficients.
#'
#' @param nsp Number of species in the community
#' @param num_ngs Number of non-generic competitors (usually, stronger competitors)
#' @param env Parameters determining the linear effect of one or more environmental variables
#' on the competitive effects of some select dynamic species whose competitive abilities vary
#' with the environment
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
#'   env = c(-1, 1), num_dynamic = 2,
#'   generic = 0.001
#' )
#'
ricker_comp_array <- function(
    nsp, num_ngs, env,
    num_dynamic, generic = 0.001,
    ng_strength = c(20, 50), intra_strength = c(80, 100)
){

  # randomly sample the species that get more or less competitive
  #  depending on the environment
  dyn_sp <- sample(1:nsp, size = num_dynamic)
  m <- length(env)

  # initialize array
  B <- array(data = 0, dim = c(nsp, m + 2, nsp))

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
          B[i, (1:m) + 2, j] <- env
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
          B[i, (1:m) + 2, j] <- env
        }
      }
    }

    if(nsp - j == 0){
      B[1:num_ngs, 2, j] <- generic * runif(num_ngs, min = 20, max = 50)
      for(i in (1:num_ngs)){
        if(i %in% dyn_sp){
          B[i, (1:m) + 2, j] <- env
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







#' Ricker competition model
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
ricker_step <- function(N_t, lambdas, A_mat, stochastic = T){

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


