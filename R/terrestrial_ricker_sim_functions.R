

#' Generate matrix of competition parameters
#'
#' This function generates a matrix of parameters that, when combined with
#' a model matrix with heterospecific abundances and environmental variables,
#' determines the strength of competition acting on the population of a focal
#' species.
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
#' @return A matrix with \eqn{S \times V + 1} rows and \eqn{S} columns, where \eqn{S} is the
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
ricker_comp_matrix <- function(
    nsp, num_ngs, env,
    num_dynamic, generic = 0.001,
    ng_strength = c(20, 50), intra_strength = c(80, 100)
){

  # randomly sample the species that get more or less competitive
  #  depending on the environment
  dyn_sp <- sample(1:nsp, size = num_dynamic)
  m <- length(env)

  # initialize matrix
  B <- matrix(data = 0, nrow = nsp * (1 + length(env)) + 1, ncol = nsp)

  # fill in the first row with the generic effect
  B[1, ] <- rep(generic, nsp)

  # proceed to fill in the remaining elements
  for(j in 1:nsp){

    # intraspecific effect
    B[(j + 1), j] <- generic * runif(1, min = intra_strength[1], max = intra_strength[2])

    # non generic effects
    if(nsp - j >= num_ngs){
      B[(j + 1) + c(1:num_ngs), j] <- generic * runif(num_ngs, min = ng_strength[1], max = ng_strength[2])
      for(i in (j + 1) + c(1:num_ngs)){
        for(v in 1:m){
          if((i - 1) %in% dyn_sp){B[nsp * v + i, j] <- env[v]}
        }
      }
    }

    # begin wrapping once we start to reach the edge
    if(nsp - j < num_ngs & nsp - j > 0){
      # get ng_sp ids
      ng_ids <- c(
        1:(num_ngs - (nsp - j)) + 1,
        (j + 1) + c(1:(nsp - j))
      )
      B[ng_ids, j] <- generic * runif(num_ngs, min = 20, max = 50)
      for(i in ng_ids){
        for(v in 1:m){
          if((i - 1) %in% dyn_sp){B[nsp * v + i, j] <- env[v]}
        }
      }
    }

    if(nsp - j == 0){
      B[1:num_ngs + 1, j] <- generic * runif(num_ngs, min = 20, max = 50)
      for(i in (1:num_ngs) + 1){
        for(v in 1:m){
          if((i - 1) %in% dyn_sp){B[nsp * v + i, j] <- env[v]}
        }
      }
    }

  }

  return(B)

}




