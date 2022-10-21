

#' Build a competition matrix
#'
#' This function builds an \eqn{N \times N} matrix,
#' where \eqn{N} is the number of species in the system,
#' of competition coefficients \eqn{A}, where
#' \eqn{A_{ij}} is the competitive effect of species \eqn{j} on \eqn{i}.
#' \eqn{A_{ij} = \alpha_i (\rho + (1 - \rho) \delta_{ij})}, where \eqn{\delta_{ij}=1} if
#' \eqn{i=j} and 0 otherwise and \eqn{\alpha_i} is the strength of self regulations for
#' species i.
#'
#' @param n_sp The number of species/populations competing
#' @param rho The strength of interspecific competition relative to intraspecific competion.
#' \eqn{\rho \in (0,1)}
#' @param alpha A scalar or vector of length \code{n_sp} with the strength(s) of
#' intraspecific competition.
#' @param num_ngs Number of strong competitors
#'
#'
#' @return \code{n_sp} X \code{n_sp} matrix
#' @export
#'
#' @examples
#' comp_matrix(3, rho=0.4, alpha=0.01, num_ngs = 1)
#'
comp_matrix <- function(n_sp, rho=0, alpha, num_ngs){
  if(rho < 0 | rho > 1){
    stop("rho must be on the unit interval (0,1)")
  }

  Id <- diag(nrow = n_sp, ncol = n_sp)
  Dalpha <- diag(alpha, nrow = n_sp, ncol = n_sp)
  mat <- Dalpha%*%(rho + (1-rho)*Id)

   if(num_ngs >= 1){
    # pull out top row (focal species)
    tr <- mat[1,]
    tr[2:(num_ngs+1)] <- alpha[1] - runif(num_ngs, min = 0, max = alpha[1]/4)

    # make closer to symmetric
    fc <- mat[, 1]
    fc[2:(num_ngs+1)] <- alpha[2:(num_ngs + 1)] - runif(num_ngs, min = alpha[2:(num_ngs + 1)]/4, max = alpha[2:(num_ngs + 1)]/3)

    mat[1,] <- tr
    mat[, 1] <- fc


   }
  return(mat)
}











#' Count number of neighbors at a given radius in a lattice
#'
#' @param M Number of rows in the lattice
#' @param J Number of columns in the lattice
#' @param r Desired radius of the neighborhood
#'
#' @return An \code{M}\eqn{\times}\code{J} matrix with the number of neighbors for each cell
#' @export
#'
#' @examples
#' count_neighbors(9, 9, r = 2)
#'
count_neighbors <- function(M, J, r){

  # neighborhood size
  n_size <- (2 * r + 1)^2 - 1

  # function to determine how many rows/columns
  #   get cut from the full neighborhood when at the
  #   edge of the lattice
  rc_cut <- function(r, i, M){
    return(
      min(c(r+1, i, M - i + 1))
    )
  }

  # function to compute the neighborhood size (even at margins)
  n_neighbors <- function(n, r, t, s){
    n - (r - t + 1)*(2*r + 1) - (r - s + 1)*(2*r + 1) + (r - t + 1)*(r - s + 1)
  }

  # put row and column indices into a dataframe
  df_long <- data.frame(
    row_index = rep(1:M, J),
    col_index = rep(1:J, each = M)
  )

  # now add t and s values
  df_long$t <- purrr::map_dbl(
    df_long$row_index,
    ~ rc_cut(r = r, i = .x, M = M)
  )
  df_long$s <- purrr::map_dbl(
    df_long$col_index,
    ~ rc_cut(r = r, i = .x, M = J)
  )
  df_long$nbsize <- purrr::map_dbl(
    1:nrow(df_long),
    ~ n_neighbors(n = n_size, r = r, t = df_long$t[.x], s = df_long$s[.x])
  )

  # reconstruct matrix with number of neighbors for each cell in radius r
  neighbors_mat <- matrix(
    df_long$nbsize,
    nrow = M, ncol = J
  )

  return(neighbors_mat)

}







#' Count neighbors of each type in a neighborhood
#'
#' @param X Matrix with species IDs in each cell in the lattice
#' @param r Radius of the neighborhood (in a lattice)
#' @param sp_list List of species IDs. Should be numeric and starting at 1 (e.g, 1--4 for a 4 species model)
#'
#' @return An array with \eqn{S} slices of matrices of equal dimension to X, where \eqn{S} is the number
#' of species in the model. Each slice represents the count of species \eqn{s=1,2,...,S} in the neighborhood
#' around each cell in X.
#' @export
#'
#' @examples
#' X <- matrix(sample(1:3, 16, replace = T), nrow = 4, ncol = 4)
#' kernel_count(X, r = 1, sp_list = 1:3)
#'
kernel_count <- function(X, r, sp_list){

  # get dimensions of lattice
  M <- nrow(X)
  J <- ncol(X)
  S <- length(sp_list)

  n <- array(dim = c(M, J, S))

  for(i in 1:M){
    for(j in 1:J){
      # range of rows and columns to use for kernel
      cols_lims <- c(
        max(1, j - r),
        min(J, j + r)
      )
      rows_lims <- c(
        max(1, i - r),
        min(M, i + r)
      )

      # create sequence of rows and columns based on limits
      rows_ij <- rows_lims[1]:rows_lims[2]
      cols_ij <- cols_lims[1]:cols_lims[2]


      # create new matrix to avoid counting focal cell
      # dummy matrix
      X_0 <- matrix(data = 0, nrow = nrow(X), ncol = ncol(X))
      X_0[i,j] <- X[i,j]

      # now get count of each species in kernel
      sp_counts <- sapply(
        sp_list,
        FUN = function(x, mat){sum(mat == x)},
        mat = X[rows_ij, cols_ij] - X_0[rows_ij, cols_ij]
      )

      # now fill slice of the array
      n[i, j, ] <- sp_counts
    }
  }

  return(n)
}










#' Fecundity for a row in the lattice
#'
#'
#' @param foc_sp Vector of species occupying each cell in the row in the lattice at the beginning of
#' the time step
#' @param lambda_t Vector of (potentially) time-varying growth rates, one for each species
#' @param alpha Matrix of competition coefficients with as many rows and columns as there are species
#' @param nbrhood Vector of sizes of the neighborhood for each cell in the row
#' @param n matrix of neighbor counts with as many columns as species and as many rows as columns
#' in the lattice
#'
#' @return vector with fecudities for each individual in one row of the lattice
#' @export
#'
#' @examples
#'
fecundity_ll <- function(foc_sp, lambda_t, alpha, nbrhood, n){
  # get vector of growth rates for each cell in row i
  lambda_it <- lambda_t[foc_sp]

  # normalizing constant for neighborhood
  norm_constant <- solve(diag(nbrhood, nrow = length(nbrhood)))

  # pull out row of alpha corresponding to focal ind. sp. id
  alpha_i <- alpha[foc_sp, ]

  # get competition vector
  comp_vec <- 1 + diag(norm_constant %*% n %*% t(alpha_i))

  # fecundity model
  return(
    lambda_it / comp_vec
  )
}





#' Truncated, discretized exponential distribution
#'
#' Probability mass function for dispersal on a grid based on discretizing an exponential distribution.
#'
#' @param d Dispersal distance (\eqn{d=\{0, 1, 2, ...\}}).
#' @param d_max Maximum dispersal distance (truncation). \eqn{d_max = \{0, 1, 2, ...\}}
#' @param lambda Exponential rate parameter (rate at which seeds fall into cells with increasing distance).
#'
#' @return Probability of a seed dispersing to a cell at distance \code{d}.
#' @export
#'
#' @examples
#' P_d_exp(1, d_max = 3, lambda = 1)
#'
P_d_exp <- function(d, d_max, lambda){
  # double check d is within support
  if(d %% 1 != 0){
    stop("d must be an integer")
  }
  # construct normalizing constant
  d_seq <- 0:d_max
  kern_seq <- exp(-lambda * d_seq) - exp(-lambda * (d_seq + 1))
  norm_const <- sum(kern_seq)

  if(d <= d_max){
    # return probability that D = d
    return(
      (exp(-lambda * d) - exp(-lambda * (d + 1)))/norm_const
    )
  } else{
    return(0)
  }
}







#' Seed rain array
#'
#' Construct an \eqn{M \times J \times S} array encoding the total seed rain into each location in the \eqn{M \times J}
#' lattice for each of \eqn{S} species.
#'
#'
#' @param F_mat Matrix of fecundities for each individual in the lattice
#' @param X Lattice with numeric species IDs (1, 2, 3, ..., S) for the individual in each cell
#' @param d_max Max distance (vertical and horizontal, not diagonal) for dispersal (currently
#' set as a global parameter and is not species-specific).
#' @param rate Exponential rate parameter(s) for each species defining the dispersal kernel. This
#' can be a vector of length \eqn{S}, where \eqn{S} is the number of species in the model, or a scalar,
#' giving all species equal dispersal distributions.
#'
#' @return An array with \eqn{S} slices made up of \eqn{M \times J} matrices, where \eqn{M \times J}
#' is the dimension of the lattice.
#' @export
#'
#' @examples
#' # lattice
#' X <- matrix(sample(1:3, 16, replace = T), nrow = 4)
#'
#' # fecundity for each individual
#' F_mat <- matrix(rgamma(16, shape = 10), nrow = 4)
#'
#' # max dispersal distance
#' d_max <- 2
#'
#' # exponential dispersal rates
#' rate <- c(0.5, 1, 2)
#'
#' # compute seed rain array
#' seed_rain_array(F_mat, X, d_max, rate)
#'
#'
seed_rain_array <- function(F_mat, X, d_max, rate = 1, nsp){

  # store some useful values
  S <- nsp                          # number of species
  M <- nrow(X)                      # number of rows in lattice
  J <- ncol(X)                      # number of columns in lattice

  # convert rate to a vector if it isn't already
  if(length(rate) == 1){
    rate <- rep(rate, S)
  }

  # now count neighbors at each distance for each cell
  nbrs_array <- lapply(
    1:d_max,
    function(x, M, J){count_neighbors(M = M, J = J, r = x)},
    M = M,
    J = J
  )

  # create seed rain into each cell by each species
  seeds <- array(data = 0, dim = c(nrow(X), ncol(X), S))
  for(i in 1:M){
    for(j in 1:J){
      # create rings of seed rain
      dummy_mat <- lapply(
        0:d_max,
        FUN = function(x, M, J){
          matrix(data = 0, nrow = M, ncol = J)
        },
        M = nrow(X),
        J = ncol(X)
      )
      dummy_mat[[1]][i,j] <- 1
      for(d in 1:d_max){
        # range of rows and columns to use for kernel
        cols_lims <- c(
          max(1, j - d),
          min(J, j + d)
        )
        rows_lims <- c(
          max(1, i - d),
          min(M, i + d)
        )

        # create sequence of rows and columns based on limits
        rows_ij <- rows_lims[1]:rows_lims[2]
        cols_ij <- cols_lims[1]:cols_lims[2]

        # add 1's to the array
        dummy_mat[[d+1]][rows_ij, cols_ij] <- 1

        # now make array of the "smaller rings"
        # i.e., the hole to remove from the doughnut
        hole <- Reduce(
          '+',
          dummy_mat[1:(d)]
        )
        dummy_mat[[d+1]] <- dummy_mat[[d+1]] - hole

      }

      # Given the list of "rings", multiply the rings by prob. dispersal,
      #  fecudity, and divide by count of cells in a ring
      disp_map_ijd <- vector(mode = "list", length = length(dummy_mat))
      disp_map_ijd[[1]] <- F_mat[i, j] * dummy_mat[[1]] * P_d_exp(d = 0, d_max, lambda = rate[X[i, j]])
      for(d in 1:d_max){
        disp_map_ijd[[d + 1]] <-
          F_mat[i, j] * dummy_mat[[d + 1]] * P_d_exp(d = d, d_max, lambda = rate[X[i, j]])/sum(dummy_mat[[d + 1]])
      }

      seeds[, , X[i, j]] <- seeds[, , X[i, j]] + Reduce('+', disp_map_ijd)

    }
  }

  return(seeds)

}









#' Death and replacement
#'
#' This function completes the final stage of the lottery model simulations by sampling
#' individuals of each species to die and be replaced by an individual of a new species
#' with probability proportional to the amount of seed rain into a given cell from each
#' species.
#'
#'
#' @param X Lattice with \eqn{M} rows and \eqn{J} columns and species ids (\eqn{1,2,...,S})
#' in each cell.
#' @param prob_death A vector of length \eqn{S} with the probability a given individual dies in
#' a given time step for each of \eqn{S} species, or a single probability for fixed death rates
#' across species.
#' @param seed_rain \eqn{M} x \eqn{J} x \eqn{S} array with the seed rain falling into each cell of
#' the lattice from each of \eqn{S} species.
#'
#' @return New lattice with some adults replaced by seedlings, potentially of a new species.
#' @export
#'
#' @examples
#' # define inputs
#' X <- matrix(sample(1:3, 16, replace = T), nrow = 4)
#' prob_death <- c(0.1, 0.05, 0.2)
#' seed_rain_list <- purrr::map(1:3, ~ matrix(rgamma(16, shape = 10), nrow = 4))
#' seed_rain <- abind::abind(seed_rain_list[[1]], seed_rain_list[[2]], seed_rain_list[[3]], along = 3)
#'
#' # simulate one round of death and replacement
#' die_replace(X = X, prob_death = prob_death, seed_rain = seed_rain)
#'
die_replace <- function(X, prob_death, seed_rain){

  # store some useful variables
  M <- nrow(X)            # number of rows in the lattice
  J <- ncol(X)            # number of columns in the lattice
  S <- dim(seed_rain)[3]  # number of competing species

  # convert prob_death to a vector if not already
  if(length(prob_death) == 1){
    prob_death <- rep(prob_death, S)
  }

  # get coordinates of each species in the lattice
  sp_coords <- lapply(
    1:S,
    FUN = function(x, lattice){
      which(lattice == x, arr.ind = T)
    },
    lattice = X
  )

  # count individuals from each species
  sp_counts <- sapply(
    sp_coords,
    FUN = nrow
  )

  # kill some off adults
  deaths <- lapply(
    1:S,
    FUN = function(x, count, p){
      rbinom(count[x], size = 1, p = p[x])
    },
    count = sp_counts,
    p = prob_death
  )

  sp_coords_deaths <- lapply(
    1:S,
    FUN = function(x, coords, deaths){
      coords[[x]][which(deaths[[x]] == 1), ]
    },
    coords = sp_coords,
    deaths = deaths
  )
  # make sure each element of sp_coords_deaths is a matrix
  sp_coords_deaths <- lapply(sp_coords_deaths, matrix, ncol = 2)

  # Replace with seeds from seed_rain array.
  # first write a function to replace the value of specific cells
  for(s in 1:S){
    if(nrow(sp_coords_deaths[[s]]) > 0){
      replacements <- vector(mode = "double", length = nrow(sp_coords_deaths[[s]]))
      for(i in 1:length(replacements)){
        replmnt_i <- sample(
          1:S,
          size = 1,
          prob = seed_rain[sp_coords_deaths[[s]][i,1], sp_coords_deaths[[s]][i,2], ] /
            sum(seed_rain[sp_coords_deaths[[s]][i,1], sp_coords_deaths[[s]][i,2], ])
        )
        # replace the individual
        X[sp_coords_deaths[[s]][i,1], sp_coords_deaths[[s]][i,2]] <- replmnt_i
      }
    }
  }

  return(X)


}








