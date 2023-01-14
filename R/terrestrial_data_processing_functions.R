#######################################
# Terrestrial data processing functions
#######################################

#' Choose focal species
#'
#' This function goes through a dataframe (\code{df}) with one column per species
#' and species abundances through time in the rows and selects a focal species given
#' a time window in the series and characteristics of the species. The focal species
#' must be present for at least 95% of the time window, have at least one non-generic
#' competitor present for 95% of the time window, and may be a rare or a common species
#' depending on the \code{rare} flag.
#'
#' @param df A dataframe of species abundance through time in rows and different species
#' across columns.
#' @param num_ngs Number of non-generic competitors.
#' @param tw Length of time window over which to assess whether the species meet the desired specs.
#' @param start Time step at which to begin looking for an appropriate species.
#' @param rare Logical (Default = \code{FALSE}). If set to \code{TRUE}, the function will look for an
#' appropriate rare species to choose as a focal. Otherwise, preference is given to the common species.
#'
#' @return Column name (integer as a character) of the focal species.
#'
#' @examples
#'
choose_focal <- function(
    df, num_ngs,
    tw = 100, start = 50,
    rare = F
){

  # store some useful variables
  nsp <- ncol(df)

  # set order for sorting based on rare vs. common
  if(isTRUE(rare)){decreasing <- F} else {decreasing <- T}

  # loop through each to see if they meet the criteria
  foc <- NULL
  l <- start + 1
  r <- start + tw

  # set initial species order
  sp_order <- (1:nsp)[
    order(apply(df[l:r, ], 2, mean), decreasing = decreasing)
  ]

  while(is.null(foc) & r <= nrow(df)){

    # start with the first species
    i <- 1

    while(i <= nsp){

      # get list of competitors
      if(nsp - sp_order[i] >= num_ngs){
        compts <- sp_order[i] + 1:num_ngs
      }
      if(nsp - sp_order[i] < num_ngs & nsp - sp_order[i] > 0){
        compts <- c(
          1:(num_ngs - nsp + sp_order[i]),
          sp_order[i] + 1:(nsp - sp_order[i])
        )
      }
      if(nsp - sp_order[i] == 0){
        compts <- 1:num_ngs
      }

      # check if the current species meets the criteria
      if(
        {mean(as.double(df[l:r, sp_order[i]]) != 0) >= 0.95} &
        {sum(apply(
          df[l:r, compts], 2, function(x){mean(x != 0) > 0.95}
        )) > 0}
      ){
        foc <- sp_order[i]
        i <- nsp + 1
      } else{
        # check the next species
        i <- i + 1
      }
    }

    # if no species was found for that time window, shift the window forward
    # and reset the species order
    if(is.null(foc)){
      l <- l + (tw + tw %% 2) / 2
      r <- l + tw - 1
      sp_order <- (1:nsp)[
        order(apply(df[l:r, ], 2, mean), decreasing = decreasing)
      ]
    }

  }

  if(is.null(foc)){
    stop("No species in the dataset fits the bill...")
  } else{
    return(
      list(foc = foc, tw = l:r)
    )
  }

}






#' Choose a focal species from the communities simulated using a Ricker model
#'
#' @param df Data frame (or matrix) with one column per species and each species' abundance
#' in each time step in the rows.
#' @param num_ngs Integer number of non-generic competitors for each species
#' @param dyn_sp Integer vector with the identities of the species that become more competitive
#' with changes in the environment
#' @param tw Integer vector with the start and end of the time window to be assessed
#' @param perc_present Vector indicating the proportion of observations for which the focal
#' (\code{perc_present[1]}) and the non-generic competing species (\code{perc_present[2]}) should
#' be present for the species to be used as the focal species.
#'
#' @return List in which \code{list$value} is the index for the focal species.
#'
#' @examples
#'
choose_focal2 <- function(
    df, num_ngs, dyn_sp,
    tw = c(101, 200),
    perc_present = c(0.95, 0.7)
){

  # store some useful variables
  nsp <- ncol(df)
  n_obs <- tw[2] - tw[1] + 1

  # loop through each to see if they meet the criteria
  foc <- NULL
  counter <- 1
  while(is.null(foc) & counter <= nsp){
    for(i in 1:nsp){
      # get list of competitors
      if(nsp - i >= num_ngs){
        compts <- i + 1:num_ngs
      }
      if(nsp - i < num_ngs & nsp - i > 0){
        compts <- c(
          1:(num_ngs - nsp + i),
          i + 1:(nsp - i)
        )
      }
      if(nsp - i == 0){
        compts <- 1:num_ngs
      }

      # check if species i has ng competitors present and one of them is a dynamic species
      if(
        {mean(df[tw[1]:tw[2], i]) > perc_present[1]} &
        {sum(df[floor(perc_present[2] * n_obs + tw[1] - 1), compts]) > 0} &
        {
          sum(which(which(df[floor(perc_present[2] * n_obs + tw[1] - 1), ] > 0) %in% compts) %in% dyn_sp) > 0
        }
      ){
        foc <- i
      }
    }
    counter <- counter + 1
  }

  # if there isn't a species with at least one dynamic competitor present,
  #   select one that has at least one ng competitor present
  if(is.null(foc)){
    counter <- 1
    while(is.null(foc) & counter <= nsp){
      for(i in 1:nsp){
        # get list of competitors
        if(nsp - i >= num_ngs){
          compts <- i + 1:num_ngs
        }
        if(nsp - i < num_ngs & nsp - i > 0){
          compts <- c(
            1:(num_ngs - nsp + i),
            i + 1:(nsp - i)
          )
        }
        if(nsp - i == 0){
          compts <- 1:num_ngs
        }

        # check if species i has ng competitors present and one of them is a dynamic species
        if(
          {mean(df[tw[1]:tw[2], i]) > perc_present[1]} &
          {sum(df[floor(perc_present[2] * n_obs + tw[1] - 1), compts]) > 0}
        ){
          foc <- i
        }
      }
      counter <- counter + 1
    }
  }

  # if there is still no species that seems good, exit
  if(is.null(foc)){
    stop("No species in the dataset fits the bill...")
  } else{
    return(
      list(foc = foc, tw = tw)
    )
  }

}







