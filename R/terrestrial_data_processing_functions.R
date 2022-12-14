#######################################
# Terrestrial data processing functions
#######################################

#' Choose focal species
#'
#' This function goes through a list of potential focal species and
#' selects one based on the desired characteristics provided as arguments.
#'
#' @param df A dataframe of species cover through time in rows and different species
#' across columns.
#' @param col_ids Vector of column indexes for the potential candidate species.
#' @param num_ngs Number of non-generic competitors.
#' @param tw Length of time window over which to assess whether the species meet the desired specs.
#' @param start Time step at which to begin looking for an appropriate species
#'
#' @return Column name (integer as a character) of the focal species.
#'
#' @examples
#'
choose_focal <- function(
    df, col_ids, num_ngs,
    exclude_names, tw = 100, start = 50,
    rare = F
){

  # get list of potential species
  pot_sp <- names(df)[col_ids]

  # get time and environment column indexes
  exclude_cols <- which(names(df) %in% exclude_names)

  # set order for sorting based on rare vs. common
  if(isTRUE(rare)){decreasing <- F} else {decreasing <- T}

  # loop through each to see if they meet the criteria
  foc <- NULL
  counter <- 1
  l <- start + 1
  r <- start + tw
  while(is.null(foc) & counter <= length(col_ids)){
    for(i in 1:length(pot_sp)){
      # get list of competitors
      compts <- as.integer(pot_sp[i]) + 1:num_ngs
      compts_ids <- which(names(df) %in% as.character(compts))
      while(r <= nrow(df)){
        if(
          {mean(as.double(df[l:r, col_ids[i]]) != 0) >= 0.7} &
          {sum(compts %in% as.integer(names(df[, -exclude_cols]))) > 0} &
          {sum(apply(
            df[l:r, compts_ids], 2, function(x){mean(x != 0) > 0.5}
          ))}
        ){
          foc <- pot_sp[i]
          r <- nrow(df) + 1
        } else{
          # slide the window forward
          l <- (r + (tw/2) %% 1) + 1
          r <- l + tw - 1
        }
      }
    }
    counter <- counter + 1
  }

  if(is.null(foc)){
    stop("No species in the dataset fits the bill...")
  } else{
    return(
      list(foc = foc, tw = l:(l + tw - 1))
    )
  }

}
