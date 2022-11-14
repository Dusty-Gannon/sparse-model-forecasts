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
#' @param t Time step at which to assess whether the species meet the desired specs.
#' @param num_ngs Number of non-generic competitors.
#' @param rare Logical telling the function whether or not to pick the rarest species
#' that meets the other criteria or the most common.
#'
#' @return Column name (integer as a character) of the focal species.
#'
#' @examples
#'
choose_focal <- function(df, col_ids, t, num_ngs, rare = F, time_colname = "t"){

  # get list of potential species
  pot_sp <- names(df)[col_ids]

  # get time column index
  time_col <- which(names(df) == time_colname)

  # set order for sorting based on rare vs. common
  if(isTRUE(rare)){decreasing <- F} else {decreasing <- T}

  # loop through each to see if they meet the criteria
  foc <- NULL
  counter <- 1
  while(is.null(foc) & counter <= length(col_ids)){
    for(i in 1:length(pot_sp)){
      # get list of competitors
      compts <- as.integer(pot_sp[i]) + 1:num_ngs
      compts_ids <- which(names(df) %in% as.character(compts))
      if(
        {as.double(df[t, col_ids[i]]) ==
            sort(as.double(df[t, col_ids]), decreasing = decreasing)[counter]} &
        {sum(compts %in% as.integer(names(df[,-time_col]))) > 0} &
        {as.double(df[t, col_ids[i]]) > 0} &
        {sum(as.double(df[t, compts_ids])) > 0}
      ){
        foc <- pot_sp[i]
      }
    }
    counter <- counter + 1
  }

  if(is.null(foc)){
    stop("No species in the dataset fits the bill...")
  } else{
    return(foc)
  }

}
