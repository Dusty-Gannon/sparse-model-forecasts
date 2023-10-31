################################################################
# This script can be used to find out which model fits failed by
# supplying a relative path to a directory with the results and
# the range of identifiers to check
################################################################

library(here)
library(stringr)

args <- commandArgs(trailingOnly = T)

file_list <- list.files(here(args[1]))

ids <- seq(
  as.numeric(args[2]), as.numeric(args[3])
)

# extract the ids that are present
ids_present <- sort(as.numeric(
  str_extract_all(file_list, pattern = "[[0-9]]+", simplify = T)
))

# print a list of ids that are not present
print(
  ids[!(ids %in% ids_present)]
)
