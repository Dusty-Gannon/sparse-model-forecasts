########################################################
# This script is meant to run on Beartooth to split
# simulation data files that contain lists of results
# from replicate simulations into separate files for each
# result. This is handy for fitting with array jobs.
# **NOTE**: The data files must contain the same number
# of results for the numbering to make sense and not
# overwrite files.
########################################################

library(here)

# Arguments include the path to the directory with the
# data files that are to be split and, optionally, an
# output directory. If no output is supplied, the files
# are stored in the input directory
args <- commandArgs(trailingOnly = T)
if(length(args) == 1){
  args <- c(args, args)
}

file_list <- list.files(here(args[1]))

# loop through and create unique numbering and save
# unique files
for(i in 1:length(file_list)){
  dat_i <- readRDS(here(paste0(args[1], file_list[i])))
  for(j in 1:length(dat_i)){
    dat_ij <- dat_i[[j]]
    name_ij <- paste(
      sub("_rep.*", "", x = file_list[i]),
      (length(dat_i) * (i - 1)) + j,
      sep = "_"
    )
    saveRDS(
      dat_ij,
      file = paste0(
        args[2],
        name_ij,
        ".rds"
      )
    )
  }
}
