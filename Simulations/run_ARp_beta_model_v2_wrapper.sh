#!/bin/bash -l

# Set the parameter combination to use and generate names of R scripts and log file
OutDir=ARp_sims_test

sbatch run_ARp_beta_model_v2.sh --export=OutDir

# define the second R script that will be run
Rscript=combine_ARp_sim_output.R

Rscript --vanilla $Rscript $1

# remove directory of individual files
rmdir -r project/modelscape/analyses/sponges/Data/aquatic_sim_data/$1
