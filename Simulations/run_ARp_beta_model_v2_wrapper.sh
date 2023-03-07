#!/bin/bash -l

# Set the parameter combination to use and generate names of R scripts and log file
OutDir=ARp_sims
nsteps=370
sigma=1

sbatch --wait run_ARp_beta_model_v2.sh $OutDir $nsteps $sigma

# define the second R script that will be run
Rscript=combine_ARp_sim_output.R

Rscript --vanilla $Rscript $OutDir $nsteps $sigma

# remove directory of individual files
rm -r /project/modelscape/analyses/sponges/Data/aquatic_sim_data/$OutDir
