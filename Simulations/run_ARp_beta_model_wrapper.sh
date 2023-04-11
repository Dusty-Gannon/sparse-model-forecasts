#!/bin/bash -l

# Set the parameter combination to use and generate names of R scripts and log file
OutDir='AR_p_sims_model_output_4_11'

module load arcc/1.0 gcc/12.2.0 r/4.2.2

sbatch --wait run_ARp_beta_model.sh $OutDir

# define the second R script that will be run
Rscript=combine_ARp_sim_output.R

Rscript --vanilla $Rscript $OutDir

# remove directory of individual files
#rm -r /project/modelscape/analyses/sponges/Data/aquatic_sim_data/$OutDir
