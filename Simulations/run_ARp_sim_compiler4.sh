#!/bin/bash -l
#SBATCH --account=modelscape
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alice.carter@flbs.umt.edu
#SBATCH --job-name=ARp_sim_compile
#SBATCH -o slurm_compile4.out
#SBATCH --mem=30G
#SBATCH --time=48:00:00


# Set the parameter combination to use and generate names of R scripts and log file
OutDir='ARp_err_sims_10_31'

module load arcc/1.0 gcc/12.2.0 r/4.4.0

# define the second R script that will be run
Rscript=combine_seasonal_ARp_sim_output.R
Rscript --vanilla $Rscript $OutDir

# remove directory of individual files
#rm -r /project/modelscape/analyses/sponges/Data/aquatic_sim_data/$OutDir
