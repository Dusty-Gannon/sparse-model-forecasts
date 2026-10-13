#!/bin/bash -l
#SBATCH --account=modelscape
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dgannon@uwyo.edu
#SBATCH --job-name=ARp_sim_compile
#SBATCH -o slurm_compile4.out
#SBATCH --mem=30G
#SBATCH --time=02:00:00


# Set the parameter combination to use and generate names of R scripts and log file
cd /project/modelscape/analyses/sparse-model-forecasts/Simulations
OutDir='ARp_err_sims_20260928'

module load arcc/1.0 gcc/14.2.0 r-mass/7.3-59

# define the second R script that will be run
Rscript=combine_seasonal_ARp_sim_output.R
Rscript --vanilla $Rscript $OutDir

# remove directory of individual files
#rm -r /project/modelscape/analyses/sponges/Data/aquatic_sim_data/$OutDir
