#!/bin/bash -l

#SBATCH --account=modelscape
#SBATCH --time=05:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dgannon@uwyo.edu
#SBATCH --job-name=lottery_50s
#SBATCH --mem=8G

# Set the parameter combination to use and generate names of R scripts and log file
Rscript=teton_lottery_model_sims.R
LogFile=lotmod_sims.log

# Change to the relevant working directory
cd /project/modelscape/analyses/sponges/Simulations

# Load R and MPI
module load gcc/7.3.0 swset/2018.05 r/3.5.3

R CMD BATCH --no-save $Rscript $LogFile
