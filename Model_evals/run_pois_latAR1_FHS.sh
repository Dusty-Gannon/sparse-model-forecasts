#!/bin/bash -l

#SBATCH --account=commbayes
#SBATCH --time=10:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=11
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dgannon@uwyo.edu
#SBATCH --job-name=pois_latAR1

# Set the parameter combination to use and generate names of R scripts and log file
Rscript=pois_latAR1_FHS.R
LogFile=pois_latAR1.log

# Change to the relevant working directory
cd /project/modelscape/analyses/sponges/Model_evals

# Load R and MPI
module load gcc/7.3.0 swset/2018.05 r/3.5.3 r-rstan/2.18.2-py27

R CMD BATCH --no-save --no-restore $Rscript $LogFile
