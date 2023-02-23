#!/bin/bash -l

#SBATCH --account=modelscape
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=21
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dgannon@uwyo.edu

# Set the parameter combination to use and generate names of R scripts and log file
Rscript=fit_Pois_Ricker_terr_comm.R
LogFile=log_ricker_n25_lambda_known

# Change to the relevant working directory
cd /project/modelscape/analyses/sponges/Analyses

# Load R and MPI
module load arcc/1.0  gcc/12.2.0 r/4.2.2

Rscript --vanilla $Rscript 25 50 "known" > $LogFile
