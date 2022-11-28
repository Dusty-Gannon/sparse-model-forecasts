#!/bin/bash -l

#SBATCH --account=modelscape
#SBATCH --time=05:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=51
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dgannon@uwyo.edu

# Set the parameter combination to use and generate names of R scripts and log file
Rscript=fit_ARMA-pq_terr_comm_sims.R
LogFile=ARMA-pq.log

# Change to the relevant working directory
cd /project/modelscape/analyses/sponges/Analyses

# Load R and MPI
module load gcc/7.3.0 swset/2018.05 r-rstan/2.18.2-py27

Rscript --vanilla $Rscript 100 300 1 1 > $LogFile
