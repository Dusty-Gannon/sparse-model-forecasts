#!/bin/bash -l

#SBATCH --account=modelscape
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dgannon@uwyo.edu

# Set the parameter combination to use and generate names of R scripts and log file
Rscript=fit_lnorm_ricker_3scenarios_tandem.R
LogFile=ricker_vrtests_tandem

# Change to the relevant working directory
cd /project/modelscape/analyses/sponges/Analyses

# Load R and MPI
module load arcc/1.0  gcc/12.2.0 r/4.2.2

Rscript --vanilla $Rscript lnorm_ricker_sims_200steps_rho0_S3_s57.rds > $LogFile
