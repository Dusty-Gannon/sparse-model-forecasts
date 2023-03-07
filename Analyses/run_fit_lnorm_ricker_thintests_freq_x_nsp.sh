#!/bin/bash -l

#SBATCH --account=modelscape
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=21
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dgannon@uwyo.edu

# Set the parameter combination to use and generate names of R scripts and log file
Rscript=fit_lnorm_ricker_mods_parallel.R
LogFile=ricker_thin2

# Change to the relevant working directory
cd /project/modelscape/analyses/sponges/Analyses

# Load R
module load arcc/1.0  gcc/12.2.0  r/4.2.2

Rscript --vanilla $Rscript lnorm_ricker_thin_freq_x_nsp_ordered_S5_s55.rds 51 500 freq_x_nsp_thintests_mfits.rds > $LogFile
