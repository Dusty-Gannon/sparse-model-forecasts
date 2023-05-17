#!/bin/bash -l

#SBATCH --account=modelscape
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=20G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dgannon@uwyo.edu

# Set the parameter combination to use and generate names of R scripts and log file
Rscript=plot_summarize_Ricker_dist_freq_x_int.R
LogFile=plot_freq_x_nsp

# Change to the relevant working directory
cd /project/modelscape/analyses/sponges/Analyses

# Load modules
module load arcc/1.0  gcc/12.2.0  r/4.2.2

Rscript --vanilla $Rscript lnorm_ricker_dist_freq_x_nsp_x_int_round1.rds round1 > $LogFile





