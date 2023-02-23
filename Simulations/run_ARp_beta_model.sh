#!/bin/bash -l

#SBATCH --account=modelscape
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=21
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=acarte26@beartooth.arcc.uwyo.edu
#SBATCH --job-name=ARp_sim1
#SBATCH --mem=10G

# Set the parameter combination to use and generate names of R scripts and log file
Rscript=AR-p_beta_p_model_sims.R
LogFile=ARpmod_sims.log

# Change to the relevant working directory
cd /project/modelscape/analyses/sponges/Simulations

# Load R and MPI
module load gcc/7.3.0 swset/2018.05 r/3.5.3

Rscript --vanilla $Rscript 2 > $LogFile
