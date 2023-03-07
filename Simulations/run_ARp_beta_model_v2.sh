#!/bin/bash -l

#SBATCH --account=modelscape
#SBATCH --nodes=1
#SBATCH --time=00:30:00
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=acarte26@beartooth.arcc.uwyo.edu
#SBATCH --job-name=ARp_sim1
#SBATCH --mem=10G
#SBATCH -o slurmlogs/slurm_%A.out
#SBATCH -e slurmlogs/slurm_%A.err
#SBATCH --array=1-2

# Set the parameter combination to use and generate names of R scripts and log file
Rscript=AR-p_beta_p_model_sims_v2.R
LogFile=ARpmod_sims.log
# OutDir=ARp_sims_test

# Change to the relevant working directory
mkdir /project/modelscape/analyses/sponges/Data/aquatic_sim_data/$1
cd /project/modelscape/analyses/sponges/Simulations

# Load R and MPI
module load arcc/1.0 gcc/12.2.0 r/4.2.2

Rscript --vanilla $Rscript $1 $SLURM_ARRAY_TASK_ID
