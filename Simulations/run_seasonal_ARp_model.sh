#!/bin/bash -l

#SBATCH --account=modelscape
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alice.carter@flbs.umt.edu
#SBATCH --job-name=ARp_sim1
#SBATCH --mem=30G
#SBATCH --time=96:00:00
#SBATCH -o slurmlogs/slurm_%A%a.out
#SBATCH -e slurmlogs/slurm_%A%a.err
#SBATCH --array=1-15

# Set the parameter combination to use and generate names of R scripts and log file
Rscript=fit_AR-p_seasonal_model_sims.R
OutDir='ARp_err_sims_9_18'

# Change to the relevant working directory
mkdir /project/modelscape/analyses/sponges/Data/aquatic_sim_data/$OutDir
cd /project/modelscape/analyses/sponges/Simulations

# Load R and MPI
module load arcc/1.0 gcc/12.2.0 r/4.2.2

Rscript --vanilla $Rscript $OutDir $SLURM_ARRAY_TASK_ID
