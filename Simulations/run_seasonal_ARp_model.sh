#!/bin/bash -l

#SBATCH --account=modelscape
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dustin.gannon@oregonstate.edu
#SBATCH --job-name=ARp_sim1
#SBATCH --mem=24G
#SBATCH --time=01:00:00
#SBATCH -o slurmlogs/slurm_%A%a.out
#SBATCH -e slurmlogs/slurm_%A%a.err
#SBATCH --array=1-2

# Set the parameter combination to use and generate names of R scripts and log file
Rscript=fit_AR-p_seasonal_model_sims.R
OutDir='ARp_err_sims_20260928'

# Change to the relevant working directory
mkdir /project/modelscape/analyses/sparse-model-forecasts/Data/fourier_sims/$OutDir
cd /project/modelscape/analyses/sparse-model-forecasts/Simulations

# Load R and MPI
module load arcc/1.0 gcc/14.2.0 r/4.4.0 r-mass/7.3-59

Rscript --vanilla $Rscript $OutDir $SLURM_ARRAY_TASK_ID
