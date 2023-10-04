#!/bin/bash -l

#SBATCH --account=modelscape
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=11
#SBATCH --cpus-per-task=1
#SBATCH --mem=20G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dgannon@uwyo.edu
#SBATCH --array=1-100%50
#SBATCH --out=./logfiles/slurm-%A_%a.out

# Set the parameter combination to use and generate names of R scripts and log file
Rscript=lnorm_Ricker_model_sims_freq_x_nsp_thin_s4t.R
LogFile_pref=./logfiles/thin_sims_rep

# Change to the relevant working directory
cd /project/modelscape/analyses/sponges/Simulations

# Load R and MPI
module load arcc/1.0  gcc/12.2.0  r/4.2.2

# create task ids
job=()
for i in ${SLURM_ARRAY_TASK_ID[@]}; do
  job+=$i
done

Rscript --vanilla $Rscript ${job} "Data/terrestrial_sim_data/lnorm_ricker/thin_sims_s4t/sim_reps/" > "${LogFile_pref}_${job}"




