#!/bin/bash -l

#SBATCH --account=modelscape
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=5
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dgannon@uwyo.edu
#SBATCH --array=1-5000%50
#SBATCH --out=./logfiles/slurm-%A_%a.out

# Set the parameter combination to use and generate names of R scripts and log file
Rscript=fit_lnorm_ricker_mod_s4t.R
LogFile_pref=logfiles/ricker_dist_mfit_

# Change to the relevant working directory
cd /project/modelscape/analyses/sponges/Analyses

# Load modules
module load arcc/1.0  gcc/12.2.0  r/4.2.2

# create arrays of arguments
job=()
for i in ${SLURM_ARRAY_TASK_ID[@]}; do
  job+=$((15000 + $i))
done

Rscript --vanilla $Rscript "Data/terrestrial_sim_data/lnorm_ricker/dist_sims_s4t/sims_separate/lnorm_ricker_dist_freq_x_nsp_x_int_${job}.rds" ${job} 21 40 "Data/terrestrial_sim_data/lnorm_ricker/dist_sims_s4t/mfits/mfit_results" dist > "$LogFile_pref${job}"



