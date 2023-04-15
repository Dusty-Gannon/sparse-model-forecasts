#!/bin/bash -l

#SBATCH --account=modelscape
#SBATCH --time=00:20:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=21
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dgannon@uwyo.edu
#SBATCH --array=1488-2000%10
#SBATCH --out=./logfiles/slurm-%A_%a.out

# Set the parameter combination to use and generate names of R scripts and log file
Rscript=fit_lnorm_ricker_mods_parallel.R
LogFile_pref=logfiles/ricker_dist1_

# Change to the relevant working directory
cd /project/modelscape/analyses/sponges/Analyses

# Load modules
module load arcc/1.0  gcc/12.2.0  r/4.2.2

# create arrays of arguments
start=()
stop=()
for i in ${SLURM_ARRAY_TASK_ID[@]}; do
  start+=$(( 10*(i-1)+1 ))
  stop+=$(( 10*i ))
done

Rscript --vanilla $Rscript lnorm_ricker_dist_freq_x_nsp_x_int_round1.rds 51 250 ${start} ${stop} "disturb_results/round1/disturb_tests_${start}_${stop}.rds" dist > "$LogFile_pref${start}_${stop}"



