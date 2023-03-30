#!/bin/bash -l

#SBATCH --account=modelscape
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=21
#SBATCH --cpus-per-task=1
#SBATCH --mem=20G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dgannon@uwyo.edu
#SBATCH --array=1-1010%10
#SBATCH --out=./logfiles/slurm-%A_%a.out

# Set the parameter combination to use and generate names of R scripts and log file
Rscript=fit_lnorm_ricker_mods_parallel.R
LogFile_pref=logfiles/ricker_thin2_

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

Rscript --vanilla $Rscript lnorm_ricker_thin_freq_x_nsp_ordered_S5_s55.rds 51 150 ${start} ${stop} "freq_x_nsp_results/freq_x_nsp_thintests_${start}_${stop}.rds" > "$LogFile_pref${start}_${stop}"



