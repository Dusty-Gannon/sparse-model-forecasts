#!/bin/bash

#SBATCH --account=modelscape
#SBATCH --nodes=1
#SBATCH --time=00:30:00
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8G
#SBATCH --mail-type=ALL
### please enter your own email address below in order to track the results
#SBATCH --mail-user=dgannon@uwyo.edu
### enter any job name that you prefer
#SBATCH --job-name=combineAICvsSTAN

# Combines the per-array-task .rds output from run_AICvsSTAN.sh with the
# config file into a single tidy .csv. Submit this AFTER the array job in
# run_AICvsSTAN.sh has finished, e.g. from the login node:
#
#   arrayjob=$(sbatch --parsable Simulations/run_AICvsSTAN.sh)
#   sbatch --dependency=afterany:${arrayjob} Simulations/run_combine_AICvsStan.sh
#
# Optional positional args (all have defaults, see combine_AICvsStan_output.R):
#   $1 = data_dir    directory holding the .rds output files (default Data/AICvsStan_data)
#   $2 = config_file config file used to launch the array job (default Simulations/AICvsStanConfig.txt)
#   $3 = out_file    path to write the combined .csv (default Simulations/AICvsStanResults.csv)

module load arcc/1.0 gcc/14.2.0 r/4.4.0 r-mass/7.3-59

cd /project/modelscape/analyses/sparse-model-forecasts

data_dir=${1:-Data/AICvsStan_data}
config_file=${2:-Simulations/AICvsStanConfig.txt}
out_file=${3:-Simulations/AICvsStanResults.csv}

Rscript Simulations/combine_AICvsStan_output.R ${data_dir} ${config_file} ${out_file}
