#!/bin/bash -l

#SBATCH --account=modelscape
#SBATCH --time=05:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=3
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dgannon@uwyo.edu
#SBATCH --array=1-10

# Set the parameter combination to use and generate names of R scripts and log file
Rscript=pois_latAR1_FHS.R
LogFile=pois_latAR1_ntest.log

# Change to the relevant working directory
cd /project/modelscape/analyses/sponges/Model_evals

# Load R and MPI
module load gcc/7.3.0 swset/2018.05 r/3.5.3 r-rstan/2.18.2-py27

# create array of arguments
args=()
for i in ${SLURM_ARRAY_TASK_ID[@]}; do
  args+=$(( 50*i ))
done


Rscript $Rscript ${args} > $LogFile

