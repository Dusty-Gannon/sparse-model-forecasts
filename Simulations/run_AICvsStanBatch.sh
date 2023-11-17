#!/bin/bash

#SBATCH --account=modelscape
#SBATCH --nodes=1
#SBATCH --time=24:00:00
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --mail-type=ALL
### please enter your own email address below in order to track the results
#SBATCH --mail-user=apatte12@uwyo.edu
### enter any job name that you prefer
#SBATCH --job-name=AICvsSTANfull
#SBATCH --array=1-4


module load arcc/1.0 gcc/12.2.0 r/4.2.2

cd /project/modelscape/analyses/sponges

config=Simulations/AICvsStanConfig.txt

numTrials=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)
n=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $3}' $config)
K=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $4}' $config)
num_strong=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $5}' $config)
prob_cycle=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $6}' $config)
trend_fraction=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $7}' $config)
freq=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $8}' $config)
sigma=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $9}' $config)
probWeakCorr=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $10}' $config)
probStrongCorr=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $11}' $config)
corrLevel=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $12}' $config)
nameID=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $13}' $config)


Rscript Simulations/AICvsStanBatch.R ${numTrials} ${n} ${K} ${num_strong} ${prob_cycle} ${trend_fraction} ${freq} ${sigma} ${probWeakCorr} ${probStrongCorr} ${corrLevel} ${nameID} > Data/AICvsStan_data/outputCompareACP1.txt

