#!/bin/bash

#SBATCH --account=modelscape
#SBATCH --nodes=1
#SBATCH --time=24:00:00
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=3
#SBATCH --mem-per-cpu=8G
#SBATCH --mail-type=ALL
### please enter your own email address below in order to track the results
#SBATCH --mail-user=apatte12@uwyo.edu
### enter any job name that you prefer
#SBATCH --job-name=AICvsSTAN


module load arcc/1.0 gcc/12.2.0 r/4.2.2

cd /project/modelscape/analyses/sponges

### please enter your own arguments below and make a new filename for the output
### Summary of arguments to be provided:
### Number of trials, time steps (n), number covars (K), number strong covars, probability of cycling, trend fraction, frequency, sigma (sd of noise), correlated (logical), rateCorr (smaller gives higher correlation values), nameID for naming files

Rscript Simulations/AICvsStanBatch.R 10 100 40 5 0.2 0.2 1 0.5 F 2 ACP1 > Data/AICvsStan_data/outputCompareACP1.txt