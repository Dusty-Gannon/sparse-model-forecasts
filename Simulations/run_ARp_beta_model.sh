#!/bin/bash -l

#SBATCH --account=modelscape
#SBATCH --nodes=1
#SBATCH --time=00:30:00
#SBATCH --ntasks-per-node=21
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=acarte26@beartooth.arcc.uwyo.edu
#SBATCH --job-name=ARp_sim1
#SBATCH --mem=10G
#SBATCH -o slurmlogs/slurm.out
#SBATCH -e slurmlogs/slurm.err

# Set the parameter combination to use and generate names of R scripts and log file
Rscript=AR-p_beta_p_model_sims.R
LogFile=ARpmod_sims.log

# Change to the relevant working directory
cd /project/modelscape/analyses/sponges/Simulations

# Load R and MPI
module load gcc/11.2.0 r/4.2.2s

echo '----------------------' > slurm.out 
echo '----------------------' > slurm.err
#Rscript --vanilla $Rscript 2 > $LogFile
R CMD BATCH --no-save --no-restore -2 $Rscript $LogFile
