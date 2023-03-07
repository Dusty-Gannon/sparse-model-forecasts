#!/bin/bash -l

#SBATCH --account=modelscape
#SBATCH --nodes=1
#SBATCH --time=00:30:00
#SBATCH --ntasks-per-node=4
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
module load arcc/1.0 gcc/12.2.0 r/4.2.2

echo '----------------------' > slurmlogs/slurm.out 
echo '----------------------' > slurmlogs/slurm.err
Rscript --vanilla $Rscript 2 > $LogFile
#R CMD BATCH --no-save --no-restore --args a=2 < $Rscript $LogFile
