#!/bin/bash
#SBATCH --job-name=Anvio_merge_profiles 
#SBATCH --output=Anvio_merge_profiles.out
#SBATCH --mem=200gb
#SBATCH --time=05:59:00
#SBATCH --cpus-per-task=24
#SBATCH --open-mode=truncate
#SBATCH --partition=himem

vOTU_DB=$1 #path to the vOTU DB
PROFILES=$2 #path to the directory with anvio per-sample profiles

echo -e '\n---- RUNNING anvio merge ----'

# Clean environment, load modules 
ml purge; ml Anaconda3; ml list
conda activate /scratch/hb-llnext/conda_envs/anvio-8; conda list

# Run anvi'o
anvi-merge \
    -c $vOTU_DB \
    ${PROFILES}/*/PROFILE.db \
    -o merged_profile \
    --sample-name merged_profile
