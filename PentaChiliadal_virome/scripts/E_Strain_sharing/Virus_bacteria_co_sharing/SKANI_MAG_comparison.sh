#!/bin/bash
#SBATCH --job-name=SKANI_MAG_comparison 
#SBATCH --output=SKANI_MAG_comparison.out
#SBATCH --mem=8gb
#SBATCH --time=00:19:00
#SBATCH --cpus-per-task=8
#SBATCH --open-mode=truncate
#SBATCH --partition=regular

MAG_dir=$1 #path to directory with MAGs to be compared

echo -e '\n---- RUNNING SKANI ----'

# Clean environment, load modules 
ml purge; ml Anaconda3; ml list
conda activate /scratch/hb-llnext/conda_envs/SKANI; conda list

mkdir -p SKANI_RESULTS

# Run SKANI comparison
skani triangle \
    ${MAG_dir}/* \
    --full-matrix \
    -t ${SLURM_CPUS_PER_TASK} \
    > SKANI_RESULTS/skani_matrix.txt

