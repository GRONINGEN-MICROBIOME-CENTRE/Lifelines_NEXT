#!/bin/bash
#SBATCH --job-name=Antidefensefinder 
#SBATCH --output=Antidefensefinder.out
#SBATCH --mem=8gb
#SBATCH --time=02:59:00
#SBATCH --cpus-per-task=24
#SBATCH --open-mode=truncate
#SBATCH --partition=regular

viral_proteins=$1 #path to FASTA file with predicted viral proteins

echo -e '\n---- RUNNING AntiDefenseFinder ----'

# Clean environment, load modules 
ml purge; ml Python; ml list
source /scratch/hb-llnext/conda_envs/defensefinder_venv/bin/activate; python -m pip list

# Run DefenseFinder
defense-finder run \
    --db-type gembase \
    --antidefensefinder-only \
    $viral_proteins \
    -o ADF_RESULTS

# Remove intermediate files
rm *idx
