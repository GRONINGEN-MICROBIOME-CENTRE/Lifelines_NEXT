#!/bin/bash
#SBATCH --job-name=Structure_prediction
#SBATCH --output=Structure_prediction_%A_%a.out
#SBATCH --mem=100gb
#SBATCH --time=7-00:00 
#SBATCH --cpus-per-task=32
#SBATCH --export=NONE
#SBATCH --open-mode=truncate
#SBATCH --partition=regular

MSAs=$1 # directory with alignments (.a3m files) generated with HHBlits
list_a3m=$2 # list of A3M (enriched cluster files)
output_dir=$3 # output directory

cluster_ID=$(sed "${SLURM_ARRAY_TASK_ID}q;d" ${list_a3m} | cut -f1 -d '.')

echo -e '\n-------------------- WORKING WITH '${cluster_ID}' PROTEIN FAMILY --------------------'

# Activate conda environment
source $JGI_SCRATCH/afernandezpato/Tools/miniconda3/bin/activate \
    $JGI_SCRATCH/afernandezpato/Tools/localcolabfold/colabfold-conda; conda list

# We run colabfold_batch to perform the structural prediction
colabfold_batch \
    $MSAs/${cluster_ID}.a3m \
    ${output_dir}/$cluster_ID \
    --num-recycle=3 \
    --recycle-early-stop-tolerance 0.5 \
    --num-seeds 3

# Keep only top-ranked model by pLDDT for each PF
cd ${output_dir}/$cluster_ID
rm *png 
rm config.json 
rm cite.bibtex

