#!/bin/bash

#SBATCH --job-name=DSP_PDC
#SBATCH --output=./out/DSP/DSP_PDC_%A.out
#SBATCH --error=./err/DSP/DSP_PDC_%A.err
#SBATCH --time=12:59:00
#SBATCH --mem=30GB
#SBATCH --cpus-per-task=8
#SBATCH --open-mode=truncate
#SBATCH --partition=regular

VIRAL_PROTEINS_FAA=$1 #path to FAA file with viral sequences
VIRAL_PROTEINS_GFF=$2 #path to GFF file with viral sequences
OUTPUT_DIR=$3 #path to output directory

echo "Viral .faa protein file: $(realpath "$VIRAL_PROTEINS_FAA")"
echo "Viral .gff protein file: $(realpath "$VIRAL_PROTEINS_GFF")"
echo "Output directory: $(realpath "$OUTPUT_DIR")"

module purge
module load Anaconda3/2024.02-1
conda activate /scratch/hb-llnext/conda_envs/padloc_env

module list
conda list

padloc \
	--faa "$VIRAL_PROTEINS_FAA" \
	--gff "$VIRAL_PROTEINS_GFF" \
	--outdir "${OUTPUT_DIR}" \
	--cpu 8

conda deactivate
module purge
