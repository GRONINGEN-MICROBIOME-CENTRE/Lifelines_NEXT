#!/bin/bash

#SBATCH --job-name=DSP_DF
#SBATCH --output=./out/DSP/DSP_DF_%A.out
#SBATCH --error=./err/DSP/DSP_DF_%A.err
#SBATCH --time=00:59:00
#SBATCH --mem=10GB
#SBATCH --cpus-per-task=8
#SBATCH --open-mode=truncate
#SBATCH --partition=regular

VIRAL_PROTEINS=$1 #path to FASTA file with viral sequences
OUTPUT_DIR=$2 #path to output directory

echo "Viral protein file: $(realpath "$VIRAL_PROTEINS")"
echo "Output directory: $(realpath "$OUTPUT_DIR")"

module purge
module load Python/3.11.5-GCCcore-13.2.0
source /scratch/hb-llnext/conda_envs/defensefinder_venv/bin/activate

module list

defense-finder run "$VIRAL_PROTEINS" \
	-o "${OUTPUT_DIR}" \
	--db-type gembase \
	-a

deactivate
module purge
