#!/bin/bash

#SBATCH --job-name=PP_1
#SBATCH --output=./out/PP/PP_1_%A.out
#SBATCH --error=./err/PP/PP_1_%A.err
#SBATCH --time=00:59:00
#SBATCH --mem=20GB
#SBATCH --cpus-per-task=1
#SBATCH --open-mode=truncate
#SBATCH --partition=regular

VIRAL_GENOMES=$1 #path to FASTA file with viral sequences
ETOF=$2 #path to vOTU metadata
OUTPUT_DIR=$3 #path to output 
SEQS_PER_BATCH=$4 #number of sequences in one batch

echo "Viral genomes file: $(realpath "$VIRAL_GENOMES")"
echo "vOTU metadata file: $(realpath "$ETOF")"
echo "Output directory: $(realpath "$OUTPUT_DIR")"
echo "Number of sequences per batch: ${SEQS_PER_BATCH}"

# Filter by quality first
mkdir -p "${OUTPUT_DIR}"

awk -F'\t' 'NR>1 && $9=="High-quality" {print $1}' "$ETOF" > "${OUTPUT_DIR}/High_quality_vOTUs.txt"

module purge
module load SeqKit

module list

seqkit grep -f "${OUTPUT_DIR}/High_quality_vOTUs.txt" "$VIRAL_GENOMES" > "${OUTPUT_DIR}/High_quality_vOTUs.fasta"

module purge
module load Anaconda3/2022.05
conda activate iphop133_env

module list
conda list

mkdir -p "${OUTPUT_DIR}/splitted_fasta"

iphop split \
        --input_file "${OUTPUT_DIR}/High_quality_vOTUs.fasta" \
        --split_dir "${OUTPUT_DIR}/splitted_fasta" \
        --n_seq "$SEQS_PER_BATCH"

ls "${OUTPUT_DIR}/splitted_fasta" > "${OUTPUT_DIR}/batch_list.txt"
sed -i 's/.fna//g' "${OUTPUT_DIR}/batch_list.txt"

conda deactivate
module purge

N_BATCHES=$(wc -l < "${OUTPUT_DIR}/batch_list.txt")
N_BATCHES=$(echo "$N_BATCHES" | xargs)  # trim spaces

echo "Number of batches: ${N_BATCHES}"

mkdir -p "${OUTPUT_DIR}/PP_PNT"
# mkdir -p "${OUTPUT_DIR}/PP_PGV"

sbatch --array=1-"$N_BATCHES" PP_PNT_2.sh "${OUTPUT_DIR}/batch_list.txt" "$OUTPUT_DIR"
# sbatch --array=1-"$N_BATCHES" PP_PGV_3.sh "${OUTPUT_DIR}/batch_list.txt" "$OUTPUT_DIR"
