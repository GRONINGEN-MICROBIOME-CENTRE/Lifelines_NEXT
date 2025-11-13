#!/bin/bash

#SBATCH --job-name=HP_1
#SBATCH --output=./out/HP/HP_1_%A.out
#SBATCH --error=./err/HP/HP_1_%A.err
#SBATCH --time=02:00:00
#SBATCH --mem=20GB
#SBATCH --cpus-per-task=1
#SBATCH --open-mode=truncate
#SBATCH --partition=regular

VIRAL_GENOMES=$1 #path to FASTA file with viral sequences
OUTPUT_DIR=$2 #path to output directory
SEQS_PER_BATCH=$3 #number of sequences in one batch 

echo "Viral genomes file: $(realpath "$VIRAL_GENOMES")"
echo "Output directory: $(realpath "$OUTPUT_DIR")"
echo "Number of sequences per batch: ${SEQS_PER_BATCH}"

module purge
module load Anaconda3/2022.05
conda activate iphop133_env

module list
conda list

mkdir -p "${OUTPUT_DIR}/splitted_fasta"

iphop split \
	--input_file "$VIRAL_GENOMES" \
	--split_dir "${OUTPUT_DIR}/splitted_fasta" \
	--n_seq "$SEQS_PER_BATCH"

ls "${OUTPUT_DIR}/splitted_fasta" > "${OUTPUT_DIR}/batch_list.txt"
sed -i 's/.fna//g' "${OUTPUT_DIR}/batch_list.txt"

conda deactivate
module purge

N_BATCHES=$(wc -l < "${OUTPUT_DIR}/batch_list.txt")
N_BATCHES=$(echo "$N_BATCHES" | xargs)  # trim spaces

echo "Number of batches: ${N_BATCHES}"

mkdir -p "${OUTPUT_DIR}/host_predictions"

sbatch --array=1-"$N_BATCHES" HP_2.sh "${OUTPUT_DIR}/batch_list.txt" "$OUTPUT_DIR"
sbatch --array=1-"$N_BATCHES" HP_3.sh "${OUTPUT_DIR}/batch_list.txt" "$OUTPUT_DIR"
