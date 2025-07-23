#!/bin/bash
#SBATCH --job-name=protein_search
#SBATCH --output=./out/14.fan/debug.out
#SBATCH --mem=2gb
#SBATCH --time=00:09:00
#SBATCH --cpus-per-task=1
#SBATCH --open-mode=truncate

INPUT=$1
OUTPUT_DIR=$2

echo "Viral genomes file: $(realpath "${INPUT}")"
echo "Output directory: $(realpath "${OUTPUT_DIR}")"

INPUT="${INPUT%.fasta}"

echo ${INPUT}

# --- LOAD MODULES --- 
module purge
module load prodigal-gv/2.11.0-GCCcore-12.2.0
module load Python/3.11.3-GCCcore-12.3.0
module list

# --- PREDICTING ORFs ---
echo "> Running parallel prodigal-gv"

# parallel-prodigal-gv.py can be found at: https://raw.githubusercontent.com/apcamargo/prodigal-gv/master/parallel-prodigal-gv.py

python /scratch/p282752/ANALYSIS_CHILIADAL/scripts/parallel-prodigal-gv.py \
        -t ${SLURM_CPUS_PER_TASK} \
        -q \
        -i ${INPUT}.fasta \
        -a ${INPUT}.faa \
        -o ${OUTPUT_DIR}/prodigal.out

# --- CLEAN ENV ---
module purge

# --- LOAD MODULES --- 
module load DIAMOND
module list

# --- MAKE DIAMOND DB ---
diamond makedb --in ${INPUT}.faa --db ${OUTPUT_DIR}/viral_proteins --threads 10

# --- Perform all-vs-all BLASTP ---
diamond blastp \
	--query ${INPUT}.faa \
	--db ${OUTPUT_DIR}/viral_proteins \
	--out ${OUTPUT_DIR}/blastp.tsv \
	--outfmt 6 \
	--evalue 1e-5 \
	--max-target-seqs 10000 \
	--query-cover 50 \
	--subject-cover 50

# --- CLEAN ENV --- 
module purge

# --- LAUNCHING CLUSTERING --- 
echo "> Launcing AAI calculation and clustering"

sbatch 14b.fa_calc_aai_mcl.sh ${INPUT}.faa ${OUTPUT_DIR}

