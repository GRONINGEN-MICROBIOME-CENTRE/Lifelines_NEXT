#!/bin/bash
#SBATCH --job-name=VIR_DB
#SBATCH --output=./out/10.vin/PD_vir_index_N_vir_genes.out
#SBATCH --mem=32gb
#SBATCH --time=02:29:00
#SBATCH --cpus-per-task=8
#SBATCH --open-mode=truncate

VIR_DB=$1

echo "Virus database: ${VIR_DB}"

VIR_DB_DIR=$(dirname ${VIR_DB})

echo "Working directory: ${VIR_DB_DIR}"

# --- LOAD MODULES ---
module purge
module load Bowtie2
module load SAMtools
module list

# --- BUILDING BOWTIE2 INDEX FOR DEREPLICATED VIRUS CONTIGS ---
mkdir -p ${VIR_DB_DIR}/VIR_DB_index

bowtie2-build \
	${VIR_DB} \
	${VIR_DB_DIR}/VIR_DB_index/NEXT_VIR_DB \
	--large-index \
	--threads ${SLURM_CPUS_PER_TASK}

# --- GENERATING BED FILE ---
samtools faidx \
	${VIR_DB}

awk 'BEGIN {FS="\t"}; {print $1 FS "0" FS $2}' \
	${VIR_DB}.fai \
	> ${VIR_DB_DIR}/VIR_DB_index/NEXT_VIR_DB.bed

# --- PREPARING FOLDERS FOR READ MAPPING ---
mkdir -p ../VIR_DB/mapping/VLP_N_vir_genes/alignment_log
mkdir -p ../VIR_DB/mapping/VLP_N_vir_genes/coverage

module purge
