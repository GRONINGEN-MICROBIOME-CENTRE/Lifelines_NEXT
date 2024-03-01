#!/bin/bash
#SBATCH --job-name=VIR_DB
#SBATCH --error=./err/10.vin/PD_w_neg_index.err
#SBATCH --output=./out/10.vin/PD_w_neg_index.out
#SBATCH --mem=32gb
#SBATCH --time=01:59:00
#SBATCH --cpus-per-task=8
#SBATCH --open-mode=truncate

# --- LOAD MODULES ---
module purge
module load Bowtie2
module load SAMtools
module list

# --- BUILDING BOWTIE2 INDEX FOR DEREPLICATED VIRUS CONTIGS ---
mkdir -p ../VIR_DB/virus_contigs/w_neg_der95_index

bowtie2-build \
	../VIR_DB/virus_contigs/NEXT_vOTU_representatives_w_neg_der95.fasta \
	../VIR_DB/virus_contigs/w_neg_der95_index/w_neg_der95 \
	--large-index \
	--threads ${SLURM_CPUS_PER_TASK}

# --- GENERATING BED FILE ---
samtools faidx \
	../VIR_DB/virus_contigs/NEXT_vOTU_representatives_w_neg_der95.fasta

awk 'BEGIN {FS="\t"}; {print $1 FS "0" FS $2}' \
	../VIR_DB/virus_contigs/NEXT_vOTU_representatives_w_neg_der95.fasta.fai \
	> ../VIR_DB/virus_contigs/w_neg_der95_index/w_neg_der95.bed

# --- PREPARING FOLDERS FOR READ MAPPING ---
mkdir -p ../VIR_DB/mapping/VLP_to_w_neg_der95/alignment_log
mkdir -p ../VIR_DB/mapping/VLP_to_w_neg_der95/coverage

module purge
