#!/bin/bash
#SBATCH --job-name=VIR_DB_NCP
#SBATCH --error=./err/08.vin/PD_index_NCP_decontam_99.err
#SBATCH --output=./out/08.vin/PD_index_NCP_decontam_99.out
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
mkdir -p ../VIR_DB/virus_contigs/der95_index_NCP_decontam_99

bowtie2-build \
	../VIR_DB/virus_contigs/NCP_vOTU_representatives_der95_wo_nc99.fasta \
	../VIR_DB/virus_contigs/der95_index_NCP_decontam_99/der95_wo_nc99_NCP \
	--large-index \
	--threads ${SLURM_CPUS_PER_TASK}

# --- GENERATING BED FILE ---
samtools faidx \
	../VIR_DB/virus_contigs/NCP_vOTU_representatives_der95_wo_nc99.fasta

awk 'BEGIN {FS="\t"}; {print $1 FS "0" FS $2}' \
	../VIR_DB/virus_contigs/NCP_vOTU_representatives_der95_wo_nc99.fasta.fai \
	> ../VIR_DB/virus_contigs/der95_index_NCP_decontam_99/der95_wo_nc99_NCP.bed

# --- PREPARING FOLDERS FOR READ MAPPING ---
mkdir -p ../VIR_DB/mapping/VLP_to_wo_nc99_der95_NCP/alignment_log_NCP
mkdir -p ../VIR_DB/mapping/VLP_to_wo_nc99_der95_NCP/coverage_NCP

module purge
