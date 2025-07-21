#!/bin/bash
#SBATCH --job-name=m_prep
#SBATCH --output=./out/00.mpr/prep_Penta_%A_%a.out
#SBATCH --mem=4gb
#SBATCH --time=00:29:00
#SBATCH --cpus-per-task=1
#SBATCH --open-mode=truncate

SAMPLE_LIST=$1

SAMPLE_DIR=$2

echo ${SAMPLE_LIST}
echo ${SAMPLE_DIR}

SAMPLE_ID=$(sed "${SLURM_ARRAY_TASK_ID}q;d" ${SAMPLE_LIST} | cut -d "_" -f1)

echo "SAMPLE_ID=${SAMPLE_ID}"

# RECONSTRUCTING USUAL DIR STRUCTURE TO RUN MGS SAMPLES THROUGH VLP VIROME DISCOVERY PIPELINE
mkdir -p ${SAMPLE_DIR}/${SAMPLE_ID}/01_sc_assembly

# --- Contig trimming ---
echo "> Trimming contigs to 1kbp"
perl filter_contigs.pl \
	1000 \
	${SAMPLE_DIR}/${SAMPLE_ID}_metaspades_contigs.fa \
	> ${SAMPLE_DIR}/${SAMPLE_ID}/01_sc_assembly/${SAMPLE_ID}_contigs.min1kbp.fasta

# --- Contig rename ---
sed -i 's/>NODE/>'${SAMPLE_ID}'_NODE/g' ${SAMPLE_DIR}/${SAMPLE_ID}/01_sc_assembly/${SAMPLE_ID}_contigs.min1kbp.fasta

