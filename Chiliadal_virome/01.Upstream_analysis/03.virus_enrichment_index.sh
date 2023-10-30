#!/bin/bash
#SBATCH --job-name=reads_QC
#SBATCH --error=reads_QC.err
#SBATCH --output=reads_QC.out
#SBATCH --mem=12gb
#SBATCH --time=24:59:00
#SBATCH --cpus-per-task=2
#SBATCH --open-mode=truncate

SAMPLE_ID=$1
echo "SAMPLE_ID=${SAMPLE_ID}"

# --- LOAD MODULES --- 
module purge
module load Python
module load Pysam
module load DIAMOND
module load Bowtie2
module load SAMtools
module load Biopython

# --- RUNNING ENRICHMENT CALCULATION ---
echo "> Running ViromeQC"
/scratch/p282752/tools/viromeqc/viromeQC.py \
        -i ../SAMPLES/${SAMPLE_ID}/clean_reads/${SAMPLE_ID}_dedup_paired_1.fastq.gz ../SAMPLES/${SAMPLE_ID}/clean_reads/${SAMPLE_ID}_dedup_paired_2.fastq.gz ../SAMPLES/${SAMPLE_ID}/clean_reads/${SAMPLE_ID}_dedup_unmatched_1.fastq.gz ../SAMPLES/${SAMPLE_ID}/clean_reads/${SAMPLE_ID}_dedup_unmatched_2.fastq.gz \
        --minlen 50 \
        --sample_name ${SAMPLE_ID} \
        -o ../SAMPLES/${SAMPLE_ID}/viromeqc_stat.txt

if [ -f ../SAMPLES/${SAMPLE_ID}/viromeqc_stat.txt ]; then
	echo "> ViromeQC has finished"
fi


module list

module purge

