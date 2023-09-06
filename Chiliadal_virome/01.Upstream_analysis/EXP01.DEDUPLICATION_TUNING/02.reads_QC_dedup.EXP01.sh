#!/bin/bash
#SBATCH --job-name=reads_QC
#SBATCH --error=reads_QC.err
#SBATCH --output=reads_QC.out
#SBATCH --mem=48gb
#SBATCH --time=4:59:00
#SBATCH --cpus-per-task=4
#SBATCH --open-mode=truncate

SAMPLE_ID=$1
echo "SAMPLE_ID=${SAMPLE_ID}"

# --- LOAD MODULES --- 
module purge
module load BBMap

#### READ DEDUPLICATION: paired mode
mkdir -p ../${SAMPLE_ID}/clean_reads/EXP01

clumpify.sh \
        in=../${SAMPLE_ID}/filtering_data/${SAMPLE_ID}_ECC_paired_1.fastq \
        in2=../${SAMPLE_ID}/filtering_data/${SAMPLE_ID}_ECC_paired_2.fastq \
        out=../${SAMPLE_ID}/clean_reads/EXP01/${SAMPLE_ID}_01_dedup_paired_1.fastq \
        out2=../${SAMPLE_ID}/clean_reads/EXP01/${SAMPLE_ID}_01_dedup_paired_2.fastq \
        dedupe=t \
        subs=0 \
        passes=2 \
        deletetemp=t 2>&1 \
       threads=${SLURM_CPUS_PER_TASK} \
       -Xmx48g | tee -a ../${SAMPLE_ID}/${SAMPLE_ID}_01_dedup.log

clumpify.sh \
        in=../${SAMPLE_ID}/filtering_data/${SAMPLE_ID}_ECC_unmatched_1.fastq \
        out=../${SAMPLE_ID}/clean_reads/EXP01/${SAMPLE_ID}_01_dedup_unmatched_1.fastq \
        dedupe=t \
        subs=0 \
        passes=2 \
        deletetemp=t 2>&1 \
       threads=${SLURM_CPUS_PER_TASK} \
       -Xmx48g | tee -a ../${SAMPLE_ID}/${SAMPLE_ID}_01_dedup.log

clumpify.sh \
        in=../${SAMPLE_ID}/filtering_data/${SAMPLE_ID}_ECC_unmatched_2.fastq \
        out=../${SAMPLE_ID}/clean_reads/EXP01/${SAMPLE_ID}_01_dedup_unmatched_2.fastq \
        dedupe=t \
        subs=0 \
        passes=2 \
        deletetemp=t 2>&1 \
       threads=${SLURM_CPUS_PER_TASK} \
       -Xmx48g | tee -a ../${SAMPLE_ID}/${SAMPLE_ID}_01_dedup.log

##### additional formatting for assembly
cat ../${SAMPLE_ID}/clean_reads/EXP01/${SAMPLE_ID}_01_dedup_unmatched_1.fastq \
	../${SAMPLE_ID}/clean_reads/EXP01/${SAMPLE_ID}_01_dedup_unmatched_2.fastq > \
	../${SAMPLE_ID}/clean_reads/EXP01/${SAMPLE_ID}_01_dedup_unmatched.fastq

module load FastQC

fastqc ../${SAMPLE_ID}/clean_reads/EXP01/${SAMPLE_ID}_01_dedup_paired_1.fastq
fastqc ../${SAMPLE_ID}/clean_reads/EXP01/${SAMPLE_ID}_01_dedup_paired_2.fastq
fastqc ../${SAMPLE_ID}/clean_reads/EXP01/${SAMPLE_ID}_01_dedup_unmatched_1.fastq
fastqc ../${SAMPLE_ID}/clean_reads/EXP01/${SAMPLE_ID}_01_dedup_unmatched_2.fastq

mv ../${SAMPLE_ID}/clean_reads/EXP01/${SAMPLE_ID}_*fastqc* ../FastQC_reports/05_DEDUP_EXP01

module list

module purge

