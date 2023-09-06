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

#### READ DEDUPLICATION: interleaved reads concatenated with unmatched
mkdir -p ../${SAMPLE_ID}/clean_reads/EXP06

#### DEDUPLICATION
clumpify.sh \
        in=../${SAMPLE_ID}/filtering_data/${SAMPLE_ID}_ECC_cat.fastq \
        out=../${SAMPLE_ID}/clean_reads/EXP06/${SAMPLE_ID}_06_dedup_cat.fastq \
        dedupe=t \
        subs=0 \
        passes=2 \
        deletetemp=t 2>&1 \
       threads=${SLURM_CPUS_PER_TASK} \
       -Xmx$((${SLURM_MEM_PER_NODE} / 1000))g | tee -a ../${SAMPLE_ID}/${SAMPLE_ID}_06_dedup.log

# REPAIRING READS
repair.sh \
        in=../${SAMPLE_ID}/clean_reads/EXP06/${SAMPLE_ID}_06_dedup_cat.fastq \
        out1=../${SAMPLE_ID}/clean_reads/EXP06/${SAMPLE_ID}_06_dedup_repaired_1.fastq \
        out2=../${SAMPLE_ID}/clean_reads/EXP06/${SAMPLE_ID}_06_dedup_repaired_2.fastq \
        outs=../${SAMPLE_ID}/clean_reads/EXP06/${SAMPLE_ID}_06_dedup_singletons.fastq \
        repair

module load FastQC

fastqc ../${SAMPLE_ID}/clean_reads/EXP06/${SAMPLE_ID}_06_dedup_cat.fastq

mv ../${SAMPLE_ID}/clean_reads/EXP06/${SAMPLE_ID}_*fastqc* ../FastQC_reports/05_DEDUP_EXP06

module list

module purge

