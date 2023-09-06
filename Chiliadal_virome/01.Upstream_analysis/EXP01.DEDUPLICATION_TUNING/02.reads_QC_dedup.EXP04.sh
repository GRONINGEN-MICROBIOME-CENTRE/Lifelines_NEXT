#!/bin/bash
#SBATCH --job-name=reads_QC
#SBATCH --error=reads_QC.err
#SBATCH --output=reads_QC.out
#SBATCH --mem=32gb
#SBATCH --time=4:59:00
#SBATCH --cpus-per-task=4
#SBATCH --open-mode=truncate

SAMPLE_ID=$1
echo "SAMPLE_ID=${SAMPLE_ID}"

# --- LOAD MODULES --- 
module purge
module load BBMap

#### READ DEDUPLICATION: interleaved reads (for paired reads and for unmatched reads)
mkdir -p ../${SAMPLE_ID}/clean_reads/EXP04

#### INTERLEAVING PAIRED READS
reformat.sh \
	in1=../${SAMPLE_ID}/filtering_data/${SAMPLE_ID}_ECC_paired_1.fastq \
	in2=../${SAMPLE_ID}/filtering_data/${SAMPLE_ID}_ECC_paired_2.fastq \
	out=../${SAMPLE_ID}/filtering_data/${SAMPLE_ID}_ECC_paired_ITR.fastq

cat ../${SAMPLE_ID}/filtering_data/${SAMPLE_ID}_ECC_unmatched_1.fastq \
	../${SAMPLE_ID}/filtering_data/${SAMPLE_ID}_ECC_unmatched_2.fastq > \
	../${SAMPLE_ID}/filtering_data/${SAMPLE_ID}_ECC_unmatched.fastq

#### DEDUPLICATION
clumpify.sh \
        in=../${SAMPLE_ID}/filtering_data/${SAMPLE_ID}_ECC_paired_ITR.fastq \
        out=../${SAMPLE_ID}/clean_reads/EXP04/${SAMPLE_ID}_04_dedup_paired_ITR.fastq \
        dedupe=t \
        subs=0 \
        passes=2 \
        deletetemp=t 2>&1 \
       threads=${SLURM_CPUS_PER_TASK} \
       -Xmx$((${SLURM_MEM_PER_NODE} / 1000))g | tee -a ../${SAMPLE_ID}/${SAMPLE_ID}_04_dedup.log

clumpify.sh \
        in=../${SAMPLE_ID}/filtering_data/${SAMPLE_ID}_ECC_unmatched.fastq \
        out=../${SAMPLE_ID}/clean_reads/EXP04/${SAMPLE_ID}_04_dedup_unmatched.fastq \
        dedupe=t \
        subs=0 \
        passes=2 \
        deletetemp=t 2>&1 \
       threads=${SLURM_CPUS_PER_TASK} \
       -Xmx$((${SLURM_MEM_PER_NODE} / 1000))g | tee -a ../${SAMPLE_ID}/${SAMPLE_ID}_04_dedup.log

module load FastQC

fastqc ../${SAMPLE_ID}/clean_reads/EXP04/${SAMPLE_ID}_04_dedup_paired_ITR.fastq

mv ../${SAMPLE_ID}/clean_reads/EXP04/${SAMPLE_ID}_*fastqc* ../FastQC_reports/05_DEDUP_EXP04

module list

module purge

