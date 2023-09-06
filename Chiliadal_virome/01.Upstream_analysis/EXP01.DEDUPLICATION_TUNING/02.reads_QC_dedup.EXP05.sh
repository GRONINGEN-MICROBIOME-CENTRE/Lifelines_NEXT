#!/bin/bash
#SBATCH --job-name=reads_QC
#SBATCH --error=reads_QC.err
#SBATCH --output=reads_QC.out
#SBATCH --mem=40gb
#SBATCH --time=4:59:00
#SBATCH --cpus-per-task=4
#SBATCH --open-mode=truncate

SAMPLE_ID=$1
echo "SAMPLE_ID=${SAMPLE_ID}"

# --- LOAD MODULES --- 
module purge
module load BBMap

#### READ DEDUPLICATION: interleaved reads concatenated with unmatched
mkdir -p ../${SAMPLE_ID}/clean_reads/EXP05

#### INTERLEAVING PAIRED READS
cat ../${SAMPLE_ID}/filtering_data/${SAMPLE_ID}_ECC_paired_ITR.fastq \
	../${SAMPLE_ID}/filtering_data/${SAMPLE_ID}_ECC_unmatched.fastq > \
	../${SAMPLE_ID}/filtering_data/${SAMPLE_ID}_ECC_ITR.fastq

#### DEDUPLICATION
clumpify.sh \
        in=../${SAMPLE_ID}/filtering_data/${SAMPLE_ID}_ECC_ITR.fastq \
        out=../${SAMPLE_ID}/clean_reads/EXP05/${SAMPLE_ID}_05_dedup_ITR.fastq \
        dedupe=t \
        subs=0 \
        passes=2 \
        deletetemp=t 2>&1 \
       threads=${SLURM_CPUS_PER_TASK} \
       -Xmx$((${SLURM_MEM_PER_NODE} / 1000))g | tee -a ../${SAMPLE_ID}/${SAMPLE_ID}_05_dedup.log


module load FastQC

fastqc ../${SAMPLE_ID}/clean_reads/EXP05/${SAMPLE_ID}_05_dedup_ITR.fastq

mv ../${SAMPLE_ID}/clean_reads/EXP05/${SAMPLE_ID}_*fastqc* ../FastQC_reports/05_DEDUP_EXP05

module list

module purge

