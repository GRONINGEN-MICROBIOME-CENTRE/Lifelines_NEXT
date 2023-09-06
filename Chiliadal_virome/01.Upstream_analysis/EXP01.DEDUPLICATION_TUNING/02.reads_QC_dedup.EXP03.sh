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

#### READ DEDUPLICATION: unpaired mode but with concatenated unmatched to respective R1 and R2
mkdir -p ../${SAMPLE_ID}/clean_reads/EXP03


clumpify.sh \
        in=../${SAMPLE_ID}/filtering_data/${SAMPLE_ID}_ECC_1.fastq \
        out=../${SAMPLE_ID}/clean_reads/EXP03/${SAMPLE_ID}_03_dedup_unpaired_1.fastq \
        dedupe=t \
        subs=0 \
        passes=2 \
        deletetemp=t 2>&1 \
       threads=${SLURM_CPUS_PER_TASK} \
       -Xmx48g | tee -a ../${SAMPLE_ID}/${SAMPLE_ID}_03_dedup.log

clumpify.sh \
        in=../${SAMPLE_ID}/filtering_data/${SAMPLE_ID}_ECC_2.fastq \
        out=../${SAMPLE_ID}/clean_reads/EXP03/${SAMPLE_ID}_03_dedup_unpaired_2.fastq \
        dedupe=t \
        subs=0 \
        passes=2 \
        deletetemp=t 2>&1 \
       threads=${SLURM_CPUS_PER_TASK} \
       -Xmx48g | tee -a ../${SAMPLE_ID}/${SAMPLE_ID}_03_dedup.log

# REPAIRING READS (those that were treated as unpaired output to deduplicating step)
repair.sh \
        in1=../${SAMPLE_ID}/clean_reads/EXP03/${SAMPLE_ID}_03_dedup_unpaired_1.fastq \
        in2=../${SAMPLE_ID}/clean_reads/EXP03/${SAMPLE_ID}_03_dedup_unpaired_2.fastq \
        out1=../${SAMPLE_ID}/clean_reads/EXP03/${SAMPLE_ID}_03_dedup_repaired_1.fastq \
        out2=../${SAMPLE_ID}/clean_reads/EXP03/${SAMPLE_ID}_03_dedup_repaired_2.fastq \
        outs=../${SAMPLE_ID}/clean_reads/EXP03/${SAMPLE_ID}_03_dedup_singletons.fastq \
        repair


module load FastQC

fastqc ../${SAMPLE_ID}/clean_reads/EXP03/${SAMPLE_ID}_03_dedup_unpaired_1.fastq
fastqc ../${SAMPLE_ID}/clean_reads/EXP03/${SAMPLE_ID}_03_dedup_unpaired_2.fastq

mv ../${SAMPLE_ID}/clean_reads/EXP03/${SAMPLE_ID}_*fastqc* ../FastQC_reports/05_DEDUP_EXP03

module list

module purge

