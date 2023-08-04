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

#### TRIMMING ADAPTASE-INTRODUCED TAILS 
bbduk.sh \
	in=../MVS_raw_reads/${SAMPLE_ID}/raw_reads/${SAMPLE_ID}_1.fastq.gz \
	out=../MVS_raw_reads/${SAMPLE_ID}/raw_reads/${SAMPLE_ID}_1.NoTail.fastq.gz \
	ftr2=11 tpe tbo 2>&1 \
        threads=${SLURM_CPUS_PER_TASK} \
	-Xmx32g | tee ../MVS_raw_reads/${SAMPLE_ID}/raw_reads/${SAMPLE_ID}_bbduk.log
bbduk.sh \
        in=../MVS_raw_reads/${SAMPLE_ID}/raw_reads/${SAMPLE_ID}_2.fastq.gz \
        out=../MVS_raw_reads/${SAMPLE_ID}/raw_reads/${SAMPLE_ID}_2.NoTail.fastq.gz \
        ftl=11 tpe tbo 2>&1 \
        threads=${SLURM_CPUS_PER_TASK} \
        -Xmx32g | tee ../MVS_raw_reads/${SAMPLE_ID}/raw_reads/${SAMPLE_ID}_bbduk.log

module load FastQC
fastqc ../MVS_raw_reads/${SAMPLE_ID}/raw_reads/${SAMPLE_ID}_1.NoTail.fastq.gz
fastqc ../MVS_raw_reads/${SAMPLE_ID}/raw_reads/${SAMPLE_ID}_2.NoTail.fastq.gz
mv ../MVS_raw_reads/${SAMPLE_ID}/raw_reads/${SAMPLE_ID}_*fastqc* ../MVS_raw_reads/FastQC_reports/NoTail/

module list

module purge

