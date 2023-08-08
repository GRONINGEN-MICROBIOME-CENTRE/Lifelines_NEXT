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
#bbduk.sh \
#	in=../MVS_raw_reads/${SAMPLE_ID}/raw_reads/${SAMPLE_ID}_1.fastq.gz \
#	out=../MVS_raw_reads/${SAMPLE_ID}/raw_reads/${SAMPLE_ID}_1.NoTail.fastq.gz \
#	ftr2=11 tpe tbo 2>&1 \
#        threads=${SLURM_CPUS_PER_TASK} \
#	-Xmx32g | tee ../MVS_raw_reads/${SAMPLE_ID}/raw_reads/${SAMPLE_ID}_bbduk.log
#bbduk.sh \
#        in=../MVS_raw_reads/${SAMPLE_ID}/raw_reads/${SAMPLE_ID}_2.fastq.gz \
#        out=../MVS_raw_reads/${SAMPLE_ID}/raw_reads/${SAMPLE_ID}_2.NoTail.fastq.gz \
#        ftl=11 tpe tbo 2>&1 \
#        threads=${SLURM_CPUS_PER_TASK} \
#        -Xmx32g | tee ../MVS_raw_reads/${SAMPLE_ID}/raw_reads/${SAMPLE_ID}_bbduk.log

#### TRIMMING ADAPTERS
#bbduk.sh \
#	in1=../MVS_raw_reads/${SAMPLE_ID}/raw_reads/${SAMPLE_ID}_1.NoTail.fastq.gz \
#	in2=../MVS_raw_reads/${SAMPLE_ID}/raw_reads/${SAMPLE_ID}_2.NoTail.fastq.gz \
#	out1=../MVS_raw_reads/${SAMPLE_ID}/raw_reads/${SAMPLE_ID}_AdaptTr_QualTr_1.fastq.gz \
#    	out2=../MVS_raw_reads/${SAMPLE_ID}/raw_reads/${SAMPLE_ID}_AdaptTr_QualTr_2.fastq.gz \
#	ref=/scratch/p282752/Data_for_HiC/adapters_UPD_IDT.fa \
#	ktrim=r k=23 mink=11 hdist=1 tpe tbo 2>&1 \
#	threads=${SLURM_CPUS_PER_TASK} \
#	-Xmx32g | tee ../MVS_raw_reads/${SAMPLE_ID}/raw_reads/${SAMPLE_ID}_bbduk.log

#### FILTERING HUMAN READS & QUALITY TRIMMING 
module load Anaconda3/2022.05
conda activate /scratch/hb-tifn/condas/conda_biobakery3/

kneaddata \
        --input ../MVS_raw_reads/${SAMPLE_ID}/raw_reads/${SAMPLE_ID}_AdaptTr_QualTr_1.fastq.gz \
        --input ../MVS_raw_reads/${SAMPLE_ID}/raw_reads/${SAMPLE_ID}_AdaptTr_QualTr_2.fastq.gz \
        --threads 4 \
        --processes 4 \
        --output-prefix ${SAMPLE_ID}_kneaddata \
        --output ../MVS_raw_reads/${SAMPLE_ID}/raw_reads/ \
        --log ../MVS_raw_reads/${SAMPLE_ID}/raw_reads/kneaddata.log \
        -db /scratch/hb-tifn/DBs/human_genomes/hg37dec_v0.1  \
        --trimmomatic /scratch/hb-tifn/condas/conda_biobakery4/share/trimmomatic-0.39-2/ \
        --run-trim-repetitive \
        --fastqc fastqc \
        --sequencer-source none \
        --trimmomatic-options "LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:50" \
        --bypass-trf \
        --reorder

module load FastQC
#fastqc ../MVS_raw_reads/${SAMPLE_ID}/raw_reads/${SAMPLE_ID}_1.NoTail.fastq.gz
#fastqc ../MVS_raw_reads/${SAMPLE_ID}/raw_reads/${SAMPLE_ID}_2.NoTail.fastq.gz
#mv ../MVS_raw_reads/${SAMPLE_ID}/raw_reads/${SAMPLE_ID}_*fastqc* ../MVS_raw_reads/FastQC_reports/NoTail/
#fastqc ../MVS_raw_reads/${SAMPLE_ID}/raw_reads/${SAMPLE_ID}_AdaptTr_QualTr_1.fastq.gz
#fastqc ../MVS_raw_reads/${SAMPLE_ID}/raw_reads/${SAMPLE_ID}_AdaptTr_QualTr_2.fastq.gz
#mv ../MVS_raw_reads/${SAMPLE_ID}/raw_reads/${SAMPLE_ID}_*AdaptTr_QualTr*fastqc* ../MVS_raw_reads/FastQC_reports/AdaptTr_QualTr/
fastqc ../MVS_raw_reads/${SAMPLE_ID}/raw_reads/${SAMPLE_ID}_kneaddata_paired_1.fastq
fastqc ../MVS_raw_reads/${SAMPLE_ID}/raw_reads/${SAMPLE_ID}_kneaddata_paired_2.fastq
mv ../MVS_raw_reads/${SAMPLE_ID}/raw_reads/${SAMPLE_ID}_kneaddata*fastqc* ../MVS_raw_reads/FastQC_reports/Kneaddata/


module list

module purge

