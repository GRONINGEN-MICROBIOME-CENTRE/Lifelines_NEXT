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

#### TRIMMING ADAPTERS
bbduk.sh \
        in1=../${SAMPLE_ID}/raw_reads/${SAMPLE_ID}_1.fastq.gz \
        in2=../${SAMPLE_ID}/raw_reads/${SAMPLE_ID}_2.fastq.gz \
        out1=../${SAMPLE_ID}/filtering_data/${SAMPLE_ID}_AdaptTr_1.fastq.gz \
        out2=../${SAMPLE_ID}/filtering_data/${SAMPLE_ID}_AdaptTr_2.fastq.gz \
        ref=/scratch/p282752/Data_for_HiC/adapters_UPD_IDT.fa \
        ktrim=r k=23 mink=11 hdist=1 tpe tbo 2>&1 \
        threads=${SLURM_CPUS_PER_TASK} \
        -Xmx32g | tee -a ../${SAMPLE_ID}/${SAMPLE_ID}_bbduk.log


#### TRIMMING ADAPTASE-INTRODUCED TAILS 
bbduk.sh \
	in=../${SAMPLE_ID}/filtering_data/${SAMPLE_ID}_AdaptTr_1.fastq.gz \
	out=../${SAMPLE_ID}/filtering_data/${SAMPLE_ID}_NoTail_1.fastq.gz \
	ftr2=11 tpe tbo 2>&1 \
        threads=${SLURM_CPUS_PER_TASK} \
	-Xmx32g | tee -a ../${SAMPLE_ID}/${SAMPLE_ID}_bbduk.log
bbduk.sh \
        in=../${SAMPLE_ID}/filtering_data/${SAMPLE_ID}_AdaptTr_2.fastq.gz \
        out=../${SAMPLE_ID}/filtering_data/${SAMPLE_ID}_NoTail_2.fastq.gz \
        ftl=11 tpe tbo 2>&1 \
        threads=${SLURM_CPUS_PER_TASK} \
        -Xmx32g | tee -a ../${SAMPLE_ID}/${SAMPLE_ID}_bbduk.log

#### FILTERING HUMAN READS & LOW QUALITY READS 
module load Anaconda3/2022.05
conda activate /scratch/hb-tifn/condas/conda_biobakery3/

kneaddata \
	--input ../${SAMPLE_ID}/filtering_data/${SAMPLE_ID}_NoTail_1.fastq.gz \
        --input ../${SAMPLE_ID}/filtering_data/${SAMPLE_ID}_NoTail_2.fastq.gz \
        --threads ${SLURM_CPUS_PER_TASK} \
        --processes 4 \
        --output-prefix ${SAMPLE_ID}_kneaddata \
        --output ../${SAMPLE_ID}/filtering_data/ \
        --log ../${SAMPLE_ID}/${SAMPLE_ID}_kneaddata.log \
        -db /scratch/hb-tifn/DBs/human_genomes/hg37dec_v0.1  \
        --trimmomatic /scratch/hb-tifn/condas/conda_biobakery4/share/trimmomatic-0.39-2/ \
        --run-trim-repetitive \
        --fastqc fastqc \
        --sequencer-source none \
        --trimmomatic-options "LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:50" \
        --bypass-trf \
        --reorder


#### READ ERROR CORRECTION FOR PAIRED
tadpole.sh \
        in=../${SAMPLE_ID}/filtering_data/${SAMPLE_ID}_kneaddata_paired_1.fastq \
        in2=../${SAMPLE_ID}/filtering_data/${SAMPLE_ID}_kneaddata_paired_2.fastq \
        out=../${SAMPLE_ID}/filtering_data/${SAMPLE_ID}_ECC_paired_1.fastq \
        out2=../${SAMPLE_ID}/filtering_data/${SAMPLE_ID}_ECC_paired_2.fastq \
        mode=correct \
        ecc=t \
        prefilter=2 2>&1 \
        threads=${SLURM_CPUS_PER_TASK} \
        -Xmx32g | tee -a ../${SAMPLE_ID}/${SAMPLE_ID}_bbduk.log

#### READ ERROR CORRECTION FOR UNMATCHED
tadpole.sh \
        in=../${SAMPLE_ID}/filtering_data/${SAMPLE_ID}_kneaddata_unmatched_1.fastq \
        out=../${SAMPLE_ID}/filtering_data/${SAMPLE_ID}_ECC_unmatched_1.fastq \
        mode=correct \
        ecc=t \
        prefilter=2 2>&1 \
        threads=${SLURM_CPUS_PER_TASK} \
        -Xmx32g | tee -a ../${SAMPLE_ID}/${SAMPLE_ID}_bbduk.log

tadpole.sh \
        in=../${SAMPLE_ID}/filtering_data/${SAMPLE_ID}_kneaddata_unmatched_2.fastq \
        out=../${SAMPLE_ID}/filtering_data/${SAMPLE_ID}_ECC_unmatched_2.fastq \
        mode=correct \
        ecc=t \
        prefilter=2 2>&1 \
       threads=${SLURM_CPUS_PER_TASK} \
        -Xmx32g | tee -a ../${SAMPLE_ID}/${SAMPLE_ID}_bbduk.log

module load FastQC
fastqc ../${SAMPLE_ID}/filtering_data/${SAMPLE_ID}_NoTail_1.fastq.gz
fastqc ../${SAMPLE_ID}/filtering_data/${SAMPLE_ID}_NoTail_2.fastq.gz
mv ../${SAMPLE_ID}/filtering_data/${SAMPLE_ID}_*fastqc* ../FastQC_reports/01_NoTail/
fastqc ../${SAMPLE_ID}/filtering_data/${SAMPLE_ID}_AdaptTr_1.fastq.gz
fastqc ../${SAMPLE_ID}/filtering_data/${SAMPLE_ID}_AdaptTr_2.fastq.gz
mv ../${SAMPLE_ID}/filtering_data/${SAMPLE_ID}_*AdaptTr_*fastqc* ../FastQC_reports/02_AdaptTr_QualTr/
fastqc ../${SAMPLE_ID}/filtering_data/${SAMPLE_ID}_kneaddata_paired_1.fastq
fastqc ../${SAMPLE_ID}/filtering_data/${SAMPLE_ID}_kneaddata_paired_2.fastq
mv ../${SAMPLE_ID}/filtering_data/${SAMPLE_ID}_kneaddata*fastqc* ../FastQC_reports/03_Kneaddata/
fastqc ../${SAMPLE_ID}/filtering_data/${SAMPLE_ID}_ECC_paired_1.fastq
fastqc ../${SAMPLE_ID}/filtering_data/${SAMPLE_ID}_ECC_paired_2.fastq
fastqc ../${SAMPLE_ID}/filtering_data/${SAMPLE_ID}_ECC_unmatched_1.fastq
fastqc ../${SAMPLE_ID}/filtering_data/${SAMPLE_ID}_ECC_unmatched_2.fastq
mv ../${SAMPLE_ID}/filtering_data/${SAMPLE_ID}_ECC_*fastqc* ../FastQC_reports/04_RECC

##### PREPARING DIFFERENT FORMATS FOR DEDUPLICATION 

# FOR EXP03
cat ../${SAMPLE_ID}/filtering_data/${SAMPLE_ID}_ECC_paired_1.fastq \
        ../${SAMPLE_ID}/filtering_data/${SAMPLE_ID}_ECC_unmatched_1.fastq > \
        ../${SAMPLE_ID}/filtering_data/${SAMPLE_ID}_ECC_1.fastq

cat ../${SAMPLE_ID}/filtering_data/${SAMPLE_ID}_ECC_paired_2.fastq \
        ../${SAMPLE_ID}/filtering_data/${SAMPLE_ID}_ECC_unmatched_2.fastq > \
        ../${SAMPLE_ID}/filtering_data/${SAMPLE_ID}_ECC_2.fastq


# FOR EXP06
cat ../${SAMPLE_ID}/filtering_data/${SAMPLE_ID}_ECC_paired_1.fastq \
        ../${SAMPLE_ID}/filtering_data/${SAMPLE_ID}_ECC_paired_2.fastq \
        ../${SAMPLE_ID}/filtering_data/${SAMPLE_ID}_ECC_unmatched_1.fastq \
        ../${SAMPLE_ID}/filtering_data/${SAMPLE_ID}_ECC_unmatched_2.fastq > \
	../${SAMPLE_ID}/filtering_data/${SAMPLE_ID}_ECC_cat.fastq


module list

module purge

