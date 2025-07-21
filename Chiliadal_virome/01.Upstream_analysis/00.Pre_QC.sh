#!/bin/bash
#SBATCH --job-name=Chiliadal_pre_rQC
#SBATCH --output=./out/00.pQC/PreQC_%A_%a.out
#SBATCH --mem=8gb
#SBATCH --time=01:30:00
#SBATCH --cpus-per-task=2

SAMPLE_LIST=$1

echo ${SAMPLE_LIST}

SAMPLE_ID=$(sed "${SLURM_ARRAY_TASK_ID}q;d" ${SAMPLE_LIST} | cut -d "_" -f1)

echo "SAMPLE_ID=${SAMPLE_ID}"

##### MODULES #####
module load FastQC/0.11.9-Java-11

## preQC
fastqc ../01.Initial_data_QC/FASTQC_preQC/${SAMPLE_ID}_*_R1_*.filt.fastq.gz
fastqc ../01.Initial_data_QC/FASTQC_preQC/${SAMPLE_ID}_*_R2_*.filt.fastq.gz

## calculate N of raw reads
## calculate N of bases

N_reads_1=$( zcat ../SAMPLES/${SAMPLE_ID}/raw_reads/${SAMPLE_ID}_1.fastq.gz | echo $((`wc -l`/4)) )
N_bases_1=$( zcat  ../SAMPLES/${SAMPLE_ID}/raw_reads/${SAMPLE_ID}_1.fastq.gz | paste - - - - | cut -f2 | wc -c )
N_reads_2=$( zcat ../SAMPLES/${SAMPLE_ID}/raw_reads/${SAMPLE_ID}_2.fastq.gz | echo $((`wc -l`/4)) )
N_bases_2=$( zcat  ../SAMPLES/${SAMPLE_ID}/raw_reads/${SAMPLE_ID}_2.fastq.gz | paste - - - - | cut -f2 | wc -c )

## removal of fastq.gz:
#rm ../01.Initial_data_QC/FASTQC_preQC/${SAMPLE_ID}_*_R*_*.filt.fastq.gz
echo "Raw Reads FQ1: ${N_reads_1}" >> ../upstream_stats/N_bases/${SAMPLE_ID}.out
echo "Raw Reads FQ2: ${N_reads_2}" >> ../upstream_stats/N_bases/${SAMPLE_ID}.out
echo "N bases FQ1: ${N_bases_1}" >> ../upstream_stats/N_bases/${SAMPLE_ID}.out
echo "N bases FQ2: ${N_bases_2}" >> ../upstream_stats/N_bases/${SAMPLE_ID}.out
