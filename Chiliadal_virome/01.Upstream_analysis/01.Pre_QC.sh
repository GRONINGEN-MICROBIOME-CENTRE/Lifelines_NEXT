#!/bin/bash
#SBATCH --job-name=QC
#SBATCH --output=QC.out
#SBATCH --error=QC.err
#SBATCH --mem=8gb
#SBATCH --time=01:30:00
#SBATCH --cpus-per-task=2

SAMPLE_ID=$1
echo "SAMPLE_ID=${SAMPLE_ID}"

##### MODULES #####
module load FastQC/0.11.9-Java-11

## preQC
fastqc ../01.Initial_data_QC/FASTQC_preQC/${SAMPLE_ID}_*_R1_*.filt.fastq.gz 
fastqc ../01.Initial_data_QC/FASTQC_preQC/${SAMPLE_ID}_*_R2_*.filt.fastq.gz

## calculate N of raw reads
## calculate N of bases

N_reads_1=$( zcat ../01.Initial_data_QC/FASTQC_preQC/${SAMPLE_ID}_*_R1_*.filt.fastq.gz | echo $((`wc -l`/4)) ) 
N_bases_1=$( zcat  ../01.Initial_data_QC/FASTQC_preQC/${SAMPLE_ID}_*_R1_*.filt.fastq.gz | paste - - - - | cut -f2 | wc -c )
N_reads_2=$( zcat ../01.Initial_data_QC/FASTQC_preQC/${SAMPLE_ID}_*_R2_*.filt.fastq.gz | echo $((`wc -l`/4)) )
N_bases_2=$( zcat  ../01.Initial_data_QC/FASTQC_preQC/${SAMPLE_ID}_*_R2_*.filt.fastq.gz | paste - - - - | cut -f2 | wc -c )

## removal of fastq.gz:
rm ../01.Initial_data_QC/FASTQC_preQC/${SAMPLE_ID}_*_R*_*.filt.fastq.gz
echo "Raw Reads FQ1: ${N_reads_1}" >> ./out/${SAMPLE_ID}.out
echo "Raw Reads FQ2: ${N_reads_2}" >> ./out/${SAMPLE_ID}.out
echo "N bases FQ1: ${N_bases_1}" >> ./out/${SAMPLE_ID}.out
echo "N bases FQ2: ${N_bases_2}" >> ./out/${SAMPLE_ID}.out
