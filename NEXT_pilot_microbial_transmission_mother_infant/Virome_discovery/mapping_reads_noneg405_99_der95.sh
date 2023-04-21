#!/bin/bash
#SBATCH --job-name=mapping_reads
#SBATCH --error=./err/mapping_reads.err
#SBATCH --output=./out/mapping_reads.out
#SBATCH --mem=8gb
#SBATCH --time=11:29:00
#SBATCH --cpus-per-task=1
#SBATCH --open-mode=truncate

SAMPLE_ID=$1
echo "SAMPLE_ID=${SAMPLE_ID}"

BATCH_ID=$2
echo "BATCH_ID=${BATCH_ID}"

if [ -f ../BATCH${BATCH_ID}/${SAMPLE_ID}/clean_reads/${SAMPLE_ID}_kneaddata_cleaned_pair_2.fastq.gz ]; then
	gunzip ../BATCH${BATCH_ID}/${SAMPLE_ID}/clean_reads/*.gz
fi

# PURGING ENVIRUMENT 
module purge 

# --- LOAD MODULES --- 
module load Python
module load Bowtie2
module load SAMtools

bowtie2 \
        -x ../clean_viral/dbs/viral_decontaminated/viral_decontaminated_db \
        -1 ../BATCH${BATCH_ID}/${SAMPLE_ID}/clean_reads/${SAMPLE_ID}_kneaddata_cleaned_pair_1.fastq \
        -2 ../BATCH${BATCH_ID}/${SAMPLE_ID}/clean_reads/${SAMPLE_ID}_kneaddata_cleaned_pair_2.fastq \
        -S ../BATCH${BATCH_ID}/${SAMPLE_ID}/${SAMPLE_ID}_all_vir_alignments.sam \
        -p ${SLURM_CPUS_PER_TASK} \
        --no-unal \
        --end-to-end \
        &> ../clean_viral/VLP_to_viral_decontaminated/alignment_log/${SAMPLE_ID}.bowtie2.log

samtools view \
	-S ../BATCH${BATCH_ID}/${SAMPLE_ID}/${SAMPLE_ID}_all_vir_alignments.sam \
	-b \
	-o ../BATCH${BATCH_ID}/${SAMPLE_ID}/${SAMPLE_ID}_all_vir_alignments.bam

samtools sort \
	../BATCH${BATCH_ID}/${SAMPLE_ID}/${SAMPLE_ID}_all_vir_alignments.bam \
	-o ../BATCH${BATCH_ID}/${SAMPLE_ID}/${SAMPLE_ID}_all_vir_alignments.sorted.bam

samtools index \
	../BATCH${BATCH_ID}/${SAMPLE_ID}/${SAMPLE_ID}_all_vir_alignments.sorted.bam

samtools mpileup \
	../BATCH${BATCH_ID}/${SAMPLE_ID}/${SAMPLE_ID}_all_vir_alignments.sorted.bam \
	-A \
	-f ../clean_viral/dbs/viral_decontaminated/viral_noneg405_99_der95_decontaminated.fasta \
	-o ../clean_viral/VLP_to_viral_decontaminated/breadth_cov/${SAMPLE_ID}.A.cov.txt

samtools idxstats \
	../BATCH${BATCH_ID}/${SAMPLE_ID}/${SAMPLE_ID}_all_vir_alignments.sorted.bam \
	> ../clean_viral/VLP_to_viral_decontaminated/counts_table/${SAMPLE_ID}.counts.idx.txt

samtools flagstat \
	../BATCH${BATCH_ID}/${SAMPLE_ID}/${SAMPLE_ID}_all_vir_alignments.sorted.bam \
	 > ../clean_viral/VLP_to_viral_decontaminated/counts_table/${SAMPLE_ID}.flagstat.txt

mv ../BATCH${BATCH_ID}/${SAMPLE_ID}/${SAMPLE_ID}_all_vir_alignments.sorted.bam ../clean_viral/for_Alex/VLP_to_viral_decontaminated/

rm ../BATCH${BATCH_ID}/${SAMPLE_ID}/${SAMPLE_ID}_all_vir_alignments*

gzip ../BATCH${BATCH_ID}/${SAMPLE_ID}/clean_reads/*.fastq

# PURGING ENVIRUMENT
module purge
