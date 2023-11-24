#!/bin/bash
#SBATCH --job-name=QC_BBDuK_KneadData
#SBATCH --output=QC_BBDuK_KneadData_%A_%a.out
#SBATCH --mem=32gb
#SBATCH --time=02:59:00
#SBATCH --cpus-per-task=24
#SBATCH --export=NONE
#SBATCH --get-user-env=L
#SBATCH --partition=regular

sample_dir=$1 #directory with the FASTQ files 
sample_list=${sample_dir}/$2 #file with the list of all samples in the directory
SAMPLE_ID=$(sed "${SLURM_ARRAY_TASK_ID}q;d" ${sample_list})

echo -e '\n-------------------- WORKING WITH '${SAMPLE_ID}' SAMPLE --------------------'

echo -e '\n---- RENAMING FASTQ files ----'

mv ${sample_dir}/${SAMPLE_ID}/${SAMPLE_ID}*1.fq.gz ${sample_dir}/${SAMPLE_ID}/${SAMPLE_ID}_1.fq.gz
mv ${sample_dir}/${SAMPLE_ID}/${SAMPLE_ID}*2.fq.gz ${sample_dir}/${SAMPLE_ID}/${SAMPLE_ID}_2.fq.gz

echo -e '\n---- INITIAL Quality Check ----'

mkdir -p QC_RESULTS/INITIAL/${SAMPLE_ID}
mkdir -p QC_RESULTS/INITIAL/OUTPUT_files

# Clean environment, load modules 
module purge; module load FastQC/0.11.9-Java-11; module list

# Run FASTQC
fastqc -o QC_RESULTS/INITIAL/${SAMPLE_ID} -t ${SLURM_CPUS_PER_TASK} ${sample_dir}/${SAMPLE_ID}/${SAMPLE_ID}_1.fq.gz 
fastqc -o QC_RESULTS/INITIAL/${SAMPLE_ID} -t ${SLURM_CPUS_PER_TASK} ${sample_dir}/${SAMPLE_ID}/${SAMPLE_ID}_2.fq.gz

# Calculate number of raw reads and bases
N_reads_1=$( zcat ${sample_dir}/${SAMPLE_ID}/${SAMPLE_ID}_1.fq.gz | echo $((`wc -l`/4)) ) 
N_bases_1=$( zcat  ${sample_dir}/${SAMPLE_ID}/${SAMPLE_ID}_1.fq.gz  | paste - - - - | cut -f2 | wc -c )
N_reads_2=$( zcat ${sample_dir}/${SAMPLE_ID}/${SAMPLE_ID}_2.fq.gz | echo $((`wc -l`/4)) )
N_bases_2=$( zcat  ${sample_dir}/${SAMPLE_ID}/${SAMPLE_ID}_2.fq.gz | paste - - - - | cut -f2 | wc -c )

echo "Raw Reads FQ1: ${N_reads_1}" > QC_RESULTS/INITIAL/OUTPUT_files/${SAMPLE_ID}.out
echo "Raw Reads FQ2: ${N_reads_2}" >> QC_RESULTS/INITIAL/OUTPUT_files/${SAMPLE_ID}.out
echo "N bases FQ1: ${N_bases_1}" >> QC_RESULTS/INITIAL/OUTPUT_files/${SAMPLE_ID}.out
echo "N bases FQ2: ${N_bases_2}" >> QC_RESULTS/INITIAL/OUTPUT_files/${SAMPLE_ID}.out


echo -e '\n---- Copying files to tmpdir ----'

mkdir -p ${TMPDIR}/${SAMPLE_ID}/filtering_data/
cp ${sample_dir}/${SAMPLE_ID}/${SAMPLE_ID}_1.fq.gz ${TMPDIR}/${SAMPLE_ID}/filtering_data/
cp ${sample_dir}/${SAMPLE_ID}/${SAMPLE_ID}_2.fq.gz ${TMPDIR}/${SAMPLE_ID}/filtering_data/

echo -e '\n---- RUNNING BBDuK ----'

mkdir -p BBDuK_RESULTS/
mkdir -p BBDuK_RESULTS/LOG_files

# Clean environment, load modules 
module purge; module load BBMap; module list

bbduk.sh \
    in1=${TMPDIR}/${SAMPLE_ID}/filtering_data/${SAMPLE_ID}_1.fq.gz \
    in2=${TMPDIR}/${SAMPLE_ID}/filtering_data/${SAMPLE_ID}_2.fq.gz \
    out1=${TMPDIR}/${SAMPLE_ID}/filtering_data/${SAMPLE_ID}_AdaptTr_1.fq.gz \
    out2=${TMPDIR}/${SAMPLE_ID}/filtering_data/${SAMPLE_ID}_AdaptTr_2.fq.gz \
    ref=/scratch/hb-tifn/tools/GMH_pipeline/adapters/adapters.fa \
    ktrim=r k=23 mink=11 hdist=1 tpe tbo 2>&1 \
    threads=${SLURM_CPUS_PER_TASK} | tee -a BBDuK_RESULTS/LOG_files/${SAMPLE_ID}_bbduk.log


echo -e '\n---- RUNNING KneadData ----'

mkdir -p KneadData_RESULTS/${SAMPLE_ID}
mkdir -p KneadData_RESULTS/LOG_files/

# Clean environment, load new conda environment
module purge; ml Anaconda3; module list
conda activate /scratch/hb-tifn/condas/conda_biobakery3/; conda list

kneaddata \
    --input ${TMPDIR}/${SAMPLE_ID}/filtering_data/${SAMPLE_ID}_AdaptTr_1.fq.gz \
    --input ${TMPDIR}/${SAMPLE_ID}/filtering_data/${SAMPLE_ID}_AdaptTr_2.fq.gz  \
    --threads ${SLURM_CPUS_PER_TASK} \
    --processes 4 \
    --output-prefix ${SAMPLE_ID}_kneaddata \
    --output ${TMPDIR}/${SAMPLE_ID}/filtering_data \
    --log KneadData_RESULTS/LOG_files/${SAMPLE_ID}_kneaddata.log \
    -db /scratch/hb-tifn/DBs/human_genomes/GRCh38p13  \
    --trimmomatic /scratch/hb-tifn/condas/conda_biobakery4/share/trimmomatic-0.39-2/ \
    --run-trim-repetitive \
    --fastqc fastqc \
    --sequencer-source none \
    --trimmomatic-options "LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:50" \
    --bypass-trf \
    --reorder

rm ${TMPDIR}/${SAMPLE_ID}/filtering_data/*bowtie2* ${TMPDIR}/${SAMPLE_ID}/filtering_data/*kneaddata.trimmed* ${TMPDIR}/${SAMPLE_ID}/filtering_data/adapters.fa
rm -r ${TMPDIR}/${SAMPLE_ID}/filtering_data/fastqc 

echo -e '\n---- Moving results to SCRATCH ----'

LOG_FILE="KneadData_RESULTS/LOG_files/${SAMPLE_ID}_kneaddata.log"

if [ -f "$LOG_FILE" ] && ! grep -qi 'error' "$LOG_FILE"; then
    echo "KneadData is done"
    mv "${TMPDIR}/${SAMPLE_ID}/filtering_data/${SAMPLE_ID}_kneaddata_paired_1.fastq" "KneadData_RESULTS/${SAMPLE_ID}/"
    mv "${TMPDIR}/${SAMPLE_ID}/filtering_data/${SAMPLE_ID}_kneaddata_paired_2.fastq" "KneadData_RESULTS/${SAMPLE_ID}/"
    mv "${TMPDIR}/${SAMPLE_ID}/filtering_data/${SAMPLE_ID}_kneaddata_unmatched_1.fastq" "KneadData_RESULTS/${SAMPLE_ID}/"
    mv "${TMPDIR}/${SAMPLE_ID}/filtering_data/${SAMPLE_ID}_kneaddata_unmatched_2.fastq" "KneadData_RESULTS/${SAMPLE_ID}/"
else
    echo "Previous STEP failed or log file does not exist"
fi

echo "> Removing data from tmpdir"
rm -r ${TMPDIR}/${SAMPLE_ID}/filtering_data

echo -e '\n---- FINAL Quality Check ----'

mkdir -p QC_RESULTS/FINAL
mkdir -p QC_RESULTS/FINAL/CLEAN_READS/${SAMPLE_ID}
mkdir -p QC_RESULTS/FINAL/CLEAN_READS/OUTPUT_files
mkdir -p QC_RESULTS/FINAL/UNMATCHED_READS/${SAMPLE_ID}

module purge; module load FastQC/0.11.9-Java-11
fastqc -o QC_RESULTS/FINAL/CLEAN_READS/${SAMPLE_ID} -t 8 KneadData_RESULTS/${SAMPLE_ID}/${SAMPLE_ID}_kneaddata_paired_1.fastq
fastqc -o QC_RESULTS/FINAL/CLEAN_READS/${SAMPLE_ID} -t 8 KneadData_RESULTS/${SAMPLE_ID}/${SAMPLE_ID}_kneaddata_paired_2.fastq
fastqc -o QC_RESULTS/FINAL/UNMATCHED_READS/${SAMPLE_ID} -t 8 KneadData_RESULTS/${SAMPLE_ID}/${SAMPLE_ID}_kneaddata_unmatched_1.fastq
fastqc -o QC_RESULTS/FINAL/UNMATCHED_READS/${SAMPLE_ID} -t 8 KneadData_RESULTS/${SAMPLE_ID}/${SAMPLE_ID}_kneaddata_unmatched_2.fastq

# Calculate number of clean reads and bases
N_reads_1=$( cat KneadData_RESULTS/${SAMPLE_ID}/${SAMPLE_ID}_kneaddata_paired_1.fastq | echo $((`wc -l`/4)) ) 
N_bases_1=$( cat KneadData_RESULTS/${SAMPLE_ID}/${SAMPLE_ID}_kneaddata_paired_2.fastq  | paste - - - - | cut -f2 | wc -c )
N_reads_2=$( cat KneadData_RESULTS/${SAMPLE_ID}/${SAMPLE_ID}_kneaddata_paired_1.fastq | echo $((`wc -l`/4)) )
N_bases_2=$( cat KneadData_RESULTS/${SAMPLE_ID}/${SAMPLE_ID}_kneaddata_paired_2.fastq | paste - - - - | cut -f2 | wc -c )

echo "Clean Reads FQ1: ${N_reads_1}" > QC_RESULTS/FINAL/CLEAN_READS/OUTPUT_files/${SAMPLE_ID}.out
echo "Clean Reads FQ2: ${N_reads_2}" >> QC_RESULTS/FINAL/CLEAN_READS/OUTPUT_files/${SAMPLE_ID}.out
echo "N bases FQ1: ${N_bases_1}" >> QC_RESULTS/FINAL/CLEAN_READS/OUTPUT_files/${SAMPLE_ID}.out
echo "N bases FQ2: ${N_bases_2}" >> QC_RESULTS/FINAL/CLEAN_READS/OUTPUT_files/${SAMPLE_ID}.out

echo -e '\n---- Gzipping clean FASTQ files ----'

if [ $(echo $(cat KneadData_RESULTS/${SAMPLE_ID}/${SAMPLE_ID}_*_paired_1.fastq | wc -l)/4|bc)==$(echo $(cat KneadData_RESULTS/${SAMPLE_ID}/${SAMPLE_ID}_*_paired_2.fastq|wc -l)/4|bc) ]; then
	echo "Forward and reverse clean FASTQs have the same number of lines. They seem to be paired."
	cat KneadData_RESULTS/${SAMPLE_ID}/${SAMPLE_ID}_kneaddata_unmatched_1.fastq \
	KneadData_RESULTS/${SAMPLE_ID}/${SAMPLE_ID}_kneaddata_unmatched_2.fastq > \
	KneadData_RESULTS/${SAMPLE_ID}/${SAMPLE_ID}_kneaddata_unmatched.fastq
  pigz -p 2 KneadData_RESULTS/${SAMPLE_ID}/${SAMPLE_ID}*.fastq
fi

echo -e '\n---- Generating md5sums ----'
md5sum KneadData_RESULTS/${SAMPLE_ID}/${SAMPLE_ID}_kneaddata_paired_1.fastq.gz > KneadData_RESULTS/${SAMPLE_ID}/MD5.txt
md5sum KneadData_RESULTS/${SAMPLE_ID}/${SAMPLE_ID}_kneaddata_paired_2.fastq.gz >> KneadData_RESULTS/${SAMPLE_ID}/MD5.txt
md5sum KneadData_RESULTS/${SAMPLE_ID}/${SAMPLE_ID}_kneaddata_unmatched_1.fastq.gz >> KneadData_RESULTS/${SAMPLE_ID}/MD5.txt
md5sum KneadData_RESULTS/${SAMPLE_ID}/${SAMPLE_ID}_kneaddata_unmatched_2.fastq.gz >> KneadData_RESULTS/${SAMPLE_ID}/MD5.txt

echo -e '\n---- Generating folders with clean and unmatched reads. Input for assembly ----'

mkdir -p metaSPAdes/CLEAN_READS/${sample_dir}
mkdir -p metaSPAdes/UNMATCHED_READS/${sample_dir}

rsync -av $(find KneadData_RESULTS/${SAMPLE_ID} -name "*paired*.fastq.gz" -type f) metaSPAdes/CLEAN_READS/${sample_dir}
rsync -av $(find KneadData_RESULTS/${SAMPLE_ID} -name "*unmatched*fastq.gz" -type f) metaSPAdes/UNMATCHED_READS/${sample_dir}
ls metaSPAdes/CLEAN_READS/${sample_dir} | cut -f1 -d "_" | uniq > metaSPAdes/CLEAN_READS/${sample_dir}/list.txt
rm -r KneadData_RESULTS/${SAMPLE_ID}

echo -e '\n---- Checking errors ----'

mkdir -p PREPROCESSING_SUMMARY_RESULTS OUTPUT_files/

BBDUK_LOG="BBDuK_RESULTS/LOG_files/${SAMPLE_ID}_bbduk.log"
KNEADDATA_LOG="KneadData_RESULTS/LOG_files/${SAMPLE_ID}_kneaddata.log"

if [ -f "$BBDUK_LOG" ] && ! grep -qi 'error' "$BBDUK_LOG" && \
   [ -f "$KNEADDATA_LOG" ] && ! grep -qi 'error' "$KNEADDATA_LOG"; then
    echo "BBDuK processing for sample ${SAMPLE_ID} completed successfully." > "PREPROCESSING_SUMMARY_RESULTS/${SAMPLE_ID}_summary.txt"
    echo "KneadData processing for sample ${SAMPLE_ID} completed successfully." >> "PREPROCESSING_SUMMARY_RESULTS/${SAMPLE_ID}_summary.txt"
else
    echo "Previous step(s) failed or log file(s) do not exist."
fi

if [ $(find metaSPAdes/CLEAN_READS/${sample_dir} -type f -name "*${SAMPLE_ID}_kneaddata_paired*.fastq.gz" | wc -l) -eq 2 ]; then
    echo "The compressed clean FASTQ for sample ${SAMPLE_ID} files are ready for the assembly step." >> PREPROCESSING_SUMMARY_RESULTS/${SAMPLE_ID}_summary.txt
fi

if [ $(find metaSPAdes/UNMATCHED_READS/${sample_dir} -type f -name "*${SAMPLE_ID}_kneaddata_unmatched*.fastq.gz" | wc -l) -eq 1 ]; then
    echo -e "The compressed unmatched FASTQ files for sample ${SAMPLE_ID} are ready for the assembly step.\n" >> PREPROCESSING_SUMMARY_RESULTS/${SAMPLE_ID}_summary.txt
fi

echo -e '\n---- Read QC, adapter/contaminant trimming and compressing DONE ----'
