#!/bin/bash
#SBATCH --job-name=waltersse_rQC
#SBATCH --error=./err/01.rQC/waltersse_%A_%a.err
#SBATCH --output=./out/01.rQC/waltersse_%A_%a.out
#SBATCH --mem=512gb
#SBATCH --time=48:59:00
#SBATCH --cpus-per-task=2
#SBATCH --open-mode=truncate

#SAMPLE_ID=$1
#echo "SAMPLE_ID=${SAMPLE_ID}"

SAMPLE_LIST=$1

echo ${SAMPLE_LIST}

SAMPLE_ID=$(sed "${SLURM_ARRAY_TASK_ID}q;d" ${SAMPLE_LIST} | cut -d "_" -f1)

echo "SAMPLE_ID=${SAMPLE_ID}"

# --- WORKING IN $TMPDIR ---
mkdir -p ${TMPDIR}/${SAMPLE_ID}/filtering_data/

echo "> copying files to tmpdir"
cp ../SAMPLES/${SAMPLE_ID}/raw_reads/${SAMPLE_ID}.fastq.gz ${TMPDIR}/${SAMPLE_ID}/filtering_data/

# --- LOADING MODULES --- 
module purge
module load BBMap

# --- TRIMMING ADAPTERS ---
echo "> Trimming adapters" 

bbduk.sh \
        in=${TMPDIR}/${SAMPLE_ID}/filtering_data/${SAMPLE_ID}.fastq.gz \
        out=${TMPDIR}/${SAMPLE_ID}/filtering_data/${SAMPLE_ID}_AdaptTr.fastq.gz \
        ref=/scratch/p309176/amg_paper/adapters_UPD_IDT.fa \
        ktrim=r k=23 mink=11 hdist=1 2>&1 \
	threads=${SLURM_CPUS_PER_TASK} \
        -Xmx$((${SLURM_MEM_PER_NODE} / 1024))g | tee -a ../SAMPLES/${SAMPLE_ID}/${SAMPLE_ID}_bbduk.log

# --- REMOVING RAW READS FASTQS---
echo "> Removing raw reads"

rm ${TMPDIR}/${SAMPLE_ID}/filtering_data/${SAMPLE_ID}.fastq.gz

# --- FILTERING HUMAN READS & LOW QUALITY READS ---
echo "> Loading Anaconda3 and conda environment"

module load Anaconda3/2022.05
conda activate /scratch/hb-tifn/condas/conda_biobakery3/

kneaddata \
	--input ${TMPDIR}/${SAMPLE_ID}/filtering_data/${SAMPLE_ID}_AdaptTr.fastq.gz \
        --threads ${SLURM_CPUS_PER_TASK} \
        --processes 4 \
        --output-prefix ${SAMPLE_ID}_kneaddata \
        --output ${TMPDIR}/${SAMPLE_ID}/filtering_data \
        --log ../SAMPLES/${SAMPLE_ID}/${SAMPLE_ID}_kneaddata.log \
        --reference-db /scratch/hb-tifn/DBs/human_genomes/GRCh38p13  \
        --trimmomatic /scratch/hb-tifn/condas/conda_biobakery4/share/trimmomatic-0.39-2/ \
        --run-trim-repetitive \
        --fastqc fastqc \
        --sequencer-source none \
        --trimmomatic-options "LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:50" \
        --bypass-trf

# --- REMOVING KNEADDATA BYPRODUCTS AND ADAPTER-TRIMMED FASTQS ---
echo "> Removing kneaddata byproducts and adapter-trimmed fastqs"

rm ${TMPDIR}/${SAMPLE_ID}/filtering_data/${SAMPLE_ID}_AdaptTr.fastq.gz
rm ${TMPDIR}/${SAMPLE_ID}/filtering_data/${SAMPLE_ID}_kneaddata_GRCh38p13_bowtie2_contam.fastq
rm ${TMPDIR}/${SAMPLE_ID}/filtering_data/${SAMPLE_ID}_kneaddata.trimmed.fastq

# --- CORRECTING READ ERRORS IN SINGLE-END READS ---

tadpole.sh \
        in=${TMPDIR}/${SAMPLE_ID}/filtering_data/${SAMPLE_ID}_kneaddata.fastq \
        out=${TMPDIR}/${SAMPLE_ID}/filtering_data/${SAMPLE_ID}_ECC.fastq \
        mode=correct \
        ecc=t \
        prefilter=2 \
        threads=${SLURM_CPUS_PER_TASK} \
        -Xmx$((${SLURM_MEM_PER_NODE} / 1024))g 2>&1 | tee -a ../SAMPLES/${SAMPLE_ID}/${SAMPLE_ID}_bbduk.log

# --- REMOVING KNEADDATA PAIRED AND UNMATCHED FASTQS ---
echo "> Removing kneaddata paired and unmatched reads"

rm ${TMPDIR}/${SAMPLE_ID}/filtering_data/${SAMPLE_ID}_kneaddata.fastq

# --- DEDUPLICATING READS ---
echo "> Deduplicating reads"

clumpify.sh \
        in=${TMPDIR}/${SAMPLE_ID}/filtering_data/${SAMPLE_ID}_ECC.fastq \
        out=${TMPDIR}/${SAMPLE_ID}/filtering_data/${SAMPLE_ID}_dedup.fastq \
        dedupe=t \
        subs=0 \
        deletetemp=t \
       	threads=${SLURM_CPUS_PER_TASK} \
       -Xmx$((${SLURM_MEM_PER_NODE} / 1024))g 2>&1 | tee -a ../SAMPLES/${SAMPLE_ID}/${SAMPLE_ID}_bbduk.log

# --- SWITCHING TO SCRATCH ---
echo "> Moving resulting clean reads to scratch"

if [ $(grep 'Done!' ../SAMPLES/${SAMPLE_ID}/${SAMPLE_ID}_bbduk.log | wc -l) == 1 ]; then
	echo "Clumpify is done"
	mkdir -p ../SAMPLES/${SAMPLE_ID}/clean_reads/
	mv ${TMPDIR}/${SAMPLE_ID}/filtering_data/${SAMPLE_ID}_dedup.fastq ../SAMPLES/${SAMPLE_ID}/clean_reads/
else
	echo "Deduplication or an earlier step is corrupted"
fi

# --- REMOVING DATA FROM TMPDIR ---
echo "> Removing data from tmpdir"
rm -r ${TMPDIR}/${SAMPLE_ID}/filtering_data


if [ -f ../SAMPLES/${SAMPLE_ID}/clean_reads/${SAMPLE_ID}_dedup.fastq ]; then
	
# --- COMPRESSING ALL FASTQS ---	
	echo ">Compressing reads"
	pigz -p ${SLURM_CPUS_PER_TASK} ../SAMPLES/${SAMPLE_ID}/clean_reads/*.fastq

# --- GENERATING MD5SUMS FOR STORAGE --- 	
	echo "> Generating md5sums"
	md5sum ../SAMPLES/${SAMPLE_ID}/clean_reads/${SAMPLE_ID}_dedup.fastq.gz > ../SAMPLES/${SAMPLE_ID}/MD5.txt

# --- CHECKING QUALITY OF CLEAN READS ---
	module load FastQC
        fastqc -o ../FastQC_reports/ -t ${SLURM_CPUS_PER_TASK} ../SAMPLES/${SAMPLE_ID}/clean_reads/${SAMPLE_ID}_dedup.fastq.gz

# --- LAUNCHING ASSEMBLY --- 	
	echo "> Launching sc assembly"
	bash runAllSamples_02_se.bash ${SAMPLE_ID}

fi

module list

module purge

