#!/bin/bash
#SBATCH --job-name=Chiliadal_rQC
#SBATCH --error=./err/Chiliadal_%A_%a.err
#SBATCH --output=./out/Chiliadal_%A_%a.out
#SBATCH --mem=48gb
#SBATCH --time=12:59:00
#SBATCH --cpus-per-task=2
#SBATCH --open-mode=truncate

#SAMPLE_ID=$1
#echo "SAMPLE_ID=${SAMPLE_ID}"

SAMPLE_LIST=$1

echo ${SAMPLE_LIST}

SAMPLE_ID=$(sed "${SLURM_ARRAY_TASK_ID}q;d" ${SAMPLE_LIST} | cut -d "_" -f1)

echo "SAMPLE_ID=${SAMPLE_ID}"


### WORKING IN $TMPDIR
mkdir -p ${TMPDIR}/${SAMPLE_ID}/filtering_data/

echo "> copying files to tmpdir"
cp ../SAMPLES/${SAMPLE_ID}/raw_reads/${SAMPLE_ID}_1.fastq.gz ${TMPDIR}/${SAMPLE_ID}/filtering_data/
cp ../SAMPLES/${SAMPLE_ID}/raw_reads/${SAMPLE_ID}_2.fastq.gz ${TMPDIR}/${SAMPLE_ID}/filtering_data/

# --- LOAD MODULES --- 
module purge
module load BBMap

#### TRIMMING ADAPTERS
echo "> Trimming adapters" 

bbduk.sh \
        in1=${TMPDIR}/${SAMPLE_ID}/filtering_data/${SAMPLE_ID}_1.fastq.gz \
        in2=${TMPDIR}/${SAMPLE_ID}/filtering_data/${SAMPLE_ID}_2.fastq.gz \
        out1=${TMPDIR}/${SAMPLE_ID}/filtering_data/${SAMPLE_ID}_AdaptTr_1.fastq.gz \
        out2=${TMPDIR}/${SAMPLE_ID}/filtering_data/${SAMPLE_ID}_AdaptTr_2.fastq.gz \
        ref=/scratch/p282752/Data_for_HiC/adapters_UPD_IDT.fa \
        ktrim=r k=23 mink=11 hdist=1 tpe tbo 2>&1 \
        threads=${SLURM_CPUS_PER_TASK} \
        -Xmx$((${SLURM_MEM_PER_NODE} / 1024))g | tee -a ../SAMPLES/${SAMPLE_ID}/${SAMPLE_ID}_bbduk.log

echo "> Removing raw reads"

rm ${TMPDIR}/${SAMPLE_ID}/filtering_data/${SAMPLE_ID}_1.fastq.gz
rm ${TMPDIR}/${SAMPLE_ID}/filtering_data/${SAMPLE_ID}_2.fastq.gz

#### TRIMMING ADAPTASE-INTRODUCED TAILS 
echo "> Trimming adaptase-introduced tails"

bbduk.sh \
	in=${TMPDIR}/${SAMPLE_ID}/filtering_data/${SAMPLE_ID}_AdaptTr_1.fastq.gz \
	out=${TMPDIR}/${SAMPLE_ID}/filtering_data/${SAMPLE_ID}_NoTail_1.fastq.gz \
	ftr2=11 tpe tbo 2>&1 \
        threads=${SLURM_CPUS_PER_TASK} \
	-Xmx$((${SLURM_MEM_PER_NODE} / 1024))g | tee -a ../SAMPLES/${SAMPLE_ID}/${SAMPLE_ID}_bbduk.log
bbduk.sh \
        in=${TMPDIR}/${SAMPLE_ID}/filtering_data/${SAMPLE_ID}_AdaptTr_2.fastq.gz \
        out=${TMPDIR}/${SAMPLE_ID}/filtering_data/${SAMPLE_ID}_NoTail_2.fastq.gz \
        ftl=11 tpe tbo 2>&1 \
        threads=${SLURM_CPUS_PER_TASK} \
        -Xmx$((${SLURM_MEM_PER_NODE} / 1024))g | tee -a ../SAMPLES/${SAMPLE_ID}/${SAMPLE_ID}_bbduk.log

echo "> Removing adapter-trimmed fastqs"

rm ${TMPDIR}/${SAMPLE_ID}/filtering_data/${SAMPLE_ID}_AdaptTr_1.fastq.gz
rm ${TMPDIR}/${SAMPLE_ID}/filtering_data/${SAMPLE_ID}_AdaptTr_2.fastq.gz

#### FILTERING HUMAN READS & LOW QUALITY READS
echo "> Loading Anaconda3 and conda environment"

module load Anaconda3/2022.05
conda activate /scratch/hb-tifn/condas/conda_biobakery3/

kneaddata \
	--input ${TMPDIR}/${SAMPLE_ID}/filtering_data/${SAMPLE_ID}_NoTail_1.fastq.gz \
        --input ${TMPDIR}/${SAMPLE_ID}/filtering_data/${SAMPLE_ID}_NoTail_2.fastq.gz \
        --threads ${SLURM_CPUS_PER_TASK} \
        --processes 4 \
        --output-prefix ${SAMPLE_ID}_kneaddata \
        --output ${TMPDIR}/${SAMPLE_ID}/filtering_data \
        --log ../SAMPLES/${SAMPLE_ID}/${SAMPLE_ID}_kneaddata.log \
        -db /scratch/hb-tifn/DBs/human_genomes/hg37dec_v0.1  \
        --trimmomatic /scratch/hb-tifn/condas/conda_biobakery4/share/trimmomatic-0.39-2/ \
        --run-trim-repetitive \
        --fastqc fastqc \
        --sequencer-source none \
        --trimmomatic-options "LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:50" \
        --bypass-trf \
        --reorder

echo "> Removing kneaddata byproducts and tail-trimmed fastqs"

rm ${TMPDIR}/${SAMPLE_ID}/filtering_data/${SAMPLE_ID}_NoTail_1.fastq.gz
rm ${TMPDIR}/${SAMPLE_ID}/filtering_data/${SAMPLE_ID}_NoTail_2.fastq.gz
rm ${TMPDIR}/${SAMPLE_ID}/filtering_data/${SAMPLE_ID}_kneaddata_hg37dec_v0.1_bowtie2_paired_contam_1.fastq
rm ${TMPDIR}/${SAMPLE_ID}/filtering_data/${SAMPLE_ID}_kneaddata_hg37dec_v0.1_bowtie2_paired_contam_2.fastq
rm ${TMPDIR}/${SAMPLE_ID}/filtering_data/${SAMPLE_ID}_kneaddata_hg37dec_v0.1_bowtie2_unmatched_1_contam.fastq
rm ${TMPDIR}/${SAMPLE_ID}/filtering_data/${SAMPLE_ID}_kneaddata_hg37dec_v0.1_bowtie2_unmatched_2_contam.fastq
rm ${TMPDIR}/${SAMPLE_ID}/filtering_data/${SAMPLE_ID}_kneaddata.trimmed.1.fastq
rm ${TMPDIR}/${SAMPLE_ID}/filtering_data/${SAMPLE_ID}_kneaddata.trimmed.2.fastq
rm ${TMPDIR}/${SAMPLE_ID}/filtering_data/${SAMPLE_ID}_kneaddata.trimmed.single.1.fastq
rm ${TMPDIR}/${SAMPLE_ID}/filtering_data/${SAMPLE_ID}_kneaddata.trimmed.single.2.fastq

#### READ ERROR CORRECTION FOR PAIRED
echo "> Correcting of read errors"

tadpole.sh \
        in=${TMPDIR}/${SAMPLE_ID}/filtering_data/${SAMPLE_ID}_kneaddata_paired_1.fastq \
        in2=${TMPDIR}/${SAMPLE_ID}/filtering_data/${SAMPLE_ID}_kneaddata_paired_2.fastq \
        out=${TMPDIR}/${SAMPLE_ID}/filtering_data/${SAMPLE_ID}_ECC_paired_1.fastq \
        out2=${TMPDIR}/${SAMPLE_ID}/filtering_data/${SAMPLE_ID}_ECC_paired_2.fastq \
        mode=correct \
        ecc=t \
        prefilter=2 2>&1 \
        threads=${SLURM_CPUS_PER_TASK} \
        -Xmx$((${SLURM_MEM_PER_NODE} / 1024))g | tee -a ../SAMPLES/${SAMPLE_ID}/${SAMPLE_ID}_bbduk.log

#### READ ERROR CORRECTION FOR UNMATCHED
tadpole.sh \
        in=${TMPDIR}/${SAMPLE_ID}/filtering_data/${SAMPLE_ID}_kneaddata_unmatched_1.fastq \
        out=${TMPDIR}/${SAMPLE_ID}/filtering_data/${SAMPLE_ID}_ECC_unmatched_1.fastq \
        mode=correct \
        ecc=t \
        prefilter=2 2>&1 \
        threads=${SLURM_CPUS_PER_TASK} \
        -Xmx$((${SLURM_MEM_PER_NODE} / 1024))g | tee -a ../SAMPLES/${SAMPLE_ID}/${SAMPLE_ID}_bbduk.log

tadpole.sh \
        in=${TMPDIR}/${SAMPLE_ID}/filtering_data/${SAMPLE_ID}_kneaddata_unmatched_2.fastq \
        out=${TMPDIR}/${SAMPLE_ID}/filtering_data/${SAMPLE_ID}_ECC_unmatched_2.fastq \
        mode=correct \
        ecc=t \
        prefilter=2 2>&1 \
        threads=${SLURM_CPUS_PER_TASK} \
        -Xmx$((${SLURM_MEM_PER_NODE} / 1024))g | tee -a ../SAMPLES/${SAMPLE_ID}/${SAMPLE_ID}_bbduk.log

echo "> Removing kneaddata paired and unmatched reads"

rm ${TMPDIR}/${SAMPLE_ID}/filtering_data/${SAMPLE_ID}_kneaddata_paired_1.fastq
rm ${TMPDIR}/${SAMPLE_ID}/filtering_data/${SAMPLE_ID}_kneaddata_paired_2.fastq
rm ${TMPDIR}/${SAMPLE_ID}/filtering_data/${SAMPLE_ID}_kneaddata_unmatched_1.fastq
rm ${TMPDIR}/${SAMPLE_ID}/filtering_data/${SAMPLE_ID}_kneaddata_unmatched_2.fastq

#### READ DEDUPLICATION
echo "> Deduplicating reads"

clumpify.sh \
        in=${TMPDIR}/${SAMPLE_ID}/filtering_data/${SAMPLE_ID}_ECC_paired_1.fastq \
        in2=${TMPDIR}/${SAMPLE_ID}/filtering_data/${SAMPLE_ID}_ECC_paired_2.fastq \
        out=${TMPDIR}/${SAMPLE_ID}/filtering_data/${SAMPLE_ID}_dedup_paired_1.fastq \
        out2=${TMPDIR}/${SAMPLE_ID}/filtering_data/${SAMPLE_ID}_dedup_paired_2.fastq \
        dedupe=t \
        subs=0 \
        passes=2 \
        deletetemp=t 2>&1 \
       threads=${SLURM_CPUS_PER_TASK} \
       -Xmx$((${SLURM_MEM_PER_NODE} / 1024))g | tee -a ../SAMPLES/${SAMPLE_ID}/${SAMPLE_ID}_bbduk.log

clumpify.sh \
        in=${TMPDIR}/${SAMPLE_ID}/filtering_data/${SAMPLE_ID}_ECC_unmatched_1.fastq \
        out=${TMPDIR}/${SAMPLE_ID}/filtering_data/${SAMPLE_ID}_dedup_unmatched_1.fastq \
        dedupe=t \
        subs=0 \
        passes=2 \
        deletetemp=t 2>&1 \
       	threads=${SLURM_CPUS_PER_TASK} \
       -Xmx$((${SLURM_MEM_PER_NODE} / 1024))g | tee -a ../SAMPLES/${SAMPLE_ID}/${SAMPLE_ID}_bbduk.log

clumpify.sh \
        in=${TMPDIR}/${SAMPLE_ID}/filtering_data/${SAMPLE_ID}_ECC_unmatched_2.fastq \
        out=${TMPDIR}/${SAMPLE_ID}/filtering_data/${SAMPLE_ID}_dedup_unmatched_2.fastq \
        dedupe=t \
        subs=0 \
        passes=2 \
        deletetemp=t 2>&1 \
        threads=${SLURM_CPUS_PER_TASK} \
       -Xmx$((${SLURM_MEM_PER_NODE} / 1024))g | tee -a ../SAMPLES/${SAMPLE_ID}/${SAMPLE_ID}_bbduk.log

echo "> Moving resulting clean reads to scratch"

if [ $(grep 'Done!' ../SAMPLES/${SAMPLE_ID}/${SAMPLE_ID}_bbduk.log | wc -l)==3 ]; then
	echo "Clumpify is done"
	mkdir -p ../SAMPLES/${SAMPLE_ID}/clean_reads/
	mv ${TMPDIR}/${SAMPLE_ID}/filtering_data/${SAMPLE_ID}_dedup_paired_1.fastq ../SAMPLES/${SAMPLE_ID}/clean_reads/
	mv ${TMPDIR}/${SAMPLE_ID}/filtering_data/${SAMPLE_ID}_dedup_paired_2.fastq ../SAMPLES/${SAMPLE_ID}/clean_reads/
	mv ${TMPDIR}/${SAMPLE_ID}/filtering_data/${SAMPLE_ID}_dedup_unmatched_1.fastq ../SAMPLES/${SAMPLE_ID}/clean_reads/
	mv ${TMPDIR}/${SAMPLE_ID}/filtering_data/${SAMPLE_ID}_dedup_unmatched_2.fastq ../SAMPLES/${SAMPLE_ID}/clean_reads/
else
	echo "Deduplication or earlier step is corrupted"
fi

echo "> Removing data from tmpdir"
rm -r ${TMPDIR}/${SAMPLE_ID}/filtering_data

##### CHECKING QUALITY OF CLEAN READS
echo "> Running FastQC on cleanreads"

if [ $(grep 'Killed' ../SAMPLES/${SAMPLE_ID}/${SAMPLE_ID}_bbduk.log | wc -l)==0 ]; then
	module load FastQC
	fastqc -o ../FastQC_reports/ -t ${SLURM_CPUS_PER_TASK} ../SAMPLES/${SAMPLE_ID}/clean_reads/${SAMPLE_ID}_dedup_paired_1.fastq
	fastqc -o ../FastQC_reports/ -t ${SLURM_CPUS_PER_TASK} ../SAMPLES/${SAMPLE_ID}/clean_reads/${SAMPLE_ID}_dedup_paired_2.fastq
	fastqc -o ../FastQC_reports/ -t ${SLURM_CPUS_PER_TASK} ../SAMPLES/${SAMPLE_ID}/clean_reads/${SAMPLE_ID}_dedup_unmatched_1.fastq
	fastqc -o ../FastQC_reports/ -t ${SLURM_CPUS_PER_TASK} ../SAMPLES/${SAMPLE_ID}/clean_reads/${SAMPLE_ID}_dedup_unmatched_2.fastq
fi

if [ $(echo $(cat ../SAMPLES/${SAMPLE_ID}/clean_reads/${SAMPLE_ID}_dedup_paired_1.fastq|wc -l)/4|bc) == $(echo $(cat ../SAMPLES/${SAMPLE_ID}/clean_reads/${SAMPLE_ID}_dedup_paired_2.fastq|wc -l)/4|bc) ]; then
	echo "Reads seems to be paired"
	cat ../SAMPLES/${SAMPLE_ID}/clean_reads/${SAMPLE_ID}_dedup_unmatched_1.fastq \
	../SAMPLES/${SAMPLE_ID}/clean_reads/${SAMPLE_ID}_dedup_unmatched_2.fastq > \
	../SAMPLES/${SAMPLE_ID}/clean_reads/${SAMPLE_ID}_dedup_unmatched.fastq
	echo ">Compressing reads"
	pigz -p 2 ../SAMPLES/${SAMPLE_ID}/clean_reads/*.fastq
fi


echo "> Generating md5sums"
md5sum ../SAMPLES/${SAMPLE_ID}/clean_reads/${SAMPLE_ID}_dedup_paired_1.fastq.gz > ../SAMPLES/${SAMPLE_ID}/MD5.txt
md5sum ../SAMPLES/${SAMPLE_ID}/clean_reads/${SAMPLE_ID}_dedup_paired_2.fastq.gz >> ../SAMPLES/${SAMPLE_ID}/MD5.txt
md5sum ../SAMPLES/${SAMPLE_ID}/clean_reads/${SAMPLE_ID}_dedup_unmatched_1.fastq.gz >> ../SAMPLES/${SAMPLE_ID}/MD5.txt
md5sum ../SAMPLES/${SAMPLE_ID}/clean_reads/${SAMPLE_ID}_dedup_unmatched_2.fastq.gz >> ../SAMPLES/${SAMPLE_ID}/MD5.txt

echo "> Launching sc assembly"
bash runAllSamples_02.bash ${SAMPLE_ID}

module list

module purge

