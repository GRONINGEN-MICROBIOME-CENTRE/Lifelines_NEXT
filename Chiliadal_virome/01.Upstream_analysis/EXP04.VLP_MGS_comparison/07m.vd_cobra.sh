#!/bin/bash
#SBATCH --job-name=ViromeDiscovery
#SBATCH --output=./out/07.mcob/VD_PentaChiliadal_%A_%a.out
#SBATCH --mem=8gb
#SBATCH --time=10:59:00
#SBATCH --cpus-per-task=4
#SBATCH --open-mode=truncate

SAMPLE_LIST=$1

SAMPLE_DIR=$2

echo ${SAMPLE_LIST}
echo ${SAMPLE_DIR}

SAMPLE_ID=$(sed "${SLURM_ARRAY_TASK_ID}q;d" ${SAMPLE_LIST} | cut -d "_" -f1)

echo "SAMPLE_ID=${SAMPLE_ID}"

# --- PREPARING TMPDIR ---
mkdir -p ${TMPDIR}/${SAMPLE_ID}/COBRA

# --- LOAD MODULES --- 
module purge
module load Bowtie2
module load SAMtools
module list

cp ${SAMPLE_DIR}/${SAMPLE_ID}/virome_discovery/tidy/${SAMPLE_ID}_all_predicted_viral.fasta ${TMPDIR}/${SAMPLE_ID}/COBRA
sed -i 's/^[^_]*_/>/' ${TMPDIR}/${SAMPLE_ID}/COBRA/${SAMPLE_ID}_all_predicted_viral.fasta

bowtie2-build \
	/scratch/hb-llnext/LLNEXT_CONTIGS/${SAMPLE_ID}_metaspades_contigs.fa \
	${TMPDIR}/${SAMPLE_ID}/COBRA/${SAMPLE_ID}

bowtie2 \
	--very-sensitive \
	-x ${TMPDIR}/${SAMPLE_ID}/COBRA/${SAMPLE_ID} \
	-1 /scratch/hb-llnext/NEXT_CLEAN_READS_2025/paired_reads/${SAMPLE_ID}_kneaddata_paired_1.fastq.gz \
	-2 /scratch/hb-llnext/NEXT_CLEAN_READS_2025/paired_reads/${SAMPLE_ID}_kneaddata_paired_2.fastq.gz \
	--no-unal \
	--threads ${SLURM_CPUS_PER_TASK} \
	-S ${TMPDIR}/${SAMPLE_ID}/COBRA/${SAMPLE_ID}.sam

samtools view \
	-bS \
	${TMPDIR}/${SAMPLE_ID}/COBRA/${SAMPLE_ID}.sam \
	> ${TMPDIR}/${SAMPLE_ID}/COBRA/${SAMPLE_ID}.bam

samtools sort \
	-o ${TMPDIR}/${SAMPLE_ID}/COBRA/${SAMPLE_ID}_sorted.bam \
	${TMPDIR}/${SAMPLE_ID}/COBRA/${SAMPLE_ID}.bam

samtools index \
	-@ $((${SLURM_CPUS_PER_TASK}-1)) \
	${TMPDIR}/${SAMPLE_ID}/COBRA/${SAMPLE_ID}_sorted.bam

# --- LOAD MODULES --- 
module purge
module load Anaconda3
module list

source activate /scratch/hb-llnext/conda_envs/CoverM_env
conda list

coverm contig \
	--coupled /scratch/hb-llnext/NEXT_CLEAN_READS_2025/paired_reads/${SAMPLE_ID}_kneaddata_paired_1.fastq.gz /scratch/hb-llnext/NEXT_CLEAN_READS_2025/paired_reads/${SAMPLE_ID}_kneaddata_paired_2.fastq.gz \
       	--reference /scratch/hb-llnext/LLNEXT_CONTIGS/${SAMPLE_ID}_metaspades_contigs.fa \
	-o ${TMPDIR}/${SAMPLE_ID}/COBRA/${SAMPLE_ID}_coverage.txt

sed -i '1d' ${TMPDIR}/${SAMPLE_ID}/COBRA/${SAMPLE_ID}_coverage.txt 

conda deactivate

# --- LOAD MODULES --- 
module purge
module load Anaconda3
module load BLAST+
module list

source activate /scratch/hb-llnext/conda_envs/COBRA_env
module load Python/3.10.4-GCCcore-11.3.0
module load Pysam/0.19.1-GCC-11.3.0
conda list

#python /scratch/hb-llnext/conda_envs/COBRA_env/cobra/cobra.py \
python cobra_test.py \
    -f /scratch/hb-llnext/LLNEXT_CONTIGS/${SAMPLE_ID}_metaspades_contigs.fa \
    -q ${TMPDIR}/${SAMPLE_ID}/COBRA/${SAMPLE_ID}_all_predicted_viral.fasta \
    -c ${TMPDIR}/${SAMPLE_ID}/COBRA/${SAMPLE_ID}_coverage.txt \
    -m ${TMPDIR}/${SAMPLE_ID}/COBRA/${SAMPLE_ID}_sorted.bam \
    -a metaspades \
    -mink 21 \
    -maxk 55 \
    -o ${TMPDIR}/${SAMPLE_ID}/COBRA/RESULT

# --- MOVING TO SCRATCH --- 
sed "s/^>/>${SAMPLE_ID}_/" ${TMPDIR}/${SAMPLE_ID}/COBRA/RESULT/*.fasta > ${SAMPLE_DIR}/${SAMPLE_ID}/virome_discovery/tidy/${SAMPLE_ID}_extended_viral.fasta

N_EXT_CONTIGS=$(grep '>' ${SAMPLE_DIR}/${SAMPLE_ID}/virome_discovery/tidy/${SAMPLE_ID}_extended_viral.fasta | wc -l)
N_VIR_CONTIGS=$(grep '>' ${SAMPLE_DIR}/${SAMPLE_ID}/virome_discovery/tidy/${SAMPLE_ID}_all_predicted_viral.fasta | wc -l)

if [ $N_EXT_CONTIGS -eq $N_VIR_CONTIGS ]; then
    echo "Number of extended and original viral contigs is equal"
else
    echo "Number of extended and original viral contigs is unequal"
fi
echo "$N_EXT_CONTIGS contigs in the resulting fasta"

mkdir -p ${SAMPLE_DIR}/${SAMPLE_ID}/virome_discovery/tidy/COBRA
cp ${TMPDIR}/${SAMPLE_ID}/COBRA/RESULT/COBRA_joining_summary.txt ${SAMPLE_DIR}/${SAMPLE_ID}/virome_discovery/tidy/COBRA/
cp ${TMPDIR}/${SAMPLE_ID}/COBRA/RESULT/log ${SAMPLE_DIR}/${SAMPLE_ID}/virome_discovery/tidy/COBRA/COBRA.log

# --- ESTIMATING EXTENSION --- 

mkdir -p ${TMPDIR}/${SAMPLE_ID}/POST_COBRA_Q

# --- LOAD MODULES ---
module purge
module load QUAST
module list

# --- Original viral fasta quality assessment ---
quast.py \
        ${SAMPLE_DIR}/${SAMPLE_ID}/virome_discovery/tidy/${SAMPLE_ID}_all_predicted_viral.fasta \
        -o ${TMPDIR}/${SAMPLE_ID}/POST_COBRA_Q/pre \
        -m $((${SLURM_MEM_PER_NODE} / 1024)) \
        --threads ${SLURM_CPUS_PER_TASK}

# --- Extension quality assessment ---
quast.py \
        ${SAMPLE_DIR}/${SAMPLE_ID}/virome_discovery/tidy/${SAMPLE_ID}_extended_viral.fasta \
        -o ${TMPDIR}/${SAMPLE_ID}/POST_COBRA_Q/post \
        -m $((${SLURM_MEM_PER_NODE} / 1024)) \
        --threads ${SLURM_CPUS_PER_TASK}

TOTAL_LENGTH_PRE=$(grep 'Total length (>= 0 bp)' ${TMPDIR}/${SAMPLE_ID}/POST_COBRA_Q/pre/report.tsv | awk -F '\t' '{print $2}')

echo "TOTAL LENGTH OF VIRUSES: ${TOTAL_LENGTH_PRE}"

TOTAL_LENGTH_POST=$(grep 'Total length (>= 0 bp)' ${TMPDIR}/${SAMPLE_ID}/POST_COBRA_Q/post/report.tsv | awk -F '\t' '{print $2}')

echo "TOTAL LENGTH AFTER EXTENSION: ${TOTAL_LENGTH_POST}"

if [ ${TOTAL_LENGTH_PRE} -le ${TOTAL_LENGTH_POST} ]; then
    echo "Length of extended contigs is larger than that of original contigs"
else
    echo "Extension error, extended contigs are shorter"
fi

if [ $(grep 'for joining details of' ${SAMPLE_DIR}/${SAMPLE_ID}/virome_discovery/tidy/COBRA/COBRA.log | wc -l) -eq 1 ]; then
	echo "COBRA finished"
else
	echo "COBRA failed"
fi

module purge
