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

mkdir -p ../SAMPLES/${SAMPLE_ID}/01_sc_assembly

echo "> copying files to tmpdir"
mkdir -p ${TMPDIR}/${SAMPLE_ID}/raw_reads/
cp ../SAMPLES/${SAMPLE_ID}/raw_reads/${SAMPLE_ID}.fastq.gz ${TMPDIR}/${SAMPLE_ID}/raw_reads/
cp ../SAMPLES/${SAMPLE_ID}/raw_reads/${SAMPLE_ID}_1.fastq.gz ${TMPDIR}/${SAMPLE_ID}/raw_reads/
cp ../SAMPLES/${SAMPLE_ID}/raw_reads/${SAMPLE_ID}_2.fastq.gz ${TMPDIR}/${SAMPLE_ID}/raw_reads/

# --- LOAD MODULES --- 
module purge
module load SPAdes/3.15.3-GCC-11.2.0

echo "> running sc assembly"
# --- SPAdes assembly in sc mode ---
spades.py \
        --sc \
        --only-assembler \
        -k 21,33,55,77,99,127 \
        --pe1-1 ${TMPDIR}/${SAMPLE_ID}/raw_reads/${SAMPLE_ID}_1.fastq.gz \
        --pe1-2 ${TMPDIR}/${SAMPLE_ID}/raw_reads/${SAMPLE_ID}_2.fastq.gz \
        --pe1-s ${TMPDIR}/${SAMPLE_ID}/raw_reads/${SAMPLE_ID}.fastq.gz \
	-o ${TMPDIR}/${SAMPLE_ID}/01_sc_assembly \
        -m $((${SLURM_MEM_PER_NODE} / 1024)) \
        --threads ${SLURM_CPUS_PER_TASK}

echo "> Copying spades.log regardless of result"
cp ${TMPDIR}/${SAMPLE_ID}/01_sc_assembly/spades.log ../SAMPLES/${SAMPLE_ID}/01_sc_assembly/${SAMPLE_ID}_spades.log 

if [ -f ${TMPDIR}/${SAMPLE_ID}/01_sc_assembly/contigs.fasta ]; then
	echo "> Assembly has finished"
	echo "> Copying contigs and scaffolds"
	cp ${TMPDIR}/${SAMPLE_ID}/01_sc_assembly/contigs.fasta ../SAMPLES/${SAMPLE_ID}/01_sc_assembly/${SAMPLE_ID}_contigs.fasta
	cp ${TMPDIR}/${SAMPLE_ID}/01_sc_assembly/scaffolds.fasta ../SAMPLES/${SAMPLE_ID}/01_sc_assembly/${SAMPLE_ID}_scaffolds.fasta
fi

echo "> Removing data from tmpdir"
rm -r ${TMPDIR}/${SAMPLE_ID}/01_sc_assembly

# --- Assembly quality assessment ---
if [ -f ../SAMPLES/${SAMPLE_ID}/01_sc_assembly/${SAMPLE_ID}_contigs.fasta ]; then
	echo "> Assessing the assembly quality"
	module load QUAST
	quast.py \
        	../SAMPLES/${SAMPLE_ID}/01_sc_assembly/${SAMPLE_ID}_contigs.fasta \
        	-o ../SAMPLES/${SAMPLE_ID}/01_sc_assembly/quast \
        	-m $((${SLURM_MEM_PER_NODE} / 1024)) \
        	--threads ${SLURM_CPUS_PER_TASK}
fi

# --- Contig trimming ---
if [ -f ../SAMPLES/${SAMPLE_ID}/01_sc_assembly/${SAMPLE_ID}_contigs.fasta ]; then
        echo "> Trimming contigs to 1kbp"
	perl filter_contigs.pl 1000 ../SAMPLES/${SAMPLE_ID}/01_sc_assembly/${SAMPLE_ID}_contigs.fasta > ../SAMPLES/${SAMPLE_ID}/01_sc_assembly/${SAMPLE_ID}_contigs.min1kbp.fasta
	sed -i 's/>NODE/>'${SAMPLE_ID}'_NODE/g' ../SAMPLES/${SAMPLE_ID}/01_sc_assembly/${SAMPLE_ID}_contigs.min1kbp.fasta
fi

echo "> Generating md5sums"
md5sum ../SAMPLES/${SAMPLE_ID}/01_sc_assembly/${SAMPLE_ID}_contigs.fasta >> ../SAMPLES/${SAMPLE_ID}/MD5.txt 
md5sum ../SAMPLES/${SAMPLE_ID}/01_sc_assembly/${SAMPLE_ID}_scaffolds.fasta >> ../SAMPLES/${SAMPLE_ID}/MD5.txt

module list

module purge

