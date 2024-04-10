#!/bin/bash
#SBATCH --job-name=PostDiscovery
#SBATCH --error=./err/08.sti/VD_Chiliadal_%A_%a.err
#SBATCH --output=./out/08.sti/VD_Chiliadal_%A_%a.out
#SBATCH --mem=16gb
#SBATCH --time=25:59:00
#SBATCH --cpus-per-task=2
#SBATCH --open-mode=truncate

##############################################
# Use-cases: in-house negative controls only #
##############################################

SAMPLE_LIST=$1

echo ${SAMPLE_LIST}

SAMPLE_ID=$(sed "${SLURM_ARRAY_TASK_ID}q;d" ${SAMPLE_LIST} | cut -d "_" -f1)

echo "SAMPLE_ID=${SAMPLE_ID}"

# --- SIMULATE TOF ---
grep '>' ../SAMPLES/${SAMPLE_ID}/01_sc_assembly/${SAMPLE_ID}_extended_contigs.fasta | \
	sed 's/>//g' > ../SAMPLES/${SAMPLE_ID}/01_sc_assembly/contigs_for_dereplication

# --- LOAD MODULES --- 
module purge
module load bioawk
module list

# --- GET LENGTH AFTER EXTENSION ---
bioawk -c fastx '{ print $name, length($seq) }' < \
	../SAMPLES/${SAMPLE_ID}/01_sc_assembly/${SAMPLE_ID}_extended_contigs.fasta > \
	../SAMPLES/${SAMPLE_ID}/01_sc_assembly/contigs_for_dereplication_length

# --- LOAD MODULES ---
module purge
module load R

# --- GENERATING EXTENDED_TOF ---
Rscript New_controls_contigs_ID_and_metadata.R /scratch/p282752/ANALYSIS_CHILIADAL/SAMPLES/${SAMPLE_ID}/01_sc_assembly/

# --- SPLITTING FASTAs ---
mkdir -p ../SAMPLES/${SAMPLE_ID}/01_sc_assembly/yasli

cd ../SAMPLES/${SAMPLE_ID}/01_sc_assembly/yasli
/scratch/p282752/ANALYSIS_CHILIADAL/scripts/fastasplitn \
	../${SAMPLE_ID}_extended_contigs.fasta \
	100

# --- SPLITTING Extended_TOF ---

for i in frag*.fa; do 
	i=${i%.fa}
	echo ${i} >> frag.list
	grep '>' ${i}.fa | sed 's/>//g' > ${i}_contigs
	Rscript /scratch/p282752/ANALYSIS_CHILIADAL/scripts/Split_ETOF.R ${i}
	rm ${i}_contigs
done

# --- RENAME CONTIGS ---
sbatch --array=1-100 /scratch/p282752/ANALYSIS_CHILIADAL/scripts/08.dIb_controls.sh frag.list ${SAMPLE_ID}

rm /scratch/p282752/ANALYSIS_CHILIADAL/SAMPLES/${SAMPLE_ID}/01_sc_assembly/contigs_for_dereplication
rm /scratch/p282752/ANALYSIS_CHILIADAL/SAMPLES/${SAMPLE_ID}/01_sc_assembly/contigs_for_dereplication_length

module purge
