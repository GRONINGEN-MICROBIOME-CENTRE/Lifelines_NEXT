#!/bin/bash
#SBATCH --job-name=PostDiscovery
#SBATCH --error=./err/08.sti/VD_Chiliadal_%A_%a.err
#SBATCH --output=./out/08.sti/VD_Chiliadal_%A_%a.out
#SBATCH --mem=16gb
#SBATCH --time=25:59:00
#SBATCH --cpus-per-task=2
#SBATCH --open-mode=truncate

#################################################
# Use-cases: BaseClear sequencing controls only #
#################################################

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

# --- RENAME CONTIGS ---
cp ../SAMPLES/${SAMPLE_ID}/01_sc_assembly/${SAMPLE_ID}_extended_contigs.fasta ../SAMPLES/${SAMPLE_ID}/01_sc_assembly/${SAMPLE_ID}_extended_renamed_contigs.fasta
awk 'NR>1' ../SAMPLES/${SAMPLE_ID}/01_sc_assembly/Extended_TOF  | awk '{print $15"\t"$17}' | while IFS=$'\t' read -r old_id new_id; do
        sed "s/>$old_id\b/>$new_id/" -i ../SAMPLES/${SAMPLE_ID}/01_sc_assembly/${SAMPLE_ID}_extended_renamed_contigs.fasta
done

rm ../SAMPLES/${SAMPLE_ID}/01_sc_assembly/contigs_for_dereplication
rm ../SAMPLES/${SAMPLE_ID}/01_sc_assembly/contigs_for_dereplication_length

module purge
