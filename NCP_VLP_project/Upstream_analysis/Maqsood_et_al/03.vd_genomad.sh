#!/bin/bash
#SBATCH --job-name=ViromeDiscovery
#SBATCH --error=./err/03.gnd/VD_AMG_%A_%a.err
#SBATCH --output=./out/03.gnd/VD_AMG_%A_%a.out
#SBATCH --mem=32gb
#SBATCH --time=05:59:00
#SBATCH --cpus-per-task=8
#SBATCH --open-mode=truncate

SAMPLE_LIST=$1

echo ${SAMPLE_LIST}

SAMPLE_ID=$(sed "${SLURM_ARRAY_TASK_ID}q;d" ${SAMPLE_LIST})

echo "SAMPLE_ID=${SAMPLE_ID}"

# --- WORKING IN $TMPDIR ---
mkdir -p ${TMPDIR}/${SAMPLE_ID}/genomad_run/

echo "> copying files to tmpdir"
cp ../SAMPLES/${SAMPLE_ID}/01_sc_assembly/${SAMPLE_ID}_contigs.min1kbp.fasta ${TMPDIR}/${SAMPLE_ID}/genomad_run/

export PATH=/scratch/hb-llnext/tools/mmseqs/bin/:$PATH

# --- LOAD MODULES --- 
module load ARAGORN/1.2.41-foss-2021b
module load Python/3.9.5-GCCcore-10.3.0
source /scratch/hb-llnext/python_venvs/geNomad/bin/activate
# --- RUNNING geNomad ---
echo "> Running geNomad"

genomad \
	end-to-end \
	--enable-score-calibration \
	--cleanup \
	${TMPDIR}/${SAMPLE_ID}/genomad_run/${SAMPLE_ID}_contigs.min1kbp.fasta \
	${TMPDIR}/${SAMPLE_ID}/genomad_run/geNomad \
	/scratch/hb-llnext/databases/genomad_db

genomad --version

deactivate

module list

module purge

# --- SAVING THE FILES ON SCRATCH ---
echo "> Moving the files to scratch"

mkdir -p ../SAMPLES/${SAMPLE_ID}/virome_discovery/geNomad/

cat ${TMPDIR}/${SAMPLE_ID}/genomad_run/geNomad/${SAMPLE_ID}_contigs.min1kbp_annotate.log \
	${TMPDIR}/${SAMPLE_ID}/genomad_run/geNomad/${SAMPLE_ID}_contigs.min1kbp_aggregated_classification.log \
	${TMPDIR}/${SAMPLE_ID}/genomad_run/geNomad/${SAMPLE_ID}_contigs.min1kbp_find_proviruses.log \
	${TMPDIR}/${SAMPLE_ID}/genomad_run/geNomad/${SAMPLE_ID}_contigs.min1kbp_marker_classification.log \
	${TMPDIR}/${SAMPLE_ID}/genomad_run/geNomad/${SAMPLE_ID}_contigs.min1kbp_nn_classification.log \
	${TMPDIR}/${SAMPLE_ID}/genomad_run/geNomad/${SAMPLE_ID}_contigs.min1kbp_score_calibration.log \
	${TMPDIR}/${SAMPLE_ID}/genomad_run/geNomad/${SAMPLE_ID}_contigs.min1kbp_summary.log > \
	${TMPDIR}/${SAMPLE_ID}/genomad_run/geNomad/${SAMPLE_ID}_contigs.min1kbp.log 

cp ${TMPDIR}/${SAMPLE_ID}/genomad_run/geNomad/${SAMPLE_ID}_contigs.min1kbp.log ../SAMPLES/${SAMPLE_ID}/virome_discovery/geNomad/
cp ${TMPDIR}/${SAMPLE_ID}/genomad_run/geNomad/${SAMPLE_ID}_contigs.min1kbp_summary/${SAMPLE_ID}_contigs.min1kbp_plasmid_summary.tsv ../SAMPLES/${SAMPLE_ID}/virome_discovery/geNomad/
cp ${TMPDIR}/${SAMPLE_ID}/genomad_run/geNomad/${SAMPLE_ID}_contigs.min1kbp_summary/${SAMPLE_ID}_contigs.min1kbp_virus_summary.tsv ../SAMPLES/${SAMPLE_ID}/virome_discovery/geNomad/

# --- REMOVING THE FOLDER FROM TMP ---
rm -r ${TMPDIR}/${SAMPLE_ID}/genomad_run

echo "> All done!"




