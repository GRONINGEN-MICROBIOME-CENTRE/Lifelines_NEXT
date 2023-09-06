#!/bin/bash
#SBATCH --job-name=reads_QC
#SBATCH --error=reads_QC.err
#SBATCH --output=reads_QC.out
#SBATCH --mem=600gb
#SBATCH --time=16:59:00
#SBATCH --cpus-per-task=4
#SBATCH --open-mode=truncate
#SBATCH --partition=himem

SAMPLE_ID=$1
echo "SAMPLE_ID=${SAMPLE_ID}"

# --- LOAD MODULES --- 
module purge
module load BBMap

mkdir -p ../${SAMPLE_ID}/assembly/EXP11

module load SPAdes/3.15.3-GCC-11.2.0

if [ -f ../${SAMPLE_ID}/assembly/EXP11/spades.log ]
then 
	metaspades.py \
		--restart-from last \
		-o ../${SAMPLE_ID}/assembly/EXP11 \
        	-m $((${SLURM_MEM_PER_NODE} / 1000)) \
        	--threads ${SLURM_CPUS_PER_TASK}
else

### Assembly of paired-processed reads only
metaspades.py \
	-1 ../${SAMPLE_ID}/filtering_data/${SAMPLE_ID}_kneaddata_paired_1.fastq \
	-2 ../${SAMPLE_ID}/filtering_data/${SAMPLE_ID}_kneaddata_paired_2.fastq \
        -o ../${SAMPLE_ID}/assembly/EXP11 \
	-m $((${SLURM_MEM_PER_NODE} / 1000)) \
        --threads ${SLURM_CPUS_PER_TASK}

fi

module list

module purge

