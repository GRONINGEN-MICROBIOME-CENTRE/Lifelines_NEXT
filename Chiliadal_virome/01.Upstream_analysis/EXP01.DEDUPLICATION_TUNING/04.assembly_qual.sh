#!/bin/bash
#SBATCH --job-name=reads_QC
#SBATCH --error=reads_QC.err
#SBATCH --output=reads_QC.out
#SBATCH --mem=12gb
#SBATCH --time=4:59:00
#SBATCH --cpus-per-task=4
#SBATCH --open-mode=truncate

SAMPLE_ID=$1
echo "SAMPLE_ID=${SAMPLE_ID}"

EXP_ID=$2
echo "EXP_ID=${EXP_ID}"

# --- LOAD MODULES --- 
module purge
module load QUAST

### Assembly quality
quast.py \
	../${SAMPLE_ID}/assembly/${EXP_ID}/contigs.fasta \
        -o ../${SAMPLE_ID}/assembly/${EXP_ID}/quast \
	-m $((${SLURM_MEM_PER_NODE} / 1000)) \
        --threads ${SLURM_CPUS_PER_TASK}


module list

module purge

