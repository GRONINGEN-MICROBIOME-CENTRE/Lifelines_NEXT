#!/bin/bash
#SBATCH --job-name=PD_deRep
#SBATCH --output=./out/09.mdrp/PD_dRep_%A_%a.out
#SBATCH --mem=16gb
#SBATCH --time=01:30:00
#SBATCH --cpus-per-task=4
#SBATCH --open-mode=truncate

FRAG_LIST=$1

echo "frag list: ${FRAG_LIST}"

FRAG_ID=$(sed "${SLURM_ARRAY_TASK_ID}q;d" ${FRAG_LIST})

echo "FRAG=${FRAG_ID}"

## --- LOAD MODULES ---
module purge
module load BLAST+/2.13.0-gompi-2022a
module list

# --- DEREPLICATION ACCORDING TO MIUViG GUIDELINES ---

# Next, use megablast from blast+ package to perform all-vs-all blastn of sequences:
blastn \
    -query ../VIR_DB/VLP_MGS_DREP/${FRAG_ID}.fa \
    -db ../VIR_DB/VLP_MGS_DREP/NEXT_VIR_DB \
    -outfmt '6 std qlen slen' \
    -max_target_seqs 10000 \
    -out ../VIR_DB/VLP_MGS_DREP/OUTPUT/NEXT_${FRAG_ID}_blast.tsv \
    -num_threads ${SLURM_CPUS_PER_TASK}

echo "${FRAG_ID} blastn done!"

module purge

# removing the fragment:

if [ ! -f ../VIR_DB/VLP_MGS_DREP/${FRAG_ID}.fa ] || [ ! -f ../VIR_DB/VLP_MGS_DREP/OUTPUT/NEXT_${FRAG_ID}_blast.tsv ]; then
    echo "Input file(s) not found!"
    exit 1
fi

N_contigs=$(grep '>' ../VIR_DB/VLP_MGS_DREP/${FRAG_ID}.fa | wc -l)
N_uniq_queries=$(awk -F '\t' '{print $1}' ../VIR_DB/VLP_MGS_DREP/OUTPUT/NEXT_${FRAG_ID}_blast.tsv | uniq | wc -l)

if [ ${N_contigs} -eq ${N_uniq_queries} ]; then
    echo "Number of contigs and unique queries is equal"
    rm ../VIR_DB/VLP_MGS_DREP/${FRAG_ID}.fa
else
    echo "Number of contigs and unique queries is unequal"
fi

