#!/bin/bash

for SAMPLE in $@; do
 	sbatch --output ./out/${SAMPLE}_01.out --error ./err/${SAMPLE}_01.err --job-name rQC_${SAMPLE} 01.reads_QC.sh ${SAMPLE}
done
