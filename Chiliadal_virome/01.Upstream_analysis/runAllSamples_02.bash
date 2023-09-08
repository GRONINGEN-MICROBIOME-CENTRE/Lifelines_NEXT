#!/bin/bash

for SAMPLE in $@; do
 	sbatch --output ./out/${SAMPLE}_02.out --error ./err/${SAMPLE}_02.err --job-name RAs_${SAMPLE} 02.sc_assembly.sh ${SAMPLE}
done
