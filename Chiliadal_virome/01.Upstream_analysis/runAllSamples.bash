#!/bin/bash

for SAMPLE in $@; do
	sbatch --output ./out/${SAMPLE}.out --error ./err/${SAMPLE}.err --job-name ${SAMPLE} 01.Pre_QC.sh ${SAMPLE}
done
