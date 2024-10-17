#!/bin/bash
#SBATCH --job-name=inStrain_profile
#SBATCH --output=./out/13.dct/IS_compare.out
#SBATCH --time=168:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=256GB

ml Python/3.10.8-GCCcore-12.2.0
source /scratch/p282752/tools/python_envs/instrain/bin/activate
ml SAMtools

inStrain compare \
	-p 8 \
	-cov 0.5 \
	--min_cov 1 \
	-i ../VIR_DB/decontamination/ALL_PROFILE/* \
	-o ../VIR_DB/decontamination/instrain_compare_ALL \
	-sc ../VIR_DB/decontamination/vOTUs_shared_by_samples_and_ncs \
	-s ../VIR_DB/decontamination/genomes.stb \
	--database_mode \
        --skip_plot_generation	
deactivate
