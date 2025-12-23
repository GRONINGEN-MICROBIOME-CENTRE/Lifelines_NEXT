#!/bin/bash
#SBATCH --job-name=DGR_analysis 
#SBATCH --output=DGR_analysis.out
#SBATCH --mem=32gb
#SBATCH --time=06:59:00
#SBATCH --cpus-per-task=32
#SBATCH --open-mode=truncate
#SBATCH --partition=regular

contig_file=$1 #path to FASTA file with predicted viral contigs
contig_file_name="$(basename "${contig_file}")" #extract filename
contig_file_name="${contig_file_name%.*}" #extract filename without the extension

echo -e '\n---- RUNNING DGR identification ----'

mkdir -p DGR_IDENTIFICATION

# Clean environment, load conda environment 
ml purge; ml Anaconda3; ml list
conda activate /scratch/hb-llnext/conda_envs/DGR_analysis_env; conda list

# Run prodigal-gv to get the protein predictions in GFF format
python parallel-prodigal-gv.py \
	-f gff \
	-i $contig_file \
	-a DGR_IDENTIFICATION/${contig_file_name}_proteins.faa \
	-o DGR_IDENTIFICATION/${contig_file_name}_proteins.gff \
	-t ${SLURM_CPUS_PER_TASK}

rsync -av $contig_file DGR_IDENTIFICATION/${contig_file_name}_contigs.fna #avoids next steps to fail

# Search for reverse transcriptases  
./Search_for_RT.pl \
    -i ${contig_file_name} \
    -d DGR_IDENTIFICATION/ \
    -t ${SLURM_CPUS_PER_TASK}

# Search for VR and TR
./Look_for_repeats_and_targets.pl \
    -i ${contig_file_name} \
    -d DGR_IDENTIFICATION/ \
    -t ${SLURM_CPUS_PER_TASK}

# Evaluate target gene changes  
./Evaluate_target_gene_changes.pl \
    -i ${contig_file_name} \
    -d DGR_IDENTIFICATION/
    
# Generate output files after filtering those DGRs with <75% A mismatches 
awk -F'\t' '$16 > 0.75' DGR_IDENTIFICATION/${contig_file_name}/All_RT_1k_DGR_detection_filtered.tsv > DGR_IDENTIFICATION/${contig_file_name}/All_RT_1k_DGR_detection_filtered_075A.tsv
awk 'NR > 1 {print $1 ":" $3}' DGR_IDENTIFICATION/${contig_file_name}/All_RT_1k_DGR_detection_filtered_075A.tsv | grep -wf - DGR_IDENTIFICATION/${contig_file_name}/All_RT_1k_DGR_target_detection_filtered.tsv > DGR_IDENTIFICATION/${contig_file_name}/All_RT_1k_DGR_target_detection_filtered_075A_temp.tsv
{ head -n 1 DGR_IDENTIFICATION/${contig_file_name}/All_RT_1k_DGR_target_detection_filtered.tsv; cat DGR_IDENTIFICATION/${contig_file_name}/All_RT_1k_DGR_target_detection_filtered_075A_temp.tsv; } > DGR_IDENTIFICATION/${contig_file_name}/All_RT_1k_DGR_target_detection_filtered_075A.tsv
rm DGR_IDENTIFICATION/${contig_file_name}/All_RT_1k_DGR_target_detection_filtered_075A_temp.tsv

echo -e '\n---- DGR identification done ----'

