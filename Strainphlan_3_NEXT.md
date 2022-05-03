# Strainphlan 3.0 analysis for Lifelines NEXT pilot

Adapted from Biobakery (StrainPhlAn 3.0) 
Authors: Ranko Gascesa, Trishla Sinha, Asier Fernandez
Description: The script shows how strain profiling was performed using Strainphlan 3 for all samples, for all species 

## Step 1: Reconstruct all species stains

(see main microbiome profiling pipeline) 

## Step 2: Profile the clades present in the samples (profileClades.sh)

```
#!/bin/bash

#SBATCH --mem=16gb
#SBATCH --time=0-00:30:00
#SBATCH --cpus-per-task=4
#SBATCH --open-mode=truncate

# NOTES:
# > $1 is clade name
# > variables should be set to:
#   MARKERS_INPUTDATA: where markers for DATASET TO PROCESS are (must be bundled in one folder)
#   OUTFOLDER: where results go
#   MPA_MARKERS_EXTRACTED: folder with extracted metaphlan markers (extracted from metaphlan pkled markers database = MPA_PKL)
#   MPA_PKL: metaphlan pkled database
#   THREADS: NR of threads for main app
#   THREADS_RAXML: NR of threads for RAXML

# purge modules
module purge

# load conda
ml Anaconda3/5.3.0

# load conda env
source activate /groups/umcg-dag3/tmp01/rgacesa_tools/conda/envs/dag3pipe_v3_conda

# run clade profiling
strainphlan -s /groups/umcg-llnext/tmp01/pilot_microbiome/pilot_april_2022/strainphlan_all_april_2022/*.pkl --print_clades_only --output_dir . > LLNEXT_pilot_april_clades.txt

```
#### Execution ######

```
sbatch ./profileClades.sh 

```
# This will generate a .txt file with the list of clades that have been detected in the samples
# We process the output file to select only the clade names

```
cat LLnext_clades.txt | grep s__ | cut -f 2 | cut -f 1 -d ':' > LLnext_clades_names.txt

```
# Step 2: Perform MSA

Example: sbatch ./doMarkerComparisonLLNext.sh s__Bifidobacterium_bifidum

# Here, we use a loop to iterate over all detected clades
```
for i in $(cat LL_next_clades_names.txt); do sbatch doMarkerComparisonLLNext.sh $i; done 
```

#This will perform MSA and create .tre files and .aln files for each of the species you feed it in. 
# MSA was performed on consensus marker presence in at least in 50 samples 

#doMarkerComparisonLLNext.sh consists of: 

```
#!/bin/bash

#SBATCH --mem=32gb
#SBATCH --time=0-07:59:00
#SBATCH --cpus-per-task=8
#SBATCH --open-mode=truncate

# NOTES:
# > $1 is clade name
module purge
ml Anaconda3/5.3.0
# load conda env
source activate /groups/umcg-dag3/tmp01/rgacesa_tools/conda/envs/dag3pipe_v3_conda

mkdir ${1}
strainphlan -s /groups/umcg-dag3/tmp01/NEXT_pilot_results/strainphlan3/*.pkl --output_dir ./${1} --clade ${1} --marker_in_n_samples 50 --sample_with_n_markers 20 --nprocs 8
doMarkerComparisonLLNext.sh (END)

```
# Step 3: Make distance matrix from MSA file

#Example: 
#bash ./makeDistMat.sh ./s__Alistipes_shahii/s__Alistipes_shahii.StrainPhlAn3_concatenated.aln

```
for i in $(find . -type f -name *.aln); do bash makeDistMat.sh $i; done 
```
```
#!/bin/bash
echo 'maker of distance matrix from multiple alignment'
echo ' feed it with .aln file (multiple alignment)'
echo 'NOTE: make sure conda is loaded'

#Create a distance matrix from a multiple sequence alignment using the EMBOSS package (https://www.bioinformatics.nl/cgi-bin/emboss/help/distmat) 
# distmat calculates the evolutionary distance between every pair of sequences in a multiple sequence alignment.
# Uses Kimura Two-Parameter distance (distances expressed in terms of the number of substitutions per 100 b.p or amino acids) 
distmat -sequence ${1} -nucmethod 2 -outfile ${1/.aln/.dmat}

```
# Step 4: Cleaning the distance matrix from MSA file 

# First, load RPlus
ml RPlus 

#Example: Rscript parseDMat_LLNext.R s__Bifidobacterium_bifidum.dmat
```
for i in $(find . -type f -name *.dmat); do Rscript parseDMat_LLNext.R $i; done 
```
#parseDMat_LLNext.R consists of: 
```
library(optparse)
# CL PARSER
help_description <- ""
args <- OptionParser(usage = "%prog clade.distmat.txt metadata.txt ordination.png",
                      add_help_option = TRUE, prog='strainphlan_ordination.R',
                      description=help_description )
args_list <- parse_args(args, positional_arguments=TRUE)

# read in the file, skipping the first 8 rows and filling in empty columns, using the tab as sep, and stripping extra white space
inFile <- args_list$args[1]
data <- read.table( inFile, skip = 8, fill = TRUE, sep="\t", strip.white = T)

# remove the first column of the data as it is blank
data[1] <- NULL

# get the header as the last column of the data as a character vector
header <- lapply(data[,ncol(data)], as.character)

# remove the last column from the data as it has been stored as a header
data[ncol(data)] <- NULL

# remove the current last column from the data as it is blank
data[ncol(data)] <- NULL

# split header by space and digit to get the sample names
samples <- unlist(header)
# fix metaphlan just by taking only first 3 thingies after _ split
ss <- c()
sNR <- 0
for (s in samples) {
   sNR <- sNR + 1
   ssplit = strsplit(s,' ')
   ss <- c(ss,ssplit[[1]][1])
}

```

