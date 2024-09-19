
# Shah et al. (2023) Data Analysis Scripts

This repository contains scripts used for the analysis of sequencing data from the Shah et al. (2023) study.

## Table of Contents
1. [Sample rename](#sample-rename)
2. [Assembly](#assembly)
3. [Virome Discovery](#virome-discovery)
4. [Contig Extension](#contig-extension)
5. [Plasmid Removal, Prophage Pruning, and Quality Assessment](#plasmid-removal-prophage-pruning-and-quality-assessment)
6. [Read Alignment and Mapping](#read-alignment-and-mapping)
7. [Auxiliary Scripts](#auxiliary-scripts)

**Note**: No (QC) was performed on the samples, as it had already been completed by the authors of the original paper before data upload to the ENA.

---

## Sample rename
Script for renaming the sample. Performed due to the lack of consistency with the other studies, to siplify the further analysis.

- **sample_rename_shah.R**

---

## Assembly
Script for assembling the sequencing data.

- **02.sc_assembly.sh**  
  Assembly script for all samples and NCs.

### Assembly Details:
- Performed using `metaSPAdes` in meta mode.

---

## Virome Discovery
Discovery of viral sequences from assembled contigs larger than 1 kb using several tools.

- **03.vd_deepvirfinder.sh**  
  Virome discovery using DeepVirFinder.

- **03.vd_virsorter.sh**  
  Virome discovery using VirSorter.  
  **Note**: Download `prodigal-gv` prior to running this script (link in the script).

- **03.vd_vibrant.sh**  
  Virome discovery using VIBRANT.

- **03.vd_genomad.sh**  
  Virome discovery using geNomad.

---

## Parsing Virome Discovery Results
- **04.vd_parsing.sh**  
  Script for parsing the results from the virome discovery tools.

---

## Contig Extension
Extension of contigs using COBRA to create longer contigs.

- **05.vd_cobra.sh**  
  Contig extension script for all samples and NCs.

---

## Plasmid Removal, Prophage Pruning, and Quality Assessment
- **06.vd_prophage_pruning_fin.sh**  
  Script to remove plasmids (identified by geNomad), prune prophages, and perform quality assessment using Check-V.  
  **Note**: Automatically launches the R script `New_contigs_ID_and_metadata_fin.R`.

---

## Read Alignment and Mapping
### Read Alignment
- **07.pd_readalign.sh**  
  Script for aligning reads to the viral contigs.

### Read Mapping to Dereplicated Contigs
- **09.at_readmapping_NCP.sh**  
  Script for mapping reads to dereplicated viral contigs for samples.

### Read Mapping to Decontaminated Contigs
- **09.at_readmapping_decontam_99_NCP.sh**  
  Script for mapping reads to dereplicated viral contigs (99% breadth of coverage) for samples.

---

## Auxiliary Scripts
- **runAllSamples_02.bash**  
  Wrapper script for running `01.reads_QC.sh` and `02.sc_assembly.sh` for all samples.

- **New_contigs_ID_and_metadata_fin.R**  
  R script to create new unified IDs and contig metadata.  
  **Note**: This script is automatically launched when running `06.vd_prophage_pruning_fin.sh`.

- **filter_contigs.pl**  
  Script to filter out contigs shorter than 1 kb.  
  **Note**: This script is automatically launched within `02.sc_assembly.sh`.

- **table_of_origin.R**  
  R script to generate contig metadata after viral discovery.
