
# Viruses in Negative Controls

This repository contains the scripts used during the data analysis for the paper **“Viruses in Negative Controls”**.

## Citation

The citation for this study will be provided after publication.

## Short Summary

This study aims to evaluate contamination patterns in viral-like particle (VLP) sequencing data by analyzing publicly available samples and negative controls (NCs) from various mother-infant cohorts. It also proposes a decontamination strategy to effectively remove contaminants from the data. Four studies were selected for this analysis:

1. **Garmaeva, Sanzhima et al.**  
   *“Transmission and dynamics of mother-infant gut viruses during pregnancy and early life.”*  
   Nature Communications, vol. 15, no. 1, 1945, 2 Mar. 2024.  
   [doi:10.1038/s41467-024-45257-4](https://doi.org/10.1038/s41467-024-45257-4)
   
2. **Liang, Guanxiang et al.**  
   *“The stepwise assembly of the neonatal virome is modulated by breastfeeding.”*  
   Nature, vol. 581, no. 7809, 2020, pp. 470-474.  
   [doi:10.1038/s41586-020-2192-1](https://doi.org/10.1038/s41586-020-2192-1)
   
3. **Maqsood, Rabia et al.**  
   *“Discordant transmission of bacteria and viruses from mothers to babies at birth.”*  
   Microbiome, vol. 7, no. 1, 156, 10 Dec. 2019.  
   [doi:10.1186/s40168-019-0766-7](https://doi.org/10.1186/s40168-019-0766-7)
   
4. **Shah, Shiraz A et al.**  
   *“Expanding known viral diversity in the healthy infant gut.”*  
   Nature Microbiology, vol. 8, no. 5, 2023, pp. 986-998.  
   [doi:10.1038/s41564-023-01345-7](https://doi.org/10.1038/s41564-023-01345-7)

Details regarding metadata and the study setups can be found in the original publications cited above.

## Scripts Overview

The scripts in this repository are divided into three categories:

### 1. Upstream Analysis
- **Location:** [Upstream Analysis Folder](https://github.com/GRONINGEN-MICROBIOME-CENTRE/Lifelines_NEXT/tree/main/NCP_VLP_project/Upstream_analysis)
- This section contains scripts used for raw data processing up to the creation of RPKM tables. The scripts are organized into four study-specific folders, named after the first authors of the publications used in this study. 
- Additionally, the folder named "NCP_all" includes scripts for pooled data analysis, including dereplication, decontamination, and RPKM table creation.
- Further information can be found in the README files provided within each folder.

### 2. Downstream Analysis
- **Location:** [Downstream Analysis Folder](https://github.com/GRONINGEN-MICROBIOME-CENTRE/Lifelines_NEXT/tree/main/NCP_VLP_project/Downstream_analysis)
- This section contains scripts used for statistical analysis and graph creation.

### 3. Metadata processing
- **Location:** [Metadata Folder](https://github.com/GRONINGEN-MICROBIOME-CENTRE/Lifelines_NEXT/tree/main/NCP_VLP_project/Metadata_processing)
- This section contains scripts used for metadata creation.
