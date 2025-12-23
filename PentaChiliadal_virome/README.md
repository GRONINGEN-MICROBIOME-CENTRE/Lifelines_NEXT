# Pregnancy and Early-Life Gut Virome in the LLNEXT Cohort  
## Origin, Persistence, Influencing Factors, and Health Implications  

Gut virome analysis from metagenomic data in the **Lifelines NEXT (LLNEXT)** mother–infant cohort.

<p align="center">
  <a href="https://umcgresearch.org/w/lifelines-next">
    <img src="images/banner_lifelines_next.png" alt="Lifelines NEXT Cohort" height="55">
  </a>
  &nbsp;&nbsp;&nbsp;
  <a href="https://umcgresearchdatacatalogue.nl/all/cohorts/LIFELINES_NEXT">
    <img src="images/banner_data_access.png" alt="Lifelines NEXT Data Catalogue" height="55">
  </a>
  &nbsp;&nbsp;&nbsp;
  <a href="https://pubmed.ncbi.nlm.nih.gov/32100173/">
    <img src="images/banner_cohort_paper.png" alt="Lifelines NEXT Cohort Publication" height="55">
  </a>
</p>

<p align="center">
  <img src="images/Overview LLNEXT.png" alt="Overview of pregnancy and early-life gut virome analysis workflow" width="900">
</p>

---

## Overview
Here, we present a comprehensive **characterization of early-life gut virome assembly alongside maternal gut viral dynamics during pregnancy**. 
- 🧬 **Large-scale virome profiling:** We analyzed the DNA virome from 4,523 fecal and 91 breastmilk metagenomes across 714 mother–infant pairs in the LLNEXT cohort, generating a catalog of >31,000 viral operational taxonomic units (vOTUs).
- 🌱 **Contrasting virome dynamics:** The maternal gut virome is highly stable over time, whereas the infant gut virome undergoes rapid diversification during early life.
- 🚼 **Drivers of infant virome development:** Delivery mode and feeding mode are the primary determinants of infant virome developmental trajectories, with additional influences from maternal parity (presence of siblings).
- 🦠 **Virome and health outcomes:** Increased viral diversity in infancy is associated with the development of food allergy.
- 🔗 **Maternal origin of the infant virome:** Strain-level analyses identify the maternal gut as the dominant source of infant gut viruses (with reduced viral sharing following cesarean delivery) and lower sharing detected from breastmilk.
- 🛡️ **Genetic determinants of persistence:** The presence of DNA adenine N6-methyltransferase *hin1523* and diversity-generating retroelements promotes long-term viral persistence in the infant gut.

Together, these findings define the origin, dynamics, and modulating factors of the infant gut virome and uncover genetic strategies that support viral persistence in the gut ecosystem.

---

## Table of Contents
- [Cohort](#cohort)
- [Methods](#methods)
  - [Assembly](#assembly)
  - [Virus Identification and Abundance Profiling](#virus-identification-and-abundance-profiling)
  - [Virus Characterization](#viral-characterization)
  - [Strain-Level Profiling and Characterization of Shared Viruses](#strain-level-profiling-and-characterization-of-shared-viruses)
- [Data Availability](#data-availability)
- [Code Structure](#code-structure)
- [Citation](#citation)


---

## Cohort

This study is based on the **Lifelines NEXT (LLNEXT)** cohort, a large, prospective, population-based birth cohort in the Netherlands designed to investigate early-life determinants of health and disease.

- **Cohort overview and data catalogue:**  
  https://umcgresearchdatacatalogue.nl/all/cohorts/LIFELINES_NEXT

- **Original cohort description:**  
  https://pubmed.ncbi.nlm.nih.gov/32100173/

---

## Methods

### Virus Identification and Abundance Profiling

Multiple complementary tools were used to identify viral sequences from assembled contigs (**metaSPAdes** assembly):

- **VirSorter2**, **DeepVirFinder**, and **geNomad** for initial viral sequence prediction.
- **COBRA** for extension of predicted viral contigs.
- **geNomad** for initial pruning of prophage regions and additional filtering of extended viral contigs.
- **CheckV** for host contamination removal (final prophage pruning) and viral genome completeness estimation and filtering.
- **Dereplication** using CheckV scripts (`anicalc.py` and `aniclust.py`):
    - Within dataset: 99% ANI / 95% AF for deduplication, 95% ANI / 85% AF for dereplication
    - Against public viral genome databases (95% ANI / 85% AF): MGV, GPD, IMG/VR, ELGV, RefSeq, Shah et al., Benler et al., and CrAss-like phage datasets (Gulyaeva et al., Yutin et al., Guerin et al.)
- **Abundance estimation**: mapping quality-controlled reads to dereplicated viral genome catalog using **Bowtie2**

---

### Virus Characterization

- **Taxonomic assignment:** Hierarchical approach combining RefSeq annotations, **geNomad**, and **VITAP**.
- **Bacteriophage lifestyle prediction:** **BACPHLIP**.
- **Host prediction:** **iPHoP**, using (A) the default reference database and (B) a database supplemented with dereplicated metagenome-assembled genomes (MAGs) generated in this study.
- **Anti-defense system identification:** **DefenseFinder**.
- **Diversity-generating retroelement (DGR) detection:** DGR scripts described by Roux *et al.* (2021)  
  https://bitbucket.org/srouxjgi/dgr_scripts/src/master/
- **DGR activity analysis:** **Anvi’o** (`anvi-gen-variability-profile`).

---

### Strain-Level Profiling and Characterization of Shared Viruses

- **Strain-level viral profiling:** **inStrain**, with strain sharing defined as popANI ≥ 99.999%.
- **Phage–bacterial host co-sharing:** Mapping of temperate phages to MAGs using **minimap2**, followed by detection in maternal–infant MAG pairs shared at the strain level (**SKANI** ANI > 99.9%).
- **Viral protein clustering:** Two-step clustering using **MMseqs2**, following the UHGV framework.
- **Protein and protein-family functional annotation:**  
  (A) Representative proteins from protein families were annotated against HMM profiles from **PHROGs**, **KOfam**, and **(Anti)DefenseFinder**.  
  (B) Unannotated protein families were further analyzed using:
  - Structural prediction with **ColabFold** (using MSAs enriched with **UniRef30** sequences)
  - Structural similarity searches using **Foldseek** against the Protein Data Bank (PDB), AlphaFold Database (AFDB), and the Big Fantastic Virus Database (BFDV).

---

## Data Availability

The LLNEXT gut viral (LLNEXT-GV) and breastmilk (LLNEXT-BMV) catalogs are freely available via **Zenodo** (DOI to be added).

| File | Description | Link |
|------|-------------|------|
| `LLNEXT_GV_vOTU_representatives.fna.gz` | LLNEXT-GV vOTU representative viral genomes | [Download](https://zenodo.org/record/XXXXXXX/files/LLNEXT_GV_representatives.fna.gz) |
| `LLNEXT_GV_metadata.tsv` | Metadata for all LLNEXT-GV species-level vOTUs | [Download](https://zenodo.org/record/XXXXXXX/files/LLNEXT_GV_metadata.tsv) |
| `LLNEXT_GV_vOTU_representatives_LLNEXT_only.fna.gz` | LLNEXT-GV vOTU representative viral genomes (excluding external DBs) | [Download](https://zenodo.org/record/XXXXXXX/files/LLNEXT_GV_representatives.fna.gz) |
| `LLNEXT_GV_viral_genomes_LLNEXT_only.fna.gz` | LLNEXT-GV deduplicated viral genomes (excluding external DBs) | [Download](https://zenodo.org/record/XXXXXXX/files/LLNEXT_GV_representatives.fna.gz) |
| `LLNEXT_BMV_representatives.fna.gz` | LLNEXT-BMV representative viral genomes (breastmilk metagenomes) | [Download](https://zenodo.org/record/XXXXXXX/files/LLNEXT_BMV_representatives.fna.gz) |
| `LLNEXT_BMV_metadata.tsv` | Metadata for all LLNEXT-BMV species-level vOTUs | [Download](https://zenodo.org/record/XXXXXXX/files/LLNEXT_BMV_metadata.tsv) |

---

## Code Structure
All scripts used in this project can be found within scripts and are organized as follows:

```text
scripts/
├── A_Virus_Identification/
│   ├── Virus_detection/                  # VirSorter2, DeepVirFinder, geNomad
│   ├── Virus_contig_extension/           # COBRA-based viral contig extension
│   ├── Virus_filtering/                  # geNomad filtering and prophage pruning
│   └── Contamination_quality_check/      # Quality control
│
├── B_Dereplication/                      # BLAST DB generation and ANI-based dereplication (aniclust.py, anicalc.py)
│
├── C_Abundance_estimation/               # Viral abundance estimation (read mapping)
│
├── D_Virus_Characterization/
│   ├── Genome_quality/                   # Viral genome completeness estimation
│   ├── ORF_prediction_and_clustering/    # ORF prediction (Prodigal-gv) and genome clustering (AAI)
│   ├── Taxonomy/                         # Viral taxonomy (geNomad, ViTAP)
│   ├── Host_and_lifestyle_prediction/    # Host prediction (iPHoP) and lifestyle (BACPHLIP)
│   └── Functional_profiling/
│       ├── Anti_defense/                 # Anti-defense system detection (AntiDefenseFinder)
│       └── DGRs/                         # DGR detection and activity
│
├── E_Strain_sharing/
│   ├── Strain_level_profiling/            # Strain-level viral profiling (inStrain)
│   ├── Virus_bacteria_co_sharing/         # Virus–MAG mapping and MAG ANI comparison
│   └── Characterization_enriched_proteins/
│       ├── Protein_clustering/            # MMseqs2 two-step clustering
│       ├── Sequence_based_annotation/     # HMM matching (AntiDefenseFinder, PHROGs, KOfam)
│       └── Structure_based_annotation/    # Structure prediction (ColabFold) and Foldseek search
│
└── Downstream_scripts/                   # Statistical analyses and figure generation

```
---

## Citation

If you use this resource in your research, please cite both the publication and the associated data resource.

### Publication

> **Pregnancy and Early-Life Gut Virome in the Lifelines NEXT Cohort: Origin, Persistence, Influencing Factors, and Health Implications**  
> Author list to be added.  
> *Journal* (Year).

---

### Data resource

> Author list to be added. (Year). **Pregnancy and Early-Life Gut Virome in the Lifelines NEXT Cohort** [Data set]. Zenodo.  
> https://doi.org/10.5281/zenodo.XXXXXXX

