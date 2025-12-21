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
  <img src="images/Project_overview.png" alt="Overview of pregnancy and early-life gut virome analysis workflow" width="900">
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
  - [Viral Detection](#viral-detection)
  - [Contig Extension & Filtering](#contig-extension--filtering)
  - [Quality Control & Host Removal](#quality-control--host-removal)
  - [Dereplication & Abundance Estimation](#dereplication--abundance-estimation)
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

### Viral Identification and abundance profiling
We used multiple complementary tools to identify viral sequences from assembled contigs:
- **VirSorter2**, **DeepVirFinder**, and **geNomad** for initial viral prediction.
- **COBRA** for extension of predicted viral contigs.
- **geNomad** for initial prunning of prophage regions and additional filtering of extended viral contigs.
- **CheckV** for host contamination removal (final prophage prunning) and viral completeness estimation/filtering.
- **CheckV** anicalc.py and aniclust.py scripts for dereplication of predcicted viral genomes. Viral genomed were dereplicated:
  (A) Internally within the dataset (deduplication using 99% ANI / 95% AF; dereplication using 95% ANI / 85% AF).
  (B) Against public viral genome databases (95% ANI / 85% AF). Databases used included: MGV, GPD, IMG/VR, ELGV, RefSeq, Shah et al, Benler et al., and CrAss-like phage databases (Gulyaeva et al., Yutin et al. and Guerin et al.)
- **Bowtie2** for abundance estimation via read mapping.

### Viral Characterization

  
---

## Data Availability
All resulting viral genome sets, abundance tables, and associated metadata will be made available through:
- [Zenodo DOI link (coming soon)](https://zenodo.org/)
- Lifelines (to be added)

---

## Code Structure
Scripts and workflows used in this project are organized as follows:

---

## Citation

