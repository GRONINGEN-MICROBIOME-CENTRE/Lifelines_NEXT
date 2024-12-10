# Demo data and code for Lifelines NEXT manuscript

Here is the emulated dataset and example codebase for the analysis used in the Lifelines NEXT manuscript that is currently on [ResearchSquare](https://www.researchsquare.com/article/rs-5334252/v1).

## Beta diversity and between-sample distances

The script [step1_dist_analyses.R](step1_dist_analyses.R) performs this type of analysis, including 
PERMANOVA analysis for the phenotypes, and comparing between-sample distances between groups. Mostly, 
the code uses ``` vegan ``` package

## Association analysis

The script [step2_associations.R](step2_associations.R) performs association between microbiome features (such as taxa) with 
infant/maternal phenotypes, including *quantitative* and *binary* variables, as well as *multilevel factors*. 
Note that phenotypes might be *static*, which don't change throughout timespan (e.g. delivery mode) or *dynamic*,
which are changing thoughout time (e.g. feeding, weight)


## Factor analysis 

The script [step3_MEFISTO.R](step3_MEFISTO.R) performs factor analysis of microbiome data and its association to XXX. We use package

## Strain transmission analysis

The script [step4_strainTransmission.R](step4_strainTransmission.R) ``` takes tree inputs from RAxML and performs 
the analysis on strain transmission. 

## PGLMM analysis

The script [step5_PGLMM.R](step5_PGLMM.R) performs Phylogenetic Generalized Linear Mixed Model 
analysis of association of phenotypes to strain phylogeny. NOTE: PGLMM doesn't use RAxML trees. Instead,
trees were rebuilt using ``` phangorn ``` package to build rooted ultrametric trees (what RAxML can't do) 
