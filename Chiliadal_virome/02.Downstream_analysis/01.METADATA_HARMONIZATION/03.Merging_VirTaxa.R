setwd('~/Desktop/Projects_2021/NEXT_virome/09.DATA_ANALYSIS/')

#############################################################
# Merging taxonomy assignment by geNomad and VITAP
#############################################################

#############################################################
# 0. Used files source
#############################################################
# ETOF_127553vOTUr_ab3kbp_in_2200_VLP_MGS.txt - Filtered 
# expanded table of origin for 127,553 vOTUs detected in 2,220
# MGS and VLP samples; taxonomy is assigned by geNomad

# best_determined_lineages.tsv - VITAP run default raw output 
# for 127,553 vOTUs

# ICTV_Master_Species_List_2024_MSL40.v2.xlsx - 
# downloaded the current release from https://ictv.global/taxonomy
# https://ictv.global/msl/current

#############################################################
# 1. Functions
#############################################################

#############################################################
# 1. Loading libraries
#############################################################
library(tidyverse)

#############################################################
# 2. Load Input Data
#############################################################
ETOF_vOTUr <- read.table('06.CLEAN_DATA/ETOF_127553vOTUr_ab3kbp_in_2200_VLP_MGS.txt', sep='\t', header=T)

taxonomy <- read.table('06.CLEAN_DATA/best_determined_lineages.tsv', sep='\t', header=T)

ictv <- readxl::read_xlsx('06.CLEAN_DATA/ICTV_Master_Species_List_2024_MSL40.v2.xlsx', sheet=2)

#############################################################
# 3. Analysis (merging)
#############################################################
taxonomy <- taxonomy %>%
  separate(
    col    = lineage,
    into   = paste0("level", 1:8),  
    sep    = ";",
    fill   = "right",               
    remove = FALSE                  
  )

colnames(taxonomy)[3:10] <- c('Species', 'Genus', 'Family', 'Order', 'Class', 'Phylum', 'Kingdom', 'Realm')

taxonomy[taxonomy$Realm=="[Realm]_Anelloviridae",]$Realm <- "Monodnaviria"

# preparing geNomad assignments:
ETOF_vOTUr[is.na(ETOF_vOTUr$taxonomy),]$taxonomy <- "Unclassified"

genomad_taxa <- ETOF_vOTUr[,c("New_CID", "taxonomy")]
genomad_taxa <- genomad_taxa %>%
  separate(
    col    = taxonomy,
    into   = paste0("level", 1:7),  
    sep    = ";",
    fill   = "right",               
    remove = FALSE                  
  )

# removing a useless column
genomad_taxa$level1 <- NULL

colnames(genomad_taxa)[3:8] <- c('RealmG', 'KingdomG', 'PhylumG', 'ClassG', 'OrderG', 'FamilyG')
genomad_taxa[is.na(genomad_taxa$RealmG), ]$RealmG <- 'Unclassified'

genomad_taxa[genomad_taxa$taxonomy=="Viruses;;;;;;Bicaudaviridae",]$RealmG <- "Unclassified"
genomad_taxa[genomad_taxa$taxonomy=="Viruses;;;;;;",]$RealmG <- "Unclassified"
genomad_taxa[genomad_taxa$taxonomy=="Viruses;;;;Naldaviricetes;Lefavirales;Baculoviridae",]$RealmG <- "Unclassified"

genomad_taxa[genomad_taxa$RealmG=="Anelloviridae",]$FamilyG <- "Anelloviridae"
genomad_taxa[is.na(genomad_taxa$FamilyG),]$FamilyG <- 'Unclassified'
genomad_taxa[genomad_taxa$FamilyG=="Anelloviridae",]$RealmG <- "Monodnaviria"

genomad_taxa[genomad_taxa$RealmG=="Fuselloviridae",]$FamilyG <- "Fuselloviridae"
genomad_taxa[genomad_taxa$FamilyG=="Fuselloviridae",]$RealmG <- "Unclassified"

### used to see if geNomad and VITAP assignments overlap:
# test_overlap <- merge(taxonomy2[,c("Genome_ID", "Realm")], genomad_taxa[,c("New_CID", "RealmG")],
#                       by.x="Genome_ID", by.y="New_CID", all.y=T)
# test_overlap[is.na(test_overlap$Realm),]$Realm <- 'Unclassified'
# test_overlap$mismatch <- test_overlap$Realm != test_overlap$RealmG
#length(test_overlap[test_overlap$Realm!="Unclassified" & 
#                      test_overlap$RealmG=="Unclassified",]$Genome_ID)
#length(test_overlap[test_overlap$Realm=="Unclassified" & 
#                      test_overlap$RealmG!="Unclassified",]$Genome_ID)

# View(test_overlap[test_overlap$mismatch==T & test_overlap$RealmG!="Unclassified",])

tax_overlap <- merge(taxonomy, genomad_taxa, 
                     by.x = "Genome_ID", by.y="New_CID", all.y=T)

tax_overlap[is.na(tax_overlap$lineage), c('Species', 'Genus', 'Family', 'Order', 'Class', 'Phylum', 'Kingdom', 'Realm')] <- '-'
tax_overlap[is.na(tax_overlap$lineage), "lineage"] <- 'Unclassified'

# adding missing tax assignments using geNomad assignments

levels <- c("Realm","Kingdom","Phylum","Class","Order","Family")

for (nm in levels) {
  
  g <- paste0(nm, 'G')
  
  tax_overlap[[g]] <- as.character(tax_overlap[[g]])
  
  x <- tax_overlap[[g]]
  
  x[is.na(x) | x == ""] <- "Unclassified" 
  
  sel <- tax_overlap$lineage == "Unclassified" & x != "Unclassified"
  
  tax_overlap[[nm]][sel] <- x[sel]
}

newtax <- tax_overlap[,colnames(tax_overlap) %in% colnames(taxonomy)]

newtax$new_lineage <- paste0(newtax$Species, ';',
                             newtax$Genus, ';',
                             newtax$Family, ';',
                             newtax$Order, ';',
                             newtax$Class, ';',
                             newtax$Phylum, ';',
                             newtax$Kingdom, ';',
                             newtax$Realm, ';')

newtax[newtax$new_lineage=="-;-;-;-;-;-;-;-;",]$new_lineage <- "Unclassified"
newtax[is.na(newtax$Confidence_level),"Confidence_level"] <- 'geNomad'

newtax$genome <- ictv$Genome[match(newtax$Realm, ictv$Realm)]

newtax[newtax$Realm=="[Realm]_Bicaudaviridae",]$genome <- "dsDNA"
newtax[newtax$Realm=="[Realm]_Cedratvirus A11",]$genome <- "dsDNA"
newtax[newtax$Realm=="[Realm]_Lake Sarah-associated circular virus-14",]$genome <- "ssDNA"
newtax[newtax$Realm=="[Realm]_Lake Sarah-associated circular virus-29",]$genome <- "ssDNA"
newtax[newtax$Realm=="[Realm]_Lake Sarah-associated circular virus-42",]$genome <- "ssDNA"
newtax[newtax$Realm=="[Realm]_Mollivirus sibericum",]$genome <- "dsDNA"
newtax[newtax$Realm=="[Realm]_Mycoplasma phage phiMFV1",]$genome <- "dsDNA"
newtax[newtax$Realm=="[Realm]_Naldaviricetes",]$genome <- "dsDNA"
newtax[newtax$Realm=="[Realm]_Pacmanvirus A23",]$genome <- "dsDNA"
newtax[newtax$Realm=="[Realm]_Pandoravirus",]$genome <- "dsDNA"
newtax[newtax$Realm=="[Realm]_Pithovirus",]$genome <- "dsDNA"
newtax[newtax$Realm=="[Realm]_Polydnaviriformidae",]$genome <- "dsDNA"

newtax$genome_simple <- newtax$genome
newtax[is.na(newtax$genome_simple),]$genome_simple <- "Unknown"
newtax[newtax$genome_simple=="ssRNA(+/-)",]$genome_simple <- "RNA"
newtax[newtax$genome_simple=="ssDNA(+/-)",]$genome_simple <- "ssDNA"

add_t_ETOF <- c("Genome_ID", "Species", "Genus", "Family",
                "Order", "Class", "Phylum", "Kingdom", "Realm",
                "lineage_score", "Confidence_level", "new_lineage", "genome_simple")

#############################################################
# 4. Output
#############################################################
write.table(newtax[,add_t_ETOF], "./06.CLEAN_DATA/MergedTaxonomy_127553vOTUr_ab3kbp_in_2200_VLP_MGS.txt", sep='\t', quote=F, row.names=F)

