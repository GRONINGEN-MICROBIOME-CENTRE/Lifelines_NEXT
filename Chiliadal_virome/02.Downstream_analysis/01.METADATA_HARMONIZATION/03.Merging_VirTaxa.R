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

genus_clusters <- read.table('06.CLEAN_DATA/genus_clusters_127553vOTUr_labeled_long_format.txt', sep = '\t', header=T)
genus_size <- read.table('06.CLEAN_DATA/genus_clusters_127553vOTUr_labeled_size.txt', sep='\t', header=T)

family_clusters <- read.table('06.CLEAN_DATA/family_clusters_127553vOTUr_labeled_long_format.txt', sep='\t', header = T)
family_size <- read.table('06.CLEAN_DATA/family_clusters_127553vOTUr_labeled_size.txt', sep='\t', header=T)
#############################################################
# 3.1 Analysis (merging)
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
# 3.2 Analysis (adding lost genera and families)
#############################################################
# genera:
N_gen_raw <- dim(genus_size)[1]

lost_genera <- as.data.frame(setdiff(ETOF_vOTUr$New_CID, genus_clusters$Cluster_member))
colnames(lost_genera) <- "Cluster_member"

N_gen_lost <- dim(lost_genera)[1]

lost_genera$Representative <- paste0("Genus_", (N_gen_raw + 1) : (N_gen_raw + N_gen_lost))

genus_clusters_UPD <- rbind(genus_clusters, lost_genera)

genus_size_UPD <- as.data.frame(table(genus_clusters_UPD$Representative))
colnames(genus_size_UPD) <- colnames(genus_size)

rm(genus_clusters, genus_size)
# save new output?

# families:
N_fam_raw <- dim(family_size)[1]

lost_families <- as.data.frame(setdiff(ETOF_vOTUr$New_CID, family_clusters$Cluster_member))
colnames(lost_families) <- "Cluster_member"

N_fam_lost <- dim(lost_families)[1]

lost_families$Representative <- paste0("Family_", (N_fam_raw + 1) : (N_fam_raw + N_fam_lost))

family_clusters_UPD <- rbind(family_clusters, lost_families)
family_size_UPD <- as.data.frame(table(family_clusters_UPD$Representative))
colnames(family_size_UPD) <- colnames(family_size)

rm(family_clusters, family_size)
# save new output?

#############################################################
# 2. Analysis: calculating novel genera
#############################################################
merged2<- vOTU_by_source %>%
  
  left_join(vOTU_cluster_size, by = "Representative") %>%
  
  rename(N_genomes = Cluster_size, vOTU_representative  = Representative) %>%
  
  left_join(genus_clusters_UPD, by = c("vOTU_representative" = "Cluster_member")) %>%
  
  rename(Genus = Representative) %>%
  
  left_join(ETOF_vOTUr %>% select(New_CID, POST_CHV_length, miuvig_quality, vOTU_cluster_type),
            by = c("vOTU_representative" = "New_CID")) %>%
  
  left_join(family_clusters_UPD,  by = c("vOTU_representative" = "Cluster_member")) %>%
  
  rename(Family = Representative) %>%
  
  mutate(
    vOTU_cluster_type = if_else(vOTU_cluster_type %in% c('NEXT_MGS', 'NEXT_VLP', 'NEXT_VLP+NEXT_MGS'), "NEXT", "Mixed"),
    miuvig_quality = factor(miuvig_quality,
                            levels = c("High-quality", "Genome-fragment"),
                            ordered = TRUE)
  )

genus_family_check <- merged2 %>%
  distinct(Genus, Family) %>%
  count(Genus) %>%
  filter(n > 1)

major_family <- merged2 %>%
  group_by(Genus, Family) %>%
  tally(name = "count") %>%
  slice_max(order_by = count, n = 1, with_ties = FALSE) %>%
  ungroup()

merged2_fixed <- merged2 %>%
  select(-Family) %>%            
  left_join(major_family, by = "Genus")

merged2 <- merged2_fixed

table(merged2$vOTU_cluster_type) # 110,796 NEXT-only, 16,757 - Mixed or other DB only

table(merged2[merged2$vOTU_representative %in% hq_votus,]$vOTU_cluster_type) # 8,580 NEXT-only, 7628 - Mixed or other DB only


top_gen <- merged2 %>%
  
  select(Genus, vOTU_representative, POST_CHV_length, miuvig_quality) %>%
  
  group_by(Genus) %>%
  
  arrange(desc(POST_CHV_length), desc(miuvig_quality)) %>%
  
  slice_head(n = 1) %>%
  
  ungroup() %>%
  
  left_join(
    merged2 %>% group_by(Genus) %>% summarise(across(c(2:15), sum), .groups="drop"),
    by = "Genus"
  ) %>%
  
  left_join(genus_size_UPD, by = c("Genus" = "Representative")) %>%
  
  rename(N_vOTUs = Cluster_size) %>%
  
  mutate(
    genus_cluster_type = if_else(N_genomes == NEXT_MGS + NEXT_VLP, "NEXT", "Mixed"),
    N_vOTUs = as.numeric(N_vOTUs)
  )

table(top_gen$genus_cluster_type) # 24,800 NEXT-only, 3783 - Mixed or other DB only

table(top_gen[top_gen$vOTU_representative %in% hq_votus,]$genus_cluster_type) # 494 NEXT-only, 1,260 - Mixed or other DB only

check <- top_gen[ top_gen$genus_cluster_type=="NEXT" & top_gen$miuvig_quality=="High-quality",]
sum(check$N_vOTUs) # 1,614
sum(check$N_vOTUs > 1)


# family
top_fam <- merged2 %>%
  
  select(Family, vOTU_representative, POST_CHV_length, miuvig_quality) %>%
  
  group_by(Family) %>%
  
  arrange(desc(POST_CHV_length), desc(miuvig_quality)) %>%
  
  slice_head(n = 1) %>%
  
  ungroup() %>%
  
  left_join(
    merged2 %>% group_by(Family) %>% summarise(across(c(2:15), sum), .groups="drop"),
    by = "Family"
  ) %>%
  
  left_join(family_size_UPD, by = c("Family" = "Representative")) %>%
  
  rename(N_vOTUs = Cluster_size) %>%
  
  mutate(
    family_cluster_type = if_else(N_genomes == NEXT_MGS + NEXT_VLP, "NEXT", "Mixed"),
    N_vOTUs = as.numeric(N_vOTUs)
  )

table(top_fam$family_cluster_type) # 6,954 NEXT-only, 952 - Mixed or other DB only

table(top_fam[top_fam$vOTU_representative %in% hq_votus,]$family_cluster_type) # 148 NEXT-only, 276 - Mixed or other DB only

check2 <- top_fam[top_fam$family_cluster_type=="NEXT" & top_fam$miuvig_quality=="High-quality",]
#sum(check2$N_vOTUs) ### here this won't be representative because the size of families changed after reassignment

check3 <- merged2_fixed[merged2_fixed$Family %in% check2$Family & merged2_fixed$miuvig_quality=="High-quality",] 
#293 novel vOTU make them up?
length(unique(check3$Genus)) # 177 genera make up 148 novel families

# creating genus_RPKM based on the vOTU binning into genus clusters 
group <- merged2_fixed$Genus[ match(rownames(cleanest),
                                    merged2_fixed$vOTU_representative) ]

genus_RPKM <- as.data.frame(rowsum(as.matrix(cleanest), group))

write.table(genus_RPKM, "06.CLEAN_DATA/02.FINAL/genus_RPKM_VLP_MGS.txt", sep='\t', quote=F)

group <- merged2_fixed$Family[ match(rownames(cleanest),
                                     merged2_fixed$vOTU_representative) ]

family_RPKM <- as.data.frame(rowsum(as.matrix(cleanest), group))

write.table(family_RPKM, "06.CLEAN_DATA/02.FINAL/family_RPKM_VLP_MGS.txt", sep='\t', quote=F)

#############################################################
# 4. Output
#############################################################
write.table(newtax[,add_t_ETOF], "./06.CLEAN_DATA/MergedTaxonomy_127553vOTUr_ab3kbp_in_2200_VLP_MGS.txt", sep='\t', quote=F, row.names=F)

