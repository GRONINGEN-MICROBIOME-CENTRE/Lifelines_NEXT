setwd('~/Desktop/Projects_2022/NEXT_pilot_FUP/')

#############################################################
# Here we create the metadata for virome by merging phenotypes
# and upstream analysis derived metrics. 
#############################################################

##############################
# Functions
##############################
# RPKM_table function requires two dfs of the certain format: 
# counts: N of reads aligned to the sequence (contig)
# coverage: N of non-zero bases of the sequence length
# Rows are contigs and columns are samples, and one column is 
# length of contigs. RPKM_table returns a table of RPKM counts. 
# Be aware that obtained RPKM counts will be representing counts 
# normalized on the number of ALIGNED reads (to the database), 
# not the initial size of quality-trimmed library. For the virome 
# analysis, please, remember, that quite some reads could be thrown 
# away after checking of coverage, i.e. RPKM transformation function 
# will NOT use "Total_reads" (output of samtools indexing).

RPKM_table <- function(counts, coverage){
  counts$coverage <- NULL # removing column with k-mer coverage of a contig
  counts <- counts[-c( (nrow(counts)-1) , nrow(counts) ),] # removing lines with total counts (as they also contain excluded contigs)
  counts <- counts[order( row.names(counts) ),]
  
  counts <- counts[row.names(counts) %in% row.names(coverage),]
  
  coverage <- coverage[order( row.names(coverage) ),]
  
  if (identical(row.names(coverage), row.names(counts))) {
    coverage$length <- counts$length 
  } else {
    print("Counts and coverage tables have different sets of viruses")
  }
  
  coverage <- coverage/coverage$length #normalizing the length coverage
  if (identical(colnames(counts), colnames(coverage))) {
    counts[coverage<=0.75] <-0
  } else {
    print("Counts and coverage tables have different sets of samples")
  }
  
  abundance_table <- as.data.frame( t( as.data.frame( t(counts) ) / (colSums(counts)/1000000) ) ) / (counts$length/1000)
  abundance_table$length <- NULL
  return(abundance_table)
}

##############################
# Loading libraries
##############################
library(dplyr)
library(stringr)
library(vegan)
##############################
# Input data & data transformation
##############################
VLP_metadata <- read.table('01.RAW_DATA/Metadata/VLP_metadata_updated_ages_01_02_2023.txt', sep='\t', header=T)

# based on the Qubit detection threshold & Novogene concentration assignment for too lows:
VLP_metadata[VLP_metadata$DNA_CONC=='too low',]$DNA_CONC <- "0.01"
VLP_metadata$DNA_CONC <- as.numeric(VLP_metadata$DNA_CONC)

# Number of reads per sample:
kneaddata_stat <- read.table('01.RAW_DATA/NEXT_pilot_kneaddata.txt', sep = ',', header=T, check.names = F)
colnames(kneaddata_stat)[1] <- 'Unified_ID'
kneaddata_stat$Raw_reads <- kneaddata_stat$`Raw reads (p1)` + kneaddata_stat$`Raw reads (p2)`
kneaddata_stat$Clean_reads <- kneaddata_stat$`Final (p1)` + kneaddata_stat$`Final (p2)`
kneaddata_stat$Human_reads <- kneaddata_stat$`Contaminants (p1)` + kneaddata_stat$`Contaminants (p2)`
kneaddata_stat <- kneaddata_stat[,c("Unified_ID","Raw_reads", "Clean_reads", "Human_reads")]

# N of scaffolds per sample predicted to be viral (before decontamination)
N_viral_contigs <- read.table('01.RAW_DATA/N_predicted_viral_scaffolds', sep=' ')
colnames(N_viral_contigs) <- c("N_viral_contigs", "Unified_ID")

# Length of all viral contigs per sample:
all_viral_length <- read.table('01.RAW_DATA/Length_contigs_all_viral_df', sep='\t')
colnames(all_viral_length) <- c('Unified_ID', 'Length_all_viral_contigs')

# Quast ouput (N of contigs assembled)
contigs_stat <- read.table('01.RAW_DATA/NEXT_pilot_quast_tmp.txt', sep='', header=T)
colnames(contigs_stat)[1] <- 'Unified_ID'
contigs_stat <- contigs_stat[,c("Unified_ID", "contigs...0bp.", "contigs...1000bp.",
                                "Total_length...1000bp.","N50","L50")]

# cpn60db alignment to bacterial genomic contamination
cpn60db_alignment <- read.table('01.RAW_DATA/cpn60_alignment_stat.txt', sep='\t', header=T)
colnames(cpn60db_alignment)[1] <- 'Unified_ID'
cpn60db_alignment <- merge(cpn60db_alignment, kneaddata_stat[,c("Unified_ID", "Clean_reads")], by='Unified_ID')
cpn60db_alignment$cpn60_fraction <- ((cpn60db_alignment$concordantly + cpn60db_alignment$concordantly1 + cpn60db_alignment$discordantly)*2 + cpn60db_alignment$single + cpn60db_alignment$single1)/(cpn60db_alignment$Clean_reads)*100
cpn60db_alignment$total_reads_aligned_cpn60db <- ((cpn60db_alignment$concordantly + cpn60db_alignment$concordantly1 + cpn60db_alignment$discordantly)*2 + cpn60db_alignment$single + cpn60db_alignment$single1)
cpn60db_alignment$bacterial_contamination_perc_reads <- cpn60db_alignment$cpn60_fraction*3500/0.555 #fraction x typical bacterial genome size in kb / size of cpn60 gene in kb

# read alignment to viral contigs
alignment_stat <- read.table('01.RAW_DATA/Counts_tables_final_noneg405_99_der95/alignment_stat.txt', sep='\t', header=T)
colnames(alignment_stat)[1] <- 'Unified_ID'

# merging to one metadata
metadata_full <- VLP_metadata %>%
  full_join(N_viral_contigs, by='Unified_ID') %>%
  full_join(all_viral_length, by='Unified_ID') %>%
  full_join(contigs_stat, by='Unified_ID') %>%
  full_join(kneaddata_stat, by='Unified_ID') %>%
  full_join(cpn60db_alignment[,c('Unified_ID','bacterial_contamination_perc_reads')], by='Unified_ID') %>%
  full_join(alignment_stat, by='Unified_ID')

# Excluding samples (Negative control, that was used for cleaning, MGS swap, mothers P3 and repeated measurement)
metadata_full <- metadata_full[ !(metadata_full$Unified_ID %in% metadata_full[ (metadata_full$TIMEPOINT=="P3" | 
                                                                                metadata_full$TIMEPOINT=="Blank" | 
                                                                                metadata_full$NUMBER==454 ) ,]$Unified_ID ) , ]
# counting number of samples per individual:
metadata_full <- metadata_full %>%
  add_count(NEXT_ID, name = 'N_timepoints')

# changing NG ids to improve the error with twins:
metadata_full[metadata_full$Lifelines_ID=="6527985", "ID_1"] <- "D"
metadata_full[metadata_full$Lifelines_ID=="6527985", "NG_ID"] <- "D01VLP04224E09"

# Short sample & Individual ID:
metadata_full$Short_sample_ID <- paste0(metadata_full$ID_1, 
                                      metadata_full$ID_02,
                                      str_pad(metadata_full$ID_4, 4, pad = "0"), 
                                      "V")
metadata_full$Individual_ID <- paste0(metadata_full$ID_1,
                                        str_pad(metadata_full$ID_4, 4, pad = "0"))

metadata_full$Universal_fecal_ID <- paste0(metadata_full$ID_1, 
                                           metadata_full$ID_02,
                                           str_pad(metadata_full$ID_4, 4, pad = "0") )

# adding continuous timepoint:
metadata_full$Timepoint_continuous <- NA
metadata_full[metadata_full$TIMEPOINT=="P7",]$Timepoint_continuous <- -2
metadata_full[metadata_full$TIMEPOINT=="B",]$Timepoint_continuous <- 0
metadata_full[metadata_full$TIMEPOINT=="M1",]$Timepoint_continuous <- 1
metadata_full[metadata_full$TIMEPOINT=="M2",]$Timepoint_continuous <- 2
metadata_full[metadata_full$TIMEPOINT=="M3",]$Timepoint_continuous <- 3
metadata_full[metadata_full$TIMEPOINT=="M6",]$Timepoint_continuous <- 6
metadata_full[metadata_full$TIMEPOINT=="M12",]$Timepoint_continuous <- 12

# Creating RPKM table for viruses in metaviromes
VLP_vircounts <- read.table("01.RAW_DATA/Counts_tables_final_noneg405_99_der95/bowtie2_read_counts_with_full.txt", sep='\t', header = T, row.names = 1)

VLP_contig_coverage <- read.table('01.RAW_DATA/Counts_tables_final_noneg405_99_der95/coverage_table.txt', sep='\t', header = T, row.names = 1)

RPKM_counts_VLP <- RPKM_table(VLP_vircounts, VLP_contig_coverage)

# removing excluded samples from RPKM table of viruses in metaviromes:
RPKM_counts_VLP <- RPKM_counts_VLP[,colnames(RPKM_counts_VLP) %in% metadata_full$Unified_ID]

RPKM_counts_VLP <- RPKM_counts_VLP[rowSums(RPKM_counts_VLP)>0,] #110 contigs did not survive dereplication

RPKM_counts_VLP <- RPKM_counts_VLP[,metadata_full$Unified_ID]

if (identical(metadata_full$Unified_ID, colnames(RPKM_counts_VLP))) {
  colnames(RPKM_counts_VLP) <- metadata_full$Short_sample_ID
}
colSums(RPKM_counts_VLP[,c('D010422V', 'C010422V')]>0)

# Total MGS metadata
MGS_metadata <- read.table("01.RAW_DATA/Metadata/MGS_metadata_updated_ages_01_02_2023.txt", sep='\t', header=T)
# Short sample & Individual ID:
MGS_metadata$Short_sample_ID <- paste0(substr(MGS_metadata$NG_ID, 1,3),
                                       substr(MGS_metadata$NG_ID, 5,8),
                                        "M")
MGS_metadata$Short_sample_ID_bact <- paste0(substr(MGS_metadata$NG_ID, 1,3),
                                       substr(MGS_metadata$NG_ID, 5,8),
                                       "B")
MGS_metadata$Universal_fecal_ID <- paste0(substr(MGS_metadata$NG_ID, 1,3),
                                            substr(MGS_metadata$NG_ID, 5,8) )
MGS_metadata$Individual_ID <- paste0(substr(MGS_metadata$NG_ID, 1,1),
                                     substr(MGS_metadata$NG_ID, 5,8))

# adding continuous timepoint:
MGS_metadata$Timepoint_continuous <- NA
MGS_metadata[MGS_metadata$Timepoint=="P3",]$Timepoint_continuous <- -6
MGS_metadata[MGS_metadata$Timepoint=="P7",]$Timepoint_continuous <- -2
MGS_metadata[MGS_metadata$Timepoint=="MB",]$Timepoint_continuous <- 0
MGS_metadata[MGS_metadata$Timepoint=="M1",]$Timepoint_continuous <- 1
MGS_metadata[MGS_metadata$Timepoint=="M2",]$Timepoint_continuous <- 2
MGS_metadata[MGS_metadata$Timepoint=="M3",]$Timepoint_continuous <- 3
MGS_metadata[MGS_metadata$Timepoint=="M6",]$Timepoint_continuous <- 6
MGS_metadata[MGS_metadata$Timepoint=="M9",]$Timepoint_continuous <- 9
MGS_metadata[MGS_metadata$Timepoint=="M12",]$Timepoint_continuous <- 12

# counting number of samples per individual:
MGS_metadata <- MGS_metadata %>%
  add_count(NEXT_ID, name = 'N_timepoints')

# % reads aligned to the virus references:
alignment_stat_MGS <- read.table('01.RAW_DATA/Counts_tables_MGS/alignment_stat.txt', sep='\t', header=T)
colnames(alignment_stat_MGS)[1] <- "NG_ID"

MGS_metadata <- MGS_metadata %>%
  left_join(alignment_stat_MGS, by='NG_ID')


# Creating RPKM table for viruses in metagenomes
MGS_vircounts <- read.table("01.RAW_DATA/Counts_tables_MGS/bowtie2_read_counts_with_full.txt", sep='\t', header = T, row.names = 1)

MGS_contig_coverage <- read.table('01.RAW_DATA/Counts_tables_MGS/coverage_table.txt', sep='\t', header = T, row.names = 1)

RPKM_counts_MGS <- RPKM_table(MGS_vircounts, MGS_contig_coverage)

RPKM_counts_MGS <- RPKM_counts_MGS[,colnames(RPKM_counts_MGS) %in% MGS_metadata$NG_ID]

RPKM_counts_MGS <- RPKM_counts_MGS[rowSums(RPKM_counts_MGS)>0,] # 23,695 contigs are not detected in total MGS samples

RPKM_counts_MGS <- RPKM_counts_MGS[,MGS_metadata$NG_ID]

if (identical(MGS_metadata$NG_ID, colnames(RPKM_counts_MGS))) {
  colnames(RPKM_counts_MGS) <- MGS_metadata$Short_sample_ID
}

RPKM_counts_MGS <- RPKM_counts_MGS[row.names(RPKM_counts_MGS) %in% row.names(RPKM_counts_VLP), ]

# Viral contigs metadata
viral_contigs_metadata <- read.table('01.RAW_DATA/table_of_origin_decontaminated', sep = '\t', header=T)

viral_taxonomy <- read.table('01.RAW_DATA/DemoVir_assignments.txt', sep='\t', header=T)
colnames(viral_taxonomy)[1] <- "V1"
viral_contigs_metadata <- merge(viral_contigs_metadata, viral_taxonomy, by="V1", all.x = T)

viral_contigs_metadata$Family[viral_contigs_metadata$CrAss==1] <- "CrAss-like"
viral_contigs_metadata$Order[viral_contigs_metadata$Family=="CrAss-like"] <- "Crassvirales"
viral_contigs_metadata$Family[is.na(viral_contigs_metadata$Family)] <- "Unassigned"
viral_contigs_metadata$Order[is.na(viral_contigs_metadata$Order)] <- "Unassigned"

excluded_order <- c("no_order_Ascoviridae", "no_order_Astroviridae", "no_order_Caliciviridae",
                    "no_order_Flaviviridae", "no_order_Iridoviridae", "no_order_Marseilleviridae",
                    "no_order_Mimiviridae", "no_order_Nudiviridae", "no_order_Phycodnaviridae",  
                    "no_order_Pithoviridae", "no_order_Poxviridae", "no_order_Spiraviridae")
excluded_family <- c("Ascoviridae", "Astroviridae", "Caliciviridae","Flaviviridae", 
                     "Iridoviridae", "Marseilleviridae", "Mimiviridae", "Nudiviridae", 
                     "Phycodnaviridae", "Pithoviridae", "Poxviridae", "Spiraviridae")
for (i in excluded_order) {
  viral_contigs_metadata$Order[viral_contigs_metadata$Order==i] <- "Unassigned"
}
for (i in excluded_family) {
  viral_contigs_metadata$Family[viral_contigs_metadata$Family==i] <- "Unassigned"
}

viral_contigs_metadata$Family_new <- viral_contigs_metadata$Family
viral_contigs_metadata$Order_new <- viral_contigs_metadata$Order

# VC-based oredr/family assignments
for (x in unique(viral_contigs_metadata$VC[!is.na(viral_contigs_metadata$VC)])) {
  
  idx <- which(viral_contigs_metadata$VC == x)
  if (length(idx) == 1) { next }
  
  
  fam <- viral_contigs_metadata$Family_new[idx]
  fam <- unique(fam)
  fam <- fam[fam != 'Unassigned']
  if (length(fam) == 1) { viral_contigs_metadata$Family_new[idx] <- fam }
  
  
  ord <- viral_contigs_metadata$Order_new[idx]
  ord <- unique(ord)
  ord <- ord[ord != 'Unassigned']
  if (length(ord) == 1) { viral_contigs_metadata$Order_new[idx] <- ord }
  
}

viral_contigs_metadata$Family_new <- as.factor(viral_contigs_metadata$Family_new)
viral_contigs_metadata$Order_new <- as.factor(viral_contigs_metadata$Order_new)
all_viral_contigs_metadata <- viral_contigs_metadata
viral_contigs_metadata <- viral_contigs_metadata[ viral_contigs_metadata$V1 %in% row.names(RPKM_counts_VLP),]

# geNomad output merging, as geNomad might have erroneously assign phage-plasmids as plasmids
geNomad_viruses <- read.table('01.RAW_DATA/additional_qc_viruses/geNomad_viral_noneg405_99_der95_decontaminated_virus_summary.tsv', sep='\t', header=T)
geNomad_viruses <- geNomad_viruses[,c("seq_name","topology","virus_score","taxonomy")]
geNomad_viruses$type <- 'virus'
colnames(geNomad_viruses)[3] <- 'score'

genomad_plasmid <- read.table('01.RAW_DATA/additional_qc_viruses/geNomad_viral_noneg405_99_der95_decontaminated_plasmid_summary.tsv', sep='\t', header=T)
genomad_plasmid <- genomad_plasmid[,c("seq_name","topology","plasmid_score")]
genomad_plasmid$taxonomy <- "Unclassified"
genomad_plasmid$type <- 'plasmid'
colnames(genomad_plasmid)[3] <- 'score'

geNomad_merged <- rbind(geNomad_viruses, genomad_plasmid)
geNomad_merged$V1 <- ifelse(grepl("|", geNomad_merged$seq_name), 
                            sapply(strsplit(geNomad_merged$seq_name, '\\|'), "[", 1), geNomad_merged$seq_name)

colnames(geNomad_merged)[c(2:5)] <- paste0('geNomad', '_', colnames(geNomad_merged)[c(2:5)])
geNomad_merged <- geNomad_merged[geNomad_merged$V1 %in% geNomad_merged$V1,]
geNomad_merged <- geNomad_merged[!duplicated(geNomad_merged$V1),]


# checkV output: 
checkV <- read.table('01.RAW_DATA/additional_qc_viruses/CheckV_quality_summary.tsv', sep='\t', header=T)
colnames(checkV)[c(3:7,10)] <- paste0(colnames(checkV)[c(3:7,10)], '_', 'checkV')
colnames(checkV)[1] <- 'V1'
checkV <- checkV[checkV$V1 %in% viral_contigs_metadata$V1,]

# new assignment of temperate bacteriophages
phatyp <- read.csv('01.RAW_DATA/additional_qc_viruses/viral_noneg405_99_der95_decontaminated_PhaTyP_prediction.csv')
colnames(phatyp)[1] <- 'V1'
colnames(phatyp)[c(2,3)] <- paste0(colnames(phatyp)[c(2,3)], '_', 'phatyp')
phatyp <- phatyp[phatyp$V1 %in% viral_contigs_metadata$V1,]

viral_contigs_metadata <- viral_contigs_metadata %>%
  full_join(checkV[,c('V1', "provirus_checkV", "proviral_length_checkV", "gene_count_checkV", "viral_genes_checkV", "host_genes_checkV",     
                     "checkv_quality","miuvig_quality","completeness_checkV")], by='V1') %>%
  full_join(geNomad_merged[,c('V1', "geNomad_topology", "geNomad_score",    "geNomad_taxonomy", "geNomad_type")], by='V1') %>%
  full_join(phatyp[,c('V1', "Pred_phatyp",  "Score_phatyp")], by='V1')

viral_contigs_metadata[is.na(viral_contigs_metadata$geNomad_taxonomy),]$geNomad_taxonomy <- 'Unclassified'
viral_contigs_metadata$geNomad_taxonomy_new <- viral_contigs_metadata$geNomad_taxonomy
# VC-based geNomad taxonomy extension
for (x in unique(viral_contigs_metadata$VC[!is.na(viral_contigs_metadata$VC)])) {
  
  idx <- which(viral_contigs_metadata$VC == x)
  if (length(idx) == 1) { next }
  
  
  fam <- viral_contigs_metadata$geNomad_taxonomy[idx]
  fam <- unique(fam)
  fam <- fam[fam != 'Unclassified']
  if (length(fam) == 1) { viral_contigs_metadata$geNomad_taxonomy_new[idx] <- fam }
  
  
}

viral_contigs_metadata$realm_geNomad <- ifelse(grepl(";", viral_contigs_metadata$geNomad_taxonomy_new), 
                               sapply(strsplit(viral_contigs_metadata$geNomad_taxonomy_new, ';'), "[", 2), NA)
viral_contigs_metadata$kingdom <- ifelse(grepl(";", viral_contigs_metadata$geNomad_taxonomy_new), 
                                 sapply(strsplit(viral_contigs_metadata$geNomad_taxonomy_new, ';'), "[", 3), NA)
viral_contigs_metadata$phylum <- ifelse(grepl(";", viral_contigs_metadata$geNomad_taxonomy_new), 
                                sapply(strsplit(viral_contigs_metadata$geNomad_taxonomy_new, ';'), "[", 4), NA)
viral_contigs_metadata$class <- ifelse(grepl(";", viral_contigs_metadata$geNomad_taxonomy_new), 
                               sapply(strsplit(viral_contigs_metadata$geNomad_taxonomy_new, ';'), "[", 5), NA)
viral_contigs_metadata$order <- ifelse(grepl(";", viral_contigs_metadata$geNomad_taxonomy_new), 
                               sapply(strsplit(viral_contigs_metadata$geNomad_taxonomy_new, ';'), "[", 6), NA)
viral_contigs_metadata$family <- ifelse(grepl(";", viral_contigs_metadata$geNomad_taxonomy_new), 
                                sapply(strsplit(viral_contigs_metadata$geNomad_taxonomy_new, ';'), "[", 7), NA)


# metaphlan 4 output for all the NEXT samples
metaphlan4_output <- read.table('01.RAW_DATA/Metaphlan4_all_samples/LLNEXT_metaphlan_4_complete_10_02_2023.txt', sep='\t', header=T)
metaphlan4_output$NCBI_tax_id <- NULL
colnames(metaphlan4_output) <- substr(colnames(metaphlan4_output), 0, 12) 
row.names(metaphlan4_output) <- metaphlan4_output$clade_name
metaphlan4_output$clade_name <- NULL
metaphlan4_output <- metaphlan4_output[,colnames(metaphlan4_output) %in% MGS_metadata$NG_ID]
metaphlan4_output <- metaphlan4_output[rowSums(metaphlan4_output)!=0,]
metaphlan4_output <- metaphlan4_output[,MGS_metadata$NG_ID]

if (identical(MGS_metadata$NG_ID, colnames(metaphlan4_output))) {
  colnames(metaphlan4_output) <- MGS_metadata$Short_sample_ID_bact
}

# bacterial phyla table
microbiome_phyla <- metaphlan4_output[grep('p__', row.names(metaphlan4_output) ),]
microbiome_phyla <- microbiome_phyla[-grep('c__', row.names(microbiome_phyla) ), ]

# bacterial species table:
microbiome_species <- metaphlan4_output[grep('s__', row.names(metaphlan4_output) ) ,]
microbiome_species <- metaphlan4_output[-grep('t__', row.names(metaphlan4_output) ) ,]

# bacterial genera table:
microbiome_genera <- metaphlan4_output[grep('g__', row.names(metaphlan4_output)),]
microbiome_genera <- microbiome_genera[grep('s__', row.names(microbiome_genera), invert = T),]



##############################
# ANALYSIS
##############################
# Enrichment of VLP metadata with metaviromes characteristics:

## Percentage of contigs that are predicted to be viral per sample (before decontamination)
metadata_full$perc_viral_contigs <- metadata_full$N_viral_contigs/metadata_full$contigs...1000bp.*100
# Percentage of the total length of all contigs > 1 kbp occupied by virus sequences
metadata_full$proportion_viral_length <- metadata_full$Length_all_viral_contigs/metadata_full$Total_length...1000bp.*100

if ( identical( metadata_full$Short_sample_ID, colnames(RPKM_counts_VLP) ) ) {
  # viral richness
  metadata_full$viral_richness <- colSums(RPKM_counts_VLP>0)
  
  # alpha diversity
  metadata_full$viral_alpha_diversity <- diversity(as.data.frame(t(RPKM_counts_VLP)), index = "shannon")
  
  # richness of temperate viruses
  metadata_full$temperate_richness <- colSums(RPKM_counts_VLP[row.names(RPKM_counts_VLP) %in% viral_contigs_metadata[viral_contigs_metadata$temperate==1, ]$V1,] > 0)
  # proportion of temperate phages
  metadata_full$temperate_RA <- colSums(RPKM_counts_VLP[row.names(RPKM_counts_VLP) %in% viral_contigs_metadata[viral_contigs_metadata$temperate==1, ]$V1,])/colSums(RPKM_counts_VLP)*100
  
  # proportion of temperate phages based on PhaType: 
  metadata_full$phatyp_temperate_RA <- colSums(RPKM_counts_VLP[row.names(RPKM_counts_VLP) %in% viral_contigs_metadata[!is.na(viral_contigs_metadata$Pred_phatyp) & viral_contigs_metadata$Pred_phatyp=='temperate' & !(viral_contigs_metadata$checkv_quality %in% c("Not-determined", "Low-quality")),]$V1,])/colSums(RPKM_counts_VLP)*100
  metadata_full$phatyp_temperate_ricnhess <- colSums(RPKM_counts_VLP[row.names(RPKM_counts_VLP) %in%  viral_contigs_metadata[!is.na(viral_contigs_metadata$Pred_phatyp) & viral_contigs_metadata$Pred_phatyp=='temperate' & !(viral_contigs_metadata$checkv_quality %in% c("Not-determined", "Low-quality")),]$V1,] > 0)
  # proportion of temperate phages based on PhaType if no cut-offs are applied for the quality of the assembled genome
  metadata_full$phatyp_temperate_RA_all <- colSums(RPKM_counts_VLP[row.names(RPKM_counts_VLP) %in% viral_contigs_metadata[!is.na(viral_contigs_metadata$Pred_phatyp) & viral_contigs_metadata$Pred_phatyp=='temperate',]$V1,])/colSums(RPKM_counts_VLP)*100
  metadata_full$phatyp_temperate_ricnhess_all <- colSums(RPKM_counts_VLP[row.names(RPKM_counts_VLP) %in%  viral_contigs_metadata[!is.na(viral_contigs_metadata$Pred_phatyp) & viral_contigs_metadata$Pred_phatyp=='temperate',]$V1,] > 0)
  
  # proportion of temperate phages based on PhaType if no cut-offs are applied for the quality of the assembled genome
  metadata_full$phatyp_temperate_RA_most <- colSums(RPKM_counts_VLP[row.names(RPKM_counts_VLP) %in% viral_contigs_metadata[!is.na(viral_contigs_metadata$Pred_phatyp) & viral_contigs_metadata$Pred_phatyp=='temperate' & viral_contigs_metadata$checkv_quality!='Not-determined',]$V1,])/colSums(RPKM_counts_VLP)*100
  metadata_full$phatyp_temperate_ricnhess_most <- colSums(RPKM_counts_VLP[row.names(RPKM_counts_VLP) %in%  viral_contigs_metadata[!is.na(viral_contigs_metadata$Pred_phatyp) & viral_contigs_metadata$Pred_phatyp=='temperate' & viral_contigs_metadata$checkv_quality!='Not-determined',]$V1,] > 0)
  
  # proportion of temperate phages based on PhaType if no cut-offs are applied for the quality of the assembled genome but geNomad recongnize it as a virus
  metadata_full$phatyp_temperate_RA_genomad <- colSums(RPKM_counts_VLP[row.names(RPKM_counts_VLP) %in% viral_contigs_metadata[!is.na(viral_contigs_metadata$Pred_phatyp) & viral_contigs_metadata$Pred_phatyp=='temperate' & !is.na(viral_contigs_metadata$geNomad_type) & viral_contigs_metadata$geNomad_type=='virus',]$V1,])/colSums(RPKM_counts_VLP)*100
  metadata_full$phatyp_temperate_ricnhess_genomad <- colSums(RPKM_counts_VLP[row.names(RPKM_counts_VLP) %in%  viral_contigs_metadata[!is.na(viral_contigs_metadata$Pred_phatyp) & viral_contigs_metadata$Pred_phatyp=='temperate' & !is.na(viral_contigs_metadata$geNomad_type) & viral_contigs_metadata$geNomad_type=='virus',]$V1,] > 0)
}


# Enrichment of MGS metadata with viromics characteristics from MGS alignment:

if ( identical( MGS_metadata$Short_sample_ID, colnames(RPKM_counts_MGS) ) ) {
  # viral richness
  MGS_metadata$viral_richness_MGS <- colSums(RPKM_counts_MGS>0)
  
  # alpha diversity
  MGS_metadata$viral_alpha_diversity_MGS <- diversity(as.data.frame(t(RPKM_counts_MGS)), index = "shannon")
  
  # richness of temperate viruses
  MGS_metadata$temperate_richness_MGS <- colSums(RPKM_counts_MGS[row.names(RPKM_counts_MGS) %in% viral_contigs_metadata[viral_contigs_metadata$temperate==1, ]$V1,] > 0)
  # proportion of temperate phages
  MGS_metadata$temperate_RA_MGS <- colSums(RPKM_counts_MGS[row.names(RPKM_counts_MGS) %in% viral_contigs_metadata[viral_contigs_metadata$temperate==1, ]$V1,])/colSums(RPKM_counts_MGS)*100
}


# bacterial diversity:

if ( identical( colnames(microbiome_species), MGS_metadata$Short_sample_ID_bact ) ) {
  
  MGS_metadata$bacterial_alpha_diversity <- diversity(as.data.frame(t(microbiome_species)), index = "shannon")

  }

# adding the proportion of unknown in metaphlan4

if ( identical( colnames(metaphlan4_output), MGS_metadata$Short_sample_ID_bact ) ) {
  
  MGS_metadata$metaphlan_unknown_perc <- t(metaphlan4_output)[,1]
  
}


# taxonomy abundance table

# FAMILY RANK
RPKM_counts_VLP_by_family <- data.frame(matrix(ncol = ncol(RPKM_counts_VLP), nrow=length(levels(viral_contigs_metadata$Family_new) ) ) )
colnames(RPKM_counts_VLP_by_family) <- colnames(RPKM_counts_VLP)
row.names(RPKM_counts_VLP_by_family) <- levels(viral_contigs_metadata$Family_new)

for (i in colnames(RPKM_counts_VLP_by_family) ) {
  
  for (j in row.names(RPKM_counts_VLP_by_family) ) {
    
    RPKM_counts_VLP_by_family[row.names(RPKM_counts_VLP_by_family)==j,i] <- sum( RPKM_counts_VLP[ c( viral_contigs_metadata[ viral_contigs_metadata$Family_new==j, ]$V1 ) , i] )
    
  }
  
}

#ORDER RANK
RPKM_counts_VLP_by_order <- data.frame(matrix(ncol = ncol(RPKM_counts_VLP), nrow=length(levels(viral_contigs_metadata$Order_new))))
colnames(RPKM_counts_VLP_by_order) <- colnames(RPKM_counts_VLP)
row.names(RPKM_counts_VLP_by_order) <- levels(viral_contigs_metadata$Order_new)

for (i in colnames(RPKM_counts_VLP_by_order) ) {
  
  for (j in row.names(RPKM_counts_VLP_by_order) ) {
    
    RPKM_counts_VLP_by_order[row.names(RPKM_counts_VLP_by_order)==j,i] <- sum( RPKM_counts_VLP[ c( viral_contigs_metadata[ viral_contigs_metadata$Order_new==j, ]$V1 ) , i] )
    
  }
  
}

# enrichment of viral metadata with bacterial characteristics:
metadata_full <- merge(metadata_full, MGS_metadata[,c("Universal_fecal_ID","bacterial_alpha_diversity")], all.x = T, by="Universal_fecal_ID")

# changing the columnnames for uniformity:
metadata_full$PseudoID <- NULL
metadata_full$Lifelines_ID <- NULL
metadata_full$NG_ID <- NULL
colnames(metadata_full)[c(5,6)] <- c("Timepoint","Type")
metadata_full[metadata_full$Type=="Baby",]$Type <- "Infant"
metadata_full$reads_lost_QC <- (metadata_full$Raw_reads - metadata_full$Clean_reads)/metadata_full$Raw_reads
metadata_full$BOX <- NULL
metadata_full$POSITION <- NULL
metadata_full$SEQUENCING <- NULL
metadata_full$ID_1 <- NULL
metadata_full$ID_02 <- NULL
metadata_full$ID_3 <- NULL
metadata_full$ID_4 <- NULL
metadata_full$ID.5 <- NULL
metadata_full$Short_sample_ID_bact <- NA

#rearranging the order of the vlp_metadata:

metadata_full <- metadata_full[,c(1:4, 11, 6,5,26,40:42,64,8,12:25, 27:32,46,65,47, 48,45,
                                  33:39,43:44,49:63,7,9,10)]
colnames(metadata_full)[grep('LABEL',colnames(metadata_full))] <- 'Old_ID'
colnames(metadata_full)[grep('Unified_ID',colnames(metadata_full))] <- 'NG_ID'

MGS_metadata$Raw_reads <- MGS_metadata$raw_reads_FQ_1 + MGS_metadata$raw_reads_FQ_2
MGS_metadata$Clean_reads <- MGS_metadata$clean_reads_FQ_1 + MGS_metadata$clean_reads_FQ_2
MGS_metadata$Human_reads <- MGS_metadata$human_reads_FQ_1 + MGS_metadata$human_reads_FQ_2
MGS_metadata <- MGS_metadata[,!( colnames(MGS_metadata) %in% c("raw_reads_FQ_1", "raw_reads_FQ_2", 
                                                               "clean_reads_FQ_1", "clean_reads_FQ_2",
                                                               "human_reads_FQ_1", "human_reads_FQ_2") ) ]
MGS_metadata$pilot_old_ID <- NULL
MGS_metadata$FAM_ID <- paste0('FAM', str_pad(MGS_metadata$FAM_ID, 4, pad = "0"))
MGS_metadata$duplicated <- NULL
MGS_metadata[MGS_metadata$Timepoint=='MB',]$Timepoint <- 'B'

keep.mgs.metadata.intact <- MGS_metadata
#rearranging the order of the mgs_metadata:
MGS_metadata <- MGS_metadata[,c(33,1,4,2,3,5:7,43:45,9,10, 11:32, 34:42, 8)]



##############################
# OUTPUT
##############################
write.table(RPKM_counts_VLP, "02.CLEAN_DATA/RPKM_counts_VLP.txt", sep = '\t', quote=F)
write.table(metadata_full, "02.CLEAN_DATA/VLP_metadata_with_phenos.txt", sep='\t', quote=F, row.names = F)

write.table(RPKM_counts_MGS, "02.CLEAN_DATA/RPKM_counts_MGS.txt", sep='\t', quote=F)
write.table(MGS_metadata, "02.CLEAN_DATA/MGS_metadata_with_phenos.txt", sep='\t', quote=F, row.names=F)

write.table(all_viral_contigs_metadata, "02.CLEAN_DATA/all_viral_contigs_metadata.txt", sep='\t', quote = F, row.names = F)
write.table(viral_contigs_metadata, "02.CLEAN_DATA/VLP_viral_contigs_metadata.txt", sep='\t', quote = F, row.names = F)

write.table(RPKM_counts_VLP_by_family, "02.CLEAN_DATA/RPKM_counts_VLP_by_tax_family.txt", sep = '\t', quote=F)
write.table(RPKM_counts_VLP_by_order, "02.CLEAN_DATA/RPKM_counts_VLP_by_tax_order.txt", sep = '\t', quote=F)

write.table(microbiome_species, "02.CLEAN_DATA/Microbiome_species_unfiltred.txt", sep='\t', quote=F)
write.table(microbiome_genera, "02.CLEAN_DATA/Microbiome_genera_unfiltred.txt", sep='\t', quote=F)
write.table(microbiome_phyla, "02.CLEAN_DATA/Microbiome_phyla_unfiltred.txt", sep='\t', quote=F)
