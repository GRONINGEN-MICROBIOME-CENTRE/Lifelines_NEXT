setwd('~/Desktop/Projects_2022/NEXT_pilot_FUP/')

#############################################################
# Here we identify viruses that are shared between maternal 
# and infant samples from the same families for the further
# analysis of their genetic identity.
#############################################################

##############################
# Functions
##############################
# MODIFIED RPKM_table FUNCTION WITH A 0.95 CUT-OFF FOR PRESENCE
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
    counts[coverage<0.95] <-0 #ATTENTION, IN THIS FUNCTION THE CUT-OFF HAS BEEN CHANGED
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
library(reshape2)
library(ggplot2)
library(dplyr)
##############################
# Input data & data transformation
##############################
metadata_VLP <- read.table("02.CLEAN_DATA/VLP_metadata_with_phenos.txt", sep='\t', header=T)
metadata_VLP$Timepoint <- factor(metadata_VLP$Timepoint, levels=c('P7', 'B', 'M1', 'M2', 'M3', 'M6', 'M12'), ordered=T)

metadata_MGS <- read.table("02.CLEAN_DATA/MGS_metadata_with_phenos.txt", sep='\t', header=T)
metadata_MGS$Timepoint <- factor(metadata_MGS$Timepoint, levels=c('P3','P7', 'B', 'M1', 'M2', 'M3', 'M6', 'M9','M12'), ordered=T)

# combined metadata for VLP and MGS samples:
metadata <- rbind(metadata_MGS[,c("Short_sample_ID", "NG_ID", "FAM_ID", "Timepoint", "Type", "Individual_ID")],
                  metadata_VLP[,c("Short_sample_ID", "NG_ID", "FAM_ID", "Timepoint", "Type", "Individual_ID")])
metadata$source <- 'MGS'
metadata[grep('V', metadata$Short_sample_ID),]$source <- 'VLP'
metadata$Type <- factor(metadata$Type, levels=c('Mother', 'Infant'), ordered=T)

# table with >0.75 length coverage cut-off for presence 
RPKM_counts_VLP <- read.table("02.CLEAN_DATA/RPKM_counts_VLP.txt", sep='\t', header=T)
# table with >0.75 length coverage cut-off for presence 
RPKM_counts_MGS <- read.table("02.CLEAN_DATA/RPKM_counts_MGS.txt", sep='\t', header=T)

# table with >=0.95 length coverage cut-off for presence
VLP_vircounts <- read.table("01.RAW_DATA/Counts_tables_final_noneg405_99_der95/bowtie2_read_counts_with_full.txt", sep='\t', header = T, row.names = 1)
VLP_contig_coverage <- read.table('01.RAW_DATA/Counts_tables_final_noneg405_99_der95/coverage_table.txt', sep='\t', header = T, row.names = 1)
RPKM_counts_VLP_0.95 <- RPKM_table(VLP_vircounts, VLP_contig_coverage)
RPKM_counts_VLP_0.95 <- RPKM_counts_VLP_0.95[,colnames(RPKM_counts_VLP_0.95) %in% metadata_VLP$NG_ID]
RPKM_counts_VLP_0.95 <- RPKM_counts_VLP_0.95[rowSums(RPKM_counts_VLP_0.95)>0,] #1,445 contigs did not survive 0.95 presence cut-off
RPKM_counts_VLP_0.95 <- RPKM_counts_VLP_0.95[,metadata_VLP$NG_ID]
if (identical(metadata_VLP$NG_ID, colnames(RPKM_counts_VLP_0.95))) {
  colnames(RPKM_counts_VLP_0.95) <- metadata_VLP$Short_sample_ID
}

# table with >=0.95 length coverage cut-off for presence
MGS_vircounts <- read.table("01.RAW_DATA/Counts_tables_MGS/bowtie2_read_counts_with_full.txt", sep='\t', header = T, row.names = 1)
MGS_contig_coverage <- read.table('01.RAW_DATA/Counts_tables_MGS/coverage_table.txt', sep='\t', header = T, row.names = 1)
RPKM_counts_MGS_0.95 <- RPKM_table(MGS_vircounts, MGS_contig_coverage)
RPKM_counts_MGS_0.95 <- RPKM_counts_MGS_0.95[,colnames(RPKM_counts_MGS_0.95) %in% metadata_MGS$NG_ID]
RPKM_counts_MGS_0.95 <- RPKM_counts_MGS_0.95[rowSums(RPKM_counts_MGS_0.95)>0,] # 36,476 contigs are not detected in total MGS samples
RPKM_counts_MGS_0.95 <- RPKM_counts_MGS_0.95[,metadata_MGS$NG_ID]
if (identical(metadata_MGS$NG_ID, colnames(RPKM_counts_MGS_0.95))) {
  colnames(RPKM_counts_MGS_0.95) <- metadata_MGS$Short_sample_ID
}
RPKM_counts_MGS_0.95 <- RPKM_counts_MGS_0.95[row.names(RPKM_counts_MGS_0.95) %in% row.names(RPKM_counts_VLP_0.95), ]

contigs_metadata <- read.table('./02.CLEAN_DATA/VLP_viral_contigs_metadata.txt', sep='\t', header=T)

# combined RPKM table for VLP and MGS samples
RPKM_combined <- merge(RPKM_counts_VLP_0.95, RPKM_counts_MGS_0.95, by='row.names', all = T)
row.names(RPKM_combined) <- RPKM_combined$Row.names
RPKM_combined$Row.names <- NULL
RPKM_combined[is.na(RPKM_combined)] <- 0
##############################
# ANALYSIS
##############################

# viruses that are present in maternal VLP or MGS and also present in infant VLP:
phages_of_interest <- list()

RPKM_combined$MOCK <- 0 #a workaround for families where only 1 sample is available for the timepoint

for (i in unique(metadata_VLP$FAM_ID)) {
    
    phages_of_interest[[i]] <- c(row.names(RPKM_combined[  rowSums( RPKM_combined[ , c( metadata[ metadata$FAM_ID==i & metadata$Type=='Mother', ]$Short_sample_ID, 'MOCK' ) ] ) != 0  &
                                                            rowSums( RPKM_combined[ , c( metadata[ metadata$FAM_ID==i & metadata$source=='VLP' & metadata$Type=='Infant', ]$Short_sample_ID, 'MOCK', 'MOCK') ]  ) != 0, ]
    ))
  
}

# collection of information about viruses that fullfil the presence/absence criteria within families
phages_of_interest_metadata <- contigs_metadata[ contigs_metadata$V1 %in% unique( unname( unlist(phages_of_interest) )), ]
phages_of_interest_metadata$prevalence_VLP <- NA
phages_of_interest_metadata$prevalence_MGS <- NA
phages_of_interest_metadata$nz_mean_ab_VLP <- NA
phages_of_interest_metadata$nz_mean_ab_MGS <- NA
phages_of_interest_metadata$N_positive_families <- NA

for (i in phages_of_interest_metadata$V1) {
  phages_of_interest_metadata[phages_of_interest_metadata$V1==i,]$prevalence_VLP <- sum(RPKM_counts_VLP_0.95[i,]!=0)
  
  phages_of_interest_metadata[phages_of_interest_metadata$V1==i,]$prevalence_MGS <- sum(RPKM_counts_MGS_0.95[i,]!=0)
  
  phages_of_interest_metadata[phages_of_interest_metadata$V1==i,]$nz_mean_ab_VLP <- rowMeans(replace(RPKM_counts_VLP_0.95[i,], RPKM_counts_VLP_0.95[i,] == 0, NA), na.rm = TRUE)
  
  phages_of_interest_metadata[phages_of_interest_metadata$V1==i,]$nz_mean_ab_MGS <- rowMeans(replace(RPKM_counts_MGS_0.95[i,], RPKM_counts_MGS_0.95[i,] == 0, NA), na.rm = TRUE)
  
  phages_of_interest_metadata[phages_of_interest_metadata$V1==i,]$N_positive_families <-  length( unique( metadata[ metadata$Short_sample_ID %in% colnames(RPKM_combined[i, colSums(RPKM_combined[i,])!=0]), ]$FAM_ID ) )
  
}



phages_of_interest_for_check <- phages_of_interest_metadata[ phages_of_interest_metadata$length >= 3000 &
                                                             phages_of_interest_metadata$N_positive_families >=5 & 
                                                            ( phages_of_interest_metadata$Topology==1 | phages_of_interest_metadata$miuvig_quality=='High-quality'),]

phages_of_interest_for_check$ContigID_easy <- paste0('L',phages_of_interest_for_check$length, '_', 'LS', phages_of_interest_for_check$temperate)

RPKM_combined <- RPKM_combined[row.names(RPKM_combined) %in% phages_of_interest_for_check$V1,]


# plotting the abundance graphs for phage candidates per family:
phages_origin_to_plot <- lapply(phages_of_interest, function(x) intersect(x, phages_of_interest_for_check$V1))
phages_origin_to_plot <- phages_origin_to_plot[sapply(phages_origin_to_plot, length)>0]


plot_list = list()

for (fam_id in names(phages_origin_to_plot)) {
  
  for (seq_id in unlist(unname(phages_origin_to_plot[fam_id])) ) {
    
    PLOT_CONTIG <- melt(RPKM_combined[seq_id, metadata[metadata$FAM_ID==fam_id,]$Short_sample_ID  ])
    colnames(PLOT_CONTIG)[1] <- 'Short_sample_ID'
    PLOT_CONTIG <- merge(PLOT_CONTIG, metadata, by='Short_sample_ID')
    PLOT_CONTIG <- PLOT_CONTIG %>% arrange(Timepoint)
    
    p <- ggplot(PLOT_CONTIG, aes(Timepoint, value, color=Individual_ID, group=Individual_ID) ) + 
      geom_point() + 
      geom_path() +
      facet_wrap(~source+Type) +
      theme_bw()+
      scale_y_log10() + 
      labs (y="log-transformed abundance", x="Timepoint") +
      ggtitle(seq_id) + 
      theme(axis.text=element_text(size=12), 
            axis.title=element_text(size=16,face="bold"),
            strip.text.x = element_text(size = 12),
            legend.text = element_text(size=10),
            legend.title = element_text(size=12, face="bold")) +
      labs(colour = "Type") 
    plot_list[[paste0(fam_id, seq_id)]] = p
  }
}

pdf("plots_phages_of_interest_check_0.95_cut-off.pdf")
for (i in names(plot_list)) {
  print(plot_list[[i]])
}
dev.off()

RPKM_combined <- RPKM_combined[,colSums(RPKM_combined)!=0]

###### OUTPUT #####
write.table(phages_of_interest_for_check, '03.SCRIPTS/NEXT_pilot_FUP_bf_origin/Data_for_Alex/Viruses_shared_min_5_families_UPD_final.txt', sep='\t', quote=F, row.names = F)
write.table(RPKM_combined, '03.SCRIPTS/NEXT_pilot_FUP_bf_origin/Data_for_Alex/RPKM_counts_combined_0.95_UPD_final.txt', sep='\t', quote=F)



