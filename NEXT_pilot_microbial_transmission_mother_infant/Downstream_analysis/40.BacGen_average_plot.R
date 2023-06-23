setwd('~/Desktop/Projects_2022/NEXT_pilot_FUP/')

#############################################################
# Here we explore the bacterial genera-level composition
# of infant and maternal gut vs vOTUs aggregated at host
# genera level
#############################################################

##############################
# Functions
##############################
get_average_per_timepoint <- function(counts_df, metadata, ID, type, timepoints){
  
  df_average <- as.data.frame(matrix(NA, 
                                     nrow = nrow(counts_df), 
                                     ncol = length(timepoints), 
                                     dimnames = list(row.names(counts_df), timepoints)))
  
  for (TP in timepoints) {
    
    if ( identical( names(rowMeans(counts_df[,metadata[metadata$Type==type & metadata$Timepoint==TP, ID]])), row.names(df_average) ) ) {
      
      df_average[,TP] <- rowMeans(counts_df[,metadata[metadata$Type==type & metadata$Timepoint==TP, ID]])
      
    }
    
  }
  
  df_average$all <- rowMeans(counts_df[,metadata[metadata$Type==type,ID]])
  df_average <- df_average[order(df_average$all, decreasing = T),]
  
  return(df_average)
  
}

wrapper_melt <- function(df_average, taxa, timepoints, type){
  
  df_melt <- df_average[taxa, timepoints]
  df_melt[,"taxa"] <- row.names(df_melt)
  df_melt <- melt(df_melt)
  df_melt$Type <- type
  
  return(df_melt)
}
##############################
# Loading libraries
##############################

##############################
# Input data
##############################

Renoir <- met.brewer("Renoir")
Kandinsky <- met.brewer("Kandinsky")
Demuth <- met.brewer("Demuth")

MGS_metadata <- read.table('02.CLEAN_DATA/MGS_metadata_final_10_05_2023.txt', sep='\t', header=T)
Infant_timepoints <- c('M1', 'M2', 'M3', 'M6', 'M9', 'M12')
Mother_timepoints <- c("P3", "P7", "B", "M1", "M2", "M3")
VLP_metadata <- read.table('02.CLEAN_DATA/VLP_metadata_final_10_05_2023.txt', sep='\t', header=T)


bac_gen <- read.table('02.CLEAN_DATA/Microbiome_genera_unfiltred.txt', sep='\t', header=T)


BacGen_Mother_average <- get_average_per_timepoint(bac_gen, MGS_metadata, "Short_sample_ID_bact", "Mother", Mother_timepoints)
BacGen_Infant_average <- get_average_per_timepoint(bac_gen, MGS_metadata, "Short_sample_ID_bact", "Infant", Infant_timepoints)

# choosing top 10 from mothers and infants:
BacGen <- unique(c(row.names(BacGen_Mother_average[1:10,]), row.names(BacGen_Infant_average[1:10,])))

BacGen_Mother <- wrapper_melt(BacGen_Mother_average, BacGen, Mother_timepoints, 'Mother')
BacGen_Infant <- wrapper_melt(BacGen_Infant_average, BacGen, Infant_timepoints, 'Infant')

BacGen_average <- rbind(BacGen_Infant, BacGen_Mother)
BacGen_average$variable <- factor(BacGen_average$variable, levels=c("P3", "P7", "B",
                                                                    "M1", "M2", "M3",
                                                                    "M6", "M9", "M12"), ordered = T)
BacGen_average$taxa <- gsub('.*g__', '', BacGen_average$taxa)


# vOTUs aggreagtion:
host_assignment <- read.table('02.CLEAN_DATA/Host_prediction_to_genus_m90_refined_taxonomy_no_generalists.txt', sep='\t', header=T)

RPKM_counts_VLP <- read.table('02.CLEAN_DATA/RPKM_counts_VLP.txt', sep='\t', header=T)
RPKM_counts_VLP$Virus <- row.names(RPKM_counts_VLP)

host_assignment_combined <- merge(host_assignment[,c("Virus", "Genus")], RPKM_counts_VLP, by='Virus', all.y = T)
host_assignment_combined[is.na(host_assignment_combined$Genus),"Genus"] <- 'Unassigned'
host_assignment_combined$Virus <- NULL

RPKM_counts_BacGen <- aggregate(.~Genus, host_assignment_combined, sum)
row.names(RPKM_counts_BacGen) <- RPKM_counts_BacGen$Genus
RPKM_counts_BacGen$Genus <- NULL
RPKM_counts_BacGen_RA <- as.data.frame( t(t( RPKM_counts_BacGen )/colSums(RPKM_counts_BacGen))*100  )

RPKM_counts_MGS <- read.table('02.CLEAN_DATA/RPKM_counts_MGS.txt', sep='\t', header=T)
RPKM_counts_MGS$Virus <- row.names(RPKM_counts_MGS)
host_assignment_combined_MGS <- merge(host_assignment[,c("Virus", "Genus")], RPKM_counts_MGS, by='Virus', all.y = T)
host_assignment_combined_MGS[is.na(host_assignment_combined_MGS$Genus),"Genus"] <- 'Unassigned'
host_assignment_combined_MGS$Virus <- NULL

RPKM_counts_MGS_BacGen <- aggregate(.~Genus, host_assignment_combined_MGS, sum)
row.names(RPKM_counts_MGS_BacGen) <- RPKM_counts_MGS_BacGen$Genus
RPKM_counts_MGS_BacGen$Genus <- NULL
RPKM_counts_MGS_BacGen_RA <- as.data.frame( t(t( RPKM_counts_MGS_BacGen )/colSums(RPKM_counts_MGS_BacGen))*100  )
# getting absent in VLP timepoints:
RPKM_counts_MGS_BacGen_RA <- RPKM_counts_MGS_BacGen_RA[,MGS_metadata[ MGS_metadata$Timepoint %in% c('P3', 'M9'),  ]$Short_sample_ID]

RPKM_counts_BacGen_all <- merge(RPKM_counts_BacGen_RA, RPKM_counts_MGS_BacGen_RA, by='row.names', all=T)
row.names(RPKM_counts_BacGen_all) <- RPKM_counts_BacGen_all$Row.names
RPKM_counts_BacGen_all$Row.names <- NULL
RPKM_counts_BacGen_all[is.na(RPKM_counts_BacGen_all)] <- 0


upd_metadata <- rbind(VLP_metadata[,intersect(colnames(MGS_metadata), colnames(VLP_metadata))],
                  MGS_metadata[MGS_metadata$Timepoint %in% c("P3", "M9"),intersect(colnames(MGS_metadata), colnames(VLP_metadata))])

HostGen_Mother_average <- get_average_per_timepoint(RPKM_counts_BacGen_all, upd_metadata, "Short_sample_ID", "Mother", Mother_timepoints)
HostGen_Infant_average <- get_average_per_timepoint(RPKM_counts_BacGen_all, upd_metadata, "Short_sample_ID", "Infant", Infant_timepoints)

# choosing top 10 from mothers and infants:
HostGen <- unique(c(row.names(HostGen_Mother_average[1:10,]), row.names(HostGen_Infant_average[1:10,])))
HostGen <- HostGen[-1]
HostGen <- HostGen[ -grep("TM7x", HostGen) ] # unlikely to be found

HostGen_Mother <- wrapper_melt(HostGen_Mother_average, HostGen, Mother_timepoints, 'Mother')
HostGen_Infant <- wrapper_melt(HostGen_Infant_average, HostGen, Infant_timepoints, 'Infant')

HostGen_average <- rbind(HostGen_Infant, HostGen_Mother)
HostGen_average$variable <- factor(HostGen_average$variable, levels=c("P3", "P7", "B",
                                                                    "M1", "M2", "M3",
                                                                    "M6", "M9", "M12"), ordered = T)
HostGen_average$Sequencing <- ifelse( HostGen_average$variable %in% c('P3', 'M9'), 'MGS', 'VLP' )



# colors: 
Gen_palette <- data.frame(sort(unique(c(HostGen, gsub('.*g__','',BacGen)))), c(Kandinsky[1], 
                                                                                Renoir[1], 
                                                                                Kandinsky[2], 
                                                                                Renoir[2:7],
                                                                                Renoir[8:9],
                                                                               Renoir[11], Renoir[12],
                                                                                Kandinsky[3:4],
                                                                                Demuth[1]))
colnames(Gen_palette) <- c("Genus", "Color")

##############################
# ANALYSIS
##############################

ggplot(BacGen_average, aes(variable, value, fill=taxa)) + 
  geom_bar(position = "fill", stat = "identity") + 
  facet_wrap(~Type, scales = 'free_x', nrow = 1) + 
  xlab(label = "Timepoint") + 
  ylab(label = "Rescaled relative abundance") + 
  labs(fill='Bacterial genera') +
  theme_bw() + 
  scale_fill_manual(values=Gen_palette[Gen_palette$Genus %in% unique(BacGen_average$taxa), ]$Color) +
  theme(axis.text.y=element_text(size=16), 
        #axis.text.x = element_blank(),
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 12))

ggplot(HostGen_average, aes(variable, value, fill=taxa, alpha=Sequencing)) + 
  geom_bar(position = "fill", stat = "identity") + 
  facet_wrap(~Type, scales = 'free_x', nrow = 1) + 
  xlab(label = "Timepoint") + 
  ylab(label = "Rescaled relative abundance") + 
  labs(fill='Host genera') +
  theme_bw() + 
  scale_fill_manual(values=Gen_palette[Gen_palette$Genus %in% unique(HostGen_average$taxa), ]$Color) +
  scale_alpha_manual(values = c(0.8, 1)) +
  theme(axis.text.y=element_text(size=16), 
        #axis.text.x = element_blank(),
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 12))

##############################
# OUTPUT
##############################
write.table(BacGen_average, '02.CLEAN_DATA/PREPARED_DATA_FOR_PLOTS/Fig2H_BacGen_average_abundance_over_time.txt', sep='\t', quote=F, row.names=F)
write.table(HostGen_average, '02.CLEAN_DATA/PREPARED_DATA_FOR_PLOTS/Fig2G_HostGen_vOTUs_average_abundance_over_time.txt', sep = '\t', quote=F, row.names=F)
write.table(Gen_palette, '02.CLEAN_DATA/PREPARED_DATA_FOR_PLOTS/Fig2GH_gen_palette_colors.txt', sep='\t', row.names=F)
