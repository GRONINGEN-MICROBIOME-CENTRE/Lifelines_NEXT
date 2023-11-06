setwd('~/Desktop/Projects_2022/NEXT_pilot_FUP/')

#############################################################
# 
#
#############################################################

##############################
# Functions
##############################

##############################
# Loading libraries
##############################

library(tidyverse)
library(gggenomes)
library(readr)
library(pafr)
library(reshape2)
library(MetBrewer)

met.brewer('Austria')
##############################
# Input data
##############################
VLP_metadata <- read.table('02.CLEAN_DATA/VLP_metadata_final_10_05_2023.txt', sep='\t', header=T)
VLP_metadata$Alex_ID <- paste0(VLP_metadata$FAM_ID, '_', VLP_metadata$Type, '_', substr(VLP_metadata$Short_sample_ID, 1,1), '_VLP_', VLP_metadata$Timepoint)
VLP_metadata$Color <- 'grey'
VLP_metadata[VLP_metadata$FAM_ID=="FAM0234",]$Color <- 'red'
MGS_metadata <- read.table('02.CLEAN_DATA/MGS_metadata_final_10_05_2023.txt', sep='\t', header=T)
MGS_metadata$Alex_ID <- paste0(MGS_metadata$FAM_ID, '_', MGS_metadata$Type, '_', substr(MGS_metadata$Short_sample_ID_bact, 1,1), '_MGS_', MGS_metadata$Timepoint)
MGS_metadata$Color <- 'grey'
MGS_metadata[MGS_metadata$FAM_ID=="FAM0234",]$Color <- 'red'

microbiome <- read.table('02.CLEAN_DATA/Microbiome_species_unfiltred.txt', sep='\t', header = T)
B_bifidum_positive <- colnames(microbiome[,microbiome[grep('Bifidobacterium_bifidum', row.names(microbiome)),] !=0])
MGS_metadata$B_bifidum_positive <- ifelse(MGS_metadata$Short_sample_ID_bact %in% B_bifidum_positive, "YES", "NO")
MGS_metadata <- MGS_metadata[order(MGS_metadata$Color),]
MGS_metadata[MGS_metadata$FAM_ID=="FAM0234" & MGS_metadata$Timepoint %in% c('P3', 'P7', 'B'),]$Color <- 'blue'
MGS_metadata[MGS_metadata$FAM_ID=="FAM0234" & MGS_metadata$Type=="Mother" & MGS_metadata$Timepoint %in% c('M1', 'M2', 'M3'),]$Color <- 'green'


RPKM_counts_VLP <- read.table('02.CLEAN_DATA/RPKM_counts_VLP.txt', sep='\t', header=T)
L34922_LS1_positive <- colnames(RPKM_counts_VLP[,RPKM_counts_VLP[grep('LN_4F02_VL_264_NODE_17_length_34922_cov_30951.540224', row.names(RPKM_counts_VLP)),] !=0])
VLP_metadata$L34922_LS1_positive <- ifelse(VLP_metadata$Short_sample_ID %in% L34922_LS1_positive, "YES", "NO")
VLP_metadata <- VLP_metadata[order(VLP_metadata$Color),]
VLP_metadata[VLP_metadata$FAM_ID=="FAM0234" & VLP_metadata$Timepoint %in% c('P3', 'P7', 'B'),]$Color <- 'blue'
VLP_metadata[VLP_metadata$FAM_ID=="FAM0234" & VLP_metadata$Type=="Mother" & VLP_metadata$Timepoint %in% c('M1', 'M2', 'M3'),]$Color <- 'green'

##############################
# ANALYSIS
##############################

ali3 <- read_paf('03.SCRIPTS/NEXT_pilot_FUP_bf_origin/Data_for_Alex/Phages_vs_bacteria_L37775_LS1.paf')
ali3 <- ali3[ali3$qname=='C03X03483E05_NODE_35_length_84704_cov_141.800057',]
ali3 <- ali3[ali3$tname=='LN_4F04_VL_266_NODE_141_length_37775_cov_36892.259438',]
ali3$qname <- "B. scardovi patched genome fragment"
ali3$tname <- "L37775_LS1"

plot_synteny(ali3, 
             q_chrom="B. scardovi patched genome fragment", 
             t_chrom="L37775_LS1", 
             centre=TRUE) + 
  theme_bw() + 
  theme_linedraw() +
  theme(panel.grid = element_blank(), 
        axis.title = element_text(face="bold"))


### B. bifidum synteny: 
ali <- read_paf('03.SCRIPTS/NEXT_pilot_FUP_bf_origin/scaffolding/alignment_to_scaffolded2.paf')
ali$qname <- 'L34922_LS1'
ali[ali$tname=='MGYG000132487_6_RagTag',]$tname <- 'B. bifidum patched genome fragment'
plot_synteny(ali, 
             q_chrom="L34922_LS1", 
             t_chrom="B. bifidum patched genome fragment", 
             centre=TRUE) + 
  theme_linedraw() +
  theme(panel.grid = element_blank(), 
        axis.title = element_text(face="bold")) +
  annotate(geom="text", x = 750000, y=1.5, label="pident > 99%\nqcov=100%\ne-value=0")



for (seq_id in c('MGYG000132487_6_RagTag')) {
  
  len <- 990913
  
  plot(NA, 
       xlim=c(0, 150000),
       ylim=c(1, 10^4),
       yaxs='i',
       log='y',
       xlab=paste0(seq_id, ', nt'),
       ylab='mean depth + 1')
  
  w_center <- seq(from=501, to=len-500, by=1001)
  w_from <- w_center-500
  w_to <- w_center+500
  w_to[ length(w_to)] <- len #expansion of the last fragment in case it is smaller than window size
  
  # df to keep mean coverage over window
  B_bifidum_coverage_MGS <- data.frame(w_center)

  for (sample_id in MGS_metadata[MGS_metadata$B_bifidum_positive=="YES",]$NG_ID) {
    
    tab <- read.table(paste0('~/Desktop/Projects_2022/NEXT_pilot_FUP/03.SCRIPTS/NEXT_pilot_FUP_bf_origin/L34922_LS1/MGS_alignments/', sample_id, '.target.cov.txt'), sep='\t', header=F, fill=T)
    
    COVERAGE <- rep(0, len)
    COVERAGE[tab$V2] <- tab$V4
    
    V <- c()
    V <- sapply(seq_along(w_center), function(i) {
      mean( COVERAGE[ w_from[i]:w_to[i]])
    })
    
    B_bifidum_coverage_MGS <- cbind(B_bifidum_coverage_MGS, V)
    colnames(B_bifidum_coverage_MGS)[ncol(B_bifidum_coverage_MGS)] <- sample_id
     lines(x=w_center, 
           y=V+1,
           col=MGS_metadata[MGS_metadata$NG_ID==sample_id,]$Color 
           )

  }
  
  # df to keep mean coverage over window
  B_bifidum_coverage_VLP <- data.frame(w_center)
  
  for (sample_id in VLP_metadata[VLP_metadata$L34922_LS1_positive=="YES",]$NG_ID) {
    
    tab <- read.table(paste0('~/Desktop/Projects_2022/NEXT_pilot_FUP/03.SCRIPTS/NEXT_pilot_FUP_bf_origin/L34922_LS1/VLP_alignments/', sample_id, '.target.cov.txt'), sep='\t', header=F, fill=T)
    
    COVERAGE <- rep(0, len)
    COVERAGE[tab$V2] <- tab$V4
    
    V <- c()
    V <- sapply(seq_along(w_center), function(i) {
      mean( COVERAGE[ w_from[i]:w_to[i]])
    })
    
    B_bifidum_coverage_VLP <- cbind(B_bifidum_coverage_VLP, V)
    colnames(B_bifidum_coverage_VLP)[ncol(B_bifidum_coverage_VLP)] <- sample_id
    lines(x=w_center, 
          y=V+1,
          col=VLP_metadata[VLP_metadata$NG_ID==sample_id,]$Color 
    )
    
  }
  
}


# Visualization with ggplot2:
melted_df <- melt(B_bifidum_coverage_MGS, id.vars = 'w_center')

# adding 1 so that we can log-scale y-axis
melted_df$value <- melted_df$value + 1
melted_df$FAM_ID <- MGS_metadata$FAM_ID[match(melted_df$variable, MGS_metadata$NG_ID)]
melted_df$Timepoint <- MGS_metadata$Timepoint[match(melted_df$variable, MGS_metadata$NG_ID)]
melted_df$Type <- MGS_metadata$Type[match(melted_df$variable, MGS_metadata$NG_ID)]
melted_df$color_factor <- paste0(melted_df$Type, ' ', melted_df$Timepoint)
melted_df[melted_df$FAM_ID!='FAM0234',]$color_factor <- "Other"
melted_df$width <- melted_df$FAM_ID
melted_df[melted_df$width!="FAM0234",]$width <- "Other"
melted_df$color_factor <- factor(melted_df$color_factor, 
                                 levels = c("Other", "Mother P7", "Mother B", 
                                            "Mother M1", "Mother M2", "Mother M3", 
                                            "Infant M1", "Infant M2", "Infant M3",
                                            "Infant M6", "Infant M9"), 
                              ordered = T)
read_align_colors <- data.frame(Type_time=unique(melted_df$color_factor), Color=c('#E9E8E8', "#11468f", "#2d5d9d", 
                                                             "#4175aa", "#518eb8", "#61a8c6", 
                                                             "#da1212", "#da4635", "#d86555",
                                                             "#d28076", "#c79898"))

prophage_region <- data.frame(coordinate=c(54887, 90884), orientation=c('start', 'end'))

ggplot(melted_df, aes(w_center, value, group=variable, color=color_factor, size=width)) +
  geom_line() + 
  xlim(0, 150000) + 
  labs(x="B. bifidum patched genome fragment", y="Mean MGS read depth (log10)", color="Timepoint") +
  scale_y_log10() + 
  scale_size_manual(values=c(1, 0.5), guide=FALSE) +
  scale_color_manual(values=read_align_colors$Color) + 
  theme_bw() + 
  theme(legend.position = "bottom") + 
  guides(color=guide_legend(nrow=3, title.position="top", title.hjust = 0.5)) + 
  ggnewscale::new_scale_color() +
  geom_vline(data=prophage_region, aes(xintercept=coordinate, color="Prophage region"), size=1, linetype="dashed") + 
  scale_color_manual(name="",
                     values=c(`Prophage region`="black"), 
                     guide = guide_legend(order = 3)) 


# Visualization with ggplot2:
melted_df_VLP <- melt(B_bifidum_coverage_VLP, id.vars = 'w_center')

# adding 1 so that we can log-scale y-axis
melted_df_VLP$value <- melted_df_VLP$value + 1
melted_df_VLP$FAM_ID <- VLP_metadata$FAM_ID[match(melted_df_VLP$variable, VLP_metadata$NG_ID)]
melted_df_VLP$Timepoint <- VLP_metadata$Timepoint[match(melted_df_VLP$variable, VLP_metadata$NG_ID)]
melted_df_VLP$Type <- VLP_metadata$Type[match(melted_df_VLP$variable, VLP_metadata$NG_ID)]
melted_df_VLP$color_factor <- paste0(melted_df_VLP$Type, '_', melted_df_VLP$Timepoint)
melted_df_VLP[melted_df_VLP$FAM_ID!='FAM0234',]$color_factor <- "Other"
melted_df_VLP$width <- melted_df_VLP$FAM_ID
melted_df_VLP[melted_df_VLP$width!="FAM0234",]$width <- "Other"
melted_df_VLP$color_factor <- factor(melted_df_VLP$color_factor, 
                                 levels = c("Other", "Mother_P7", "Mother_B", 
                                            "Mother_M1", "Mother_M2", "Mother_M3", 
                                            "Infant_M2", "Infant_M3"), 
                                 ordered = T)
read_align_colors <- data.frame(Type_time=unique(melted_df$color_factor), Color=c('#E9E8E8', "#11468f", "#2d5d9d", 
                                                                                  "#4175aa", "#518eb8", "#61a8c6", 
                                                                                  "#da1212", "#da4635", "#d86555",
                                                                                  "#d28076", "#c79898"))


ggplot(melted_df_VLP, aes(w_center, value, group=variable, color=color_factor, size=width)) +
  geom_line() + 
  xlim(0, 150000) + 
  labs(x="B. bifidum patched genome fragment", y="Mean VLP read depth (log10)", color="Timepoint") +
  scale_y_log10() + 
  scale_size_manual(values=c(1, 0.5), guide=FALSE) +
  scale_color_manual(values=read_align_colors$Color) + 
  theme_bw() + 
  theme(legend.position = "bottom") + 
  guides(color=guide_legend(nrow=3, title.position="top", title.hjust = 0.5)) + 
  ggnewscale::new_scale_color() +
  geom_vline(data=prophage_region, aes(xintercept=coordinate, color="Prophage region"), size=1, linetype="dashed") + 
  scale_color_manual(name="",
                     values=c(`Prophage region`="black"), 
                     guide = guide_legend(order = 3)) 

# detect samples where prophage is absent at the indicated region (MGS): 
for (seq_id in c('MGYG000132487_6_RagTag')) {
  
  # df to keep all bases of the prophage region per MGS sample
  prophage_absence_MGS <- data.frame(V2=c(54887:90884))
  
  for (sample_id in MGS_metadata[MGS_metadata$B_bifidum_positive=="YES",]$NG_ID) {
    
    tab <- read.table(paste0('~/Desktop/Projects_2022/NEXT_pilot_FUP/03.SCRIPTS/NEXT_pilot_FUP_bf_origin/L34922_LS1/MGS_alignments/', sample_id, '.target.cov.txt'), sep='\t', header=F, fill=T)
    
    prophage_absence_MGS <- merge(prophage_absence_MGS, tab[,c('V2', 'V4')], by='V2', all.x=T)
    
    colnames(prophage_absence_MGS)[ncol(prophage_absence_MGS)] <- sample_id

  }
  
  # df to keep all bases of the prophage region per VLP sample
  prophage_absence_VLP <- data.frame(V2=c(54887:90884))
  
  for (sample_id in VLP_metadata[VLP_metadata$L34922_LS1_positive=="YES",]$NG_ID) {
    
    tab <- read.table(paste0('~/Desktop/Projects_2022/NEXT_pilot_FUP/03.SCRIPTS/NEXT_pilot_FUP_bf_origin/L34922_LS1/VLP_alignments/', sample_id, '.target.cov.txt'), sep='\t', header=F, fill=T)
    
    prophage_absence_VLP <- merge(prophage_absence_VLP, tab[,c('V2', 'V4')], by='V2', all.x=T)
    
    colnames(prophage_absence_VLP)[ncol(prophage_absence_VLP)] <- sample_id
    
  }
  
}

prophage_absence_MGS[is.na(prophage_absence_MGS)] <- 0
prophage_absence_VLP[is.na(prophage_absence_VLP)] <- 0


length(which(colSums(prophage_absence_MGS!=0)/35997 < 0.75))/213
which(MGS_metadata[MGS_metadata$Color!='grey' & MGS_metadata$B_bifidum_positive=='YES',"NG_ID"]  %in% names(which(colSums(prophage_absence_MGS!=0)/35997 >= 0.75)))

prophage_absence_MGS[,colSums(prophage_absence_MGS!=0)/35997 < 0.75 ]

length(which(colSums(prophage_absence_VLP!=0)/35997 < 0.7))
which(VLP_metadata[VLP_metadata$Color!='grey' & VLP_metadata$L34922_LS1_positive=='YES',"NG_ID"]  %in% names(which(colSums(prophage_absence_VLP!=0)/35997 >= 0.75)))

prophage_absence_MGS[,colSums(prophage_absence_MGS!=0)/35997 < 0.75 ]

###### FOR VISUALIZATION #####
#write.table(ali, "02.CLEAN_DATA/PREPARED_DATA_FOR_PLOTS/Fig6C_L34922_LS1_to_B.bifidum_patch.paf")
write.table(melted_df, "02.CLEAN_DATA/PREPARED_DATA_FOR_PLOTS/Fig6E_MGS_read_alignment_patched_B.bifidum.txt", sep='\t', row.names=F)
write.table(read_align_colors, "02.CLEAN_DATA/PREPARED_DATA_FOR_PLOTS/Fig6EF_read_alignment_colors.txt", sep='\t', row.names=F)
write.table(melted_df_VLP, "02.CLEAN_DATA/PREPARED_DATA_FOR_PLOTS/Fig6F_VLP_read_alignment_patched_B.bifidum.txt", sep='\t', row.names=F)
write.table(prophage_region, '02.CLEAN_DATA/PREPARED_DATA_FOR_PLOTS/Fig6EF_prophage_region_in_B.bifidum.txt', sep='\t', row.names = F)
