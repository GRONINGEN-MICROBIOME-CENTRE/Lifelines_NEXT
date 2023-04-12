setwd('~/Desktop/Projects_2022/NEXT_pilot_FUP/')

#############################################################
# Here we explore how viruses are shared between mothers and
# their infants
#############################################################

##############################
# Functions
##############################

# calculate_sharedness finds the % of shared entities between 
# the timepoint samples of infants and post and pre timepoints
# of mothers. It is tuned exactly for the format of dataframes
# generated for the NEXT pilot.

calculate_sharedness <- function(counts_df, metadata_df, preB_timepoints, postB_timepoints) {
  
  # Create empty data frame to store sharedness values
  sharedness_timepoints <- as.data.frame(matrix(NA, 
                                                nrow = length(unique(metadata_df[metadata_df$Type=='Infant',]$Individual_ID)), 
                                                ncol = length(unique(metadata_df[metadata_df$Type=='Infant',]$Timepoint))*2 ))
  
  row.names(sharedness_timepoints) <- unique(metadata_df[metadata_df$Type=='Infant',]$Individual_ID)
  
  colnames(sharedness_timepoints) <- c(paste0('to_preB', '_', unique(metadata_df[metadata_df$Type=='Infant',]$Timepoint)), 
                                       paste0('to_postB', '_', unique(metadata_df[metadata_df$Type=='Infant',]$Timepoint)))
  
  for (j in unique(metadata_df[metadata_df$Type=='Infant',]$Timepoint) ) {
    
    preB <- preB_timepoints # pre-birth timepoints
    postB <- postB_timepoints # post-birth timepoints
    
    for (i in row.names(sharedness_timepoints)) {
      
      FAMILY <- unique(metadata_df[metadata_df$Individual_ID==i,]$FAM_ID)
      
      N_preB_mother <- length(metadata_df[metadata_df$FAM_ID==FAMILY & metadata_df$Type=='Mother' & metadata_df$Timepoint %in% preB,]$Short_sample_ID)
      
      if (N_preB_mother!=0 & !is.na(any(metadata_df==metadata_df[metadata_df$Timepoint==j & metadata_df$Individual_ID==i,]$Short_sample_ID))) {
        A <- rowSums(counts_df[,c(metadata_df[metadata_df$FAM_ID==FAMILY &
                                                metadata_df$Type=='Mother' & 
                                                metadata_df$Timepoint %in% preB, ]$Short_sample_ID,'MOCK')])>0
        B <- rowSums(counts_df[,c(metadata_df[metadata_df$Individual_ID==i & 
                                                metadata_df$Type=='Infant' &
                                                metadata_df$Timepoint==j,]$Short_sample_ID,'MOCK')])>0
        sharedness_timepoints[i,paste0('to_preB', '_', j)] <- sum( rowSums(cbind(A,B)) ==2)/metadata_df[metadata_df$Individual_ID==i & metadata_df$Timepoint==j,]$viral_richness
      } else {
        sharedness_timepoints[i,paste0('to_preB', '_', j)] <- NA
      }
      
      
      N_postB_mother <- length(metadata_df[metadata_df$FAM_ID==FAMILY & metadata_df$Type=='Mother' & metadata_df$Timepoint %in% postB,]$Short_sample_ID)
      
      if (N_postB_mother!=0 & !is.na(any(metadata_df==metadata_df[metadata_df$Timepoint==j & metadata_df$Individual_ID==i,]$Short_sample_ID)) ) {
        A <- rowSums(counts_df[,c(metadata_df[metadata_df$FAM_ID==FAMILY &
                                                metadata_df$Type=='Mother' & 
                                                metadata_df$Timepoint %in% postB, ]$Short_sample_ID,'MOCK')])>0
        B <- rowSums(counts_df[,c(metadata_df[metadata_df$Individual_ID==i & 
                                                metadata_df$Type=='Infant' &
                                                metadata_df$Timepoint==j,]$Short_sample_ID, 'MOCK')])>0
        sharedness_timepoints[i,paste0('to_postB', '_', j)] <- sum( rowSums(cbind(A,B)) ==2)/metadata_df[metadata_df$Individual_ID==i & metadata_df$Timepoint==j,]$viral_richness
      } else {
        sharedness_timepoints[i,paste0('to_postB', '_', j)] <- NA
      }
      
    }
  }
  
  return(sharedness_timepoints)
}


##############################
# Loading libraries
##############################
library(stringr)
library(reshape2)
library(ggplot2)
library(ggforce)

library(lme4)
library(RLRsim)
library(lmerTest)
##############################
# Input data
##############################

VLP_metadata <- read.table("02.CLEAN_DATA/VLP_metadata_with_phenos.txt", sep='\t', header=T)

MGS_metadata <- read.table("02.CLEAN_DATA/MGS_metadata_with_phenos.txt", sep='\t', header=T)
colnames(MGS_metadata)[grep('viral_richness_MGS', colnames(MGS_metadata))] <- 'viral_richness'

combo_metadata <- rbind( VLP_metadata[,c("Type","Timepoint","FAM_ID","Short_sample_ID","Individual_ID","Universal_fecal_ID", "viral_richness")],
                         MGS_metadata[,c("Type","Timepoint","FAM_ID","Short_sample_ID","Individual_ID","Universal_fecal_ID", "viral_richness")] )

combo_metadata <- combo_metadata[ (combo_metadata$Short_sample_ID %in% VLP_metadata[VLP_metadata$Type=='Infant',]$Short_sample_ID) |
                                    ((combo_metadata$Short_sample_ID %in% MGS_metadata[MGS_metadata$Type=='Mother',]$Short_sample_ID) &
                                       (combo_metadata$Universal_fecal_ID %in% VLP_metadata[VLP_metadata$Type=='Mother',]$Universal_fecal_ID)),
  
]

RPKM_counts_VLP <- read.table("02.CLEAN_DATA/RPKM_counts_VLP.txt", sep='\t', header=T)


RPKM_counts_MGS <- read.table("02.CLEAN_DATA/RPKM_counts_MGS.txt", sep='\t', header=T)



RPKM_counts_combo <- merge(RPKM_counts_VLP, RPKM_counts_MGS, by='row.names', all = T)
row.names(RPKM_counts_combo) <- RPKM_counts_combo$Row.names
RPKM_counts_combo$Row.names <- NULL
RPKM_counts_combo[is.na(RPKM_counts_combo)] <- 0
RPKM_counts_combo <- RPKM_counts_combo[rowSums(RPKM_counts_combo)>0,]

contigs_metadata <- read.table('02.CLEAN_DATA/VLP_viral_contigs_metadata.txt', sep='\t', header=T)

##############################
# ANALYSIS
##############################

######### PURE VLP STAT #########
RPKM_counts_VLP$MOCK <- 0

sharedness_timepoints_VLP <- calculate_sharedness(RPKM_counts_VLP, VLP_metadata, c('P7','B'), c('M1','M3'))

sharedness_timepoints_VLP$ID <- row.names(sharedness_timepoints_VLP)
sharedness_timepoints_melt <- melt(sharedness_timepoints_VLP)
sharedness_timepoints_melt$Timepoint <- gsub('.*B_', '', sharedness_timepoints_melt$variable)
sharedness_timepoints_melt$Timepoint <- factor(sharedness_timepoints_melt$Timepoint, levels=c('M1', 'M2', 'M3', 'M6', "M12"), ordered = T)
sharedness_timepoints_melt$variable <- sub('B_.*', '', sharedness_timepoints_melt$variable)
sharedness_timepoints_melt$variable <- factor(sharedness_timepoints_melt$variable, levels=c('to_pre', 'to_post'), ordered=T)

pdf('./04.PLOTS/Infant_Viruses_shared_VLP_P7B_M1M3.pdf', width=12/2.54, height=9/2.54)
ggplot(sharedness_timepoints_melt, aes(Timepoint, value, fill=variable)) +
  geom_boxplot(outlier.shape = NA) +
  labs (y="% Infant viruses shared", x="Infant Timepoint") + 
  geom_sina(aes(fill=variable), size=0.6,alpha=0.5) +
  theme_bw()+
  theme(axis.text=element_text(size=12), 
        axis.title=element_text(size=16,face="bold"),
        strip.text.x = element_text(size = 12),
        legend.text = element_text(size=10),
        legend.title = element_text(size=12, face="bold")) +
  scale_fill_manual(name = "Maternal\nTimepoints", 
                    labels=c('Pre birth', 'Post birth'),
                    values=c("#FAE3D2", "#BBDED2"))
dev.off()

##### FORMAL TEST ######
sharedness_timepoints_stat <- sharedness_timepoints_melt[!is.na(sharedness_timepoints_melt$value),]
sharedness_timepoints_stat$coded_tp <- sharedness_timepoints_stat$Timepoint
sharedness_timepoints_stat$coded_tp <- str_pad(sub('M','',sharedness_timepoints_stat$coded_tp), 2, pad=0)
sharedness_timepoints_stat$Short_sample_ID <- paste0(substr(sharedness_timepoints_stat$ID, 1,1), 
                                                     sharedness_timepoints_stat$coded_tp, 
                                                     substr(sharedness_timepoints_stat$ID, 2,5),
                                                     'V')
sharedness_timepoints_stat <- merge(sharedness_timepoints_stat, VLP_metadata, by='Short_sample_ID')


t.test(sharedness_timepoints_melt[sharedness_timepoints_melt$Timepoint %in% c('M1',"M2",'M3') & 
                                    sharedness_timepoints_melt$variable=='to_pre',]$value,
       sharedness_timepoints_melt[sharedness_timepoints_melt$Timepoint %in% c('M1',"M2",'M3') & 
                                    sharedness_timepoints_melt$variable=='to_post',]$value, 
       paired=T)

# dependent on time?
mod0 <- lm(value ~ Age_days + DNA_CONC + Clean_reads, data = sharedness_timepoints_stat)
mod1  = lmer(value ~ Age_days + DNA_CONC + Clean_reads + (1|NEXT_ID), REML = F, data = sharedness_timepoints_stat)
BIC(mod0, mod1)
exactLRT(mod1,mod0)
summary(mod1) #Age_days 2.390e-04  1.084e-04  1.523e+02   2.204   0.0290 *
as.data.frame(summary(mod1)$coefficients)[,1:5]

# different between post and pre?
btmod0 <- lm(value ~ variable + Age_days + DNA_CONC + Clean_reads, data = sharedness_timepoints_stat)
btmod1  = lmer(value ~ variable + Age_days + DNA_CONC + Clean_reads + (1|NEXT_ID), REML = F, data = sharedness_timepoints_stat)
BIC(btmod0, btmod1)
exactLRT(btmod1,btmod0)
summary(btmod0) 
summary(btmod1) #variable    3.346e-02  1.637e-02  1.362e+02   2.043   0.0429 *

######### PURE MGS STAT #########
RPKM_counts_MGS$MOCK <- 0

sharedness_timepoints_MGS <- calculate_sharedness(RPKM_counts_MGS, MGS_metadata, c('P3','P7', 'B'), c('M1', 'M2','M3'))

sharedness_timepoints_MGS$ID <- row.names(sharedness_timepoints_MGS)
sharedness_timepoints_MGS_melt <- melt(sharedness_timepoints_MGS)
sharedness_timepoints_MGS_melt$Timepoint <- gsub('.*B_', '', sharedness_timepoints_MGS_melt$variable)
sharedness_timepoints_MGS_melt$Timepoint <- factor(sharedness_timepoints_MGS_melt$Timepoint, levels=c('M1', 'M2', 'M3', 'M6', 'M9', "M12"), ordered = T)
sharedness_timepoints_MGS_melt$variable <- sub('B_.*', '', sharedness_timepoints_MGS_melt$variable)
sharedness_timepoints_MGS_melt$variable <- factor(sharedness_timepoints_MGS_melt$variable, levels=c('to_pre', 'to_post'), ordered=T)

pdf('./04.PLOTS/Infant_Viruses_shared_MGS.pdf', width=12/2.54, height=9/2.54)
ggplot(sharedness_timepoints_MGS_melt, aes(Timepoint, value, fill=variable)) +
  geom_boxplot(outlier.shape = NA) +
  labs (y="% Infant viruses shared", x="Infant Timepoint") + 
  geom_sina(aes(fill=variable), size=0.6, alpha=0.5) +
  theme_bw()+
  theme(axis.text=element_text(size=12), 
        axis.title=element_text(size=16,face="bold"),
        strip.text.x = element_text(size = 12),
        legend.text = element_text(size=10),
        legend.title = element_text(size=12, face="bold")) +
  scale_fill_manual(name = "Maternal\nTimepoints", 
                    labels=c('Pre birth', 'Post birth'),
                    values=c("#FAE3D2", "#BBDED2"))
dev.off()

######### INFANT VLP TO MATERNAL MGS STAT #########
RPKM_counts_combo$MOCK <- 0

sharedness_timepoints_iV_mM <- calculate_sharedness(RPKM_counts_combo, combo_metadata, c('P7', 'B'), c('M1', 'M3'))

sharedness_timepoints_iV_mM$ID <- row.names(sharedness_timepoints_iV_mM)
sharedness_timepoints_iV_mM_melt <- melt(sharedness_timepoints_iV_mM)
sharedness_timepoints_iV_mM_melt$Timepoint <- gsub('.*B_', '', sharedness_timepoints_iV_mM_melt$variable)
sharedness_timepoints_iV_mM_melt$Timepoint <- factor(sharedness_timepoints_iV_mM_melt$Timepoint, levels=c('M1', 'M2', 'M3', 'M6',  "M12"), ordered = T)
sharedness_timepoints_iV_mM_melt$variable <- sub('B_.*', '', sharedness_timepoints_iV_mM_melt$variable)
sharedness_timepoints_iV_mM_melt$variable <- factor(sharedness_timepoints_iV_mM_melt$variable, levels=c('to_pre', 'to_post'), ordered=T)

pdf('./04.PLOTS/Infant_Viruses_shared_iVLP_mMGS.pdf', width=12/2.54, height=9/2.54)
ggplot(sharedness_timepoints_iV_mM_melt, aes(Timepoint, value, fill=variable)) +
  geom_boxplot(outlier.shape = NA) +
  labs (y="% Infant viruses shared", x="Infant Timepoint") + 
  geom_sina(aes(fill=variable), size=0.6,alpha=0.5) +
  theme_bw()+
  theme(axis.text=element_text(size=12), 
        axis.title=element_text(size=16,face="bold"),
        strip.text.x = element_text(size = 12),
        legend.text = element_text(size=10),
        legend.title = element_text(size=12, face="bold")) +
  scale_fill_manual(name = "Maternal\nTimepoints", 
                    labels=c('Pre birth', 'Post birth'),
                    values=c("#FAE3D2", "#BBDED2"))
dev.off()

######### INFANT VLP TO MATERNAL VLP ENRICHNED FOR MGS-DETECTED PROPHAGES STAT #########
RPKM_counts_VLP_plus_prophages <- RPKM_counts_VLP


RPKM_counts_MGS_prophages_only <- RPKM_counts_MGS[row.names(RPKM_counts_MGS) %in% contigs_metadata[contigs_metadata$temperate==1,]$V1,]


for (i in VLP_metadata[VLP_metadata$Type=='Mother',]$Universal_fecal_ID) {
  
  if ( any(MGS_metadata$Short_sample_ID == paste0(i,'M')) ) {
    prophages_detected_in_MGS <- row.names(RPKM_counts_MGS_prophages_only[ RPKM_counts_MGS_prophages_only[,paste0(i,'M')]!=0 , ])
    
    RPKM_counts_VLP_plus_prophages[prophages_detected_in_MGS,paste0(i,'V')] <- 1
  }
  
}

sharedness_timepoints_iV_mVM <- calculate_sharedness(RPKM_counts_VLP_plus_prophages, VLP_metadata, c('P7', 'B'), c('M1', 'M3'))

sharedness_timepoints_iV_mVM$ID <- row.names(sharedness_timepoints_iV_mVM)
sharedness_timepoints_iV_mVM_melt <- melt(sharedness_timepoints_iV_mVM)
sharedness_timepoints_iV_mVM_melt$Timepoint <- gsub('.*B_', '', sharedness_timepoints_iV_mVM_melt$variable)
sharedness_timepoints_iV_mVM_melt$Timepoint <- factor(sharedness_timepoints_iV_mVM_melt$Timepoint, levels=c('M1', 'M2', 'M3', 'M6',  "M12"), ordered = T)
sharedness_timepoints_iV_mVM_melt$variable <- sub('B_.*', '', sharedness_timepoints_iV_mVM_melt$variable)
sharedness_timepoints_iV_mVM_melt$variable <- factor(sharedness_timepoints_iV_mVM_melt$variable, levels=c('to_pre', 'to_post'), ordered=T)


pdf('./04.PLOTS/Infant_Viruses_shared_iVLP_mVLP_PRO.pdf', width=12/2.54, height=9/2.54)
ggplot(sharedness_timepoints_iV_mVM_melt, aes(Timepoint, value, fill=variable)) +
  geom_boxplot(outlier.shape = NA) +
  labs (y="% Infant viruses shared", x="Infant Timepoint") + 
  geom_sina(aes(fill=variable), size=0.6,alpha=0.5) +
  theme_bw()+
  theme(axis.text=element_text(size=12), 
        axis.title=element_text(size=16,face="bold"),
        strip.text.x = element_text(size = 12),
        legend.text = element_text(size=10),
        legend.title = element_text(size=12, face="bold")) +
  scale_fill_manual(name = "Maternal\nTimepoints", 
                    labels=c('Pre birth', 'Post birth'),
                    values=c("#FAE3D2", "#BBDED2"))
dev.off()

##### DOES TAKING INTO ACCOUNT PROPHAGES INCREASE SHAREDNESS?
comparison_VLP_VLP_pro <- sharedness_timepoints_melt
comparison_VLP_VLP_pro$source <- 'VLP'
buffer <- sharedness_timepoints_iV_mVM_melt
buffer$source <- 'VLP_plus_prophages'

comparison_VLP_VLP_pro <- rbind(comparison_VLP_VLP_pro, buffer)

t.test(comparison_VLP_VLP_pro$value ~ comparison_VLP_VLP_pro$source, paired=T)

pdf('./04.PLOTS/Difference_infant_virus_shared_VLP_vs_VLP_PRO.pdf', width=12/2.54, height=9/2.54)
ggplot(comparison_VLP_VLP_pro, aes(source, value)) +
  geom_boxplot() + 
  labs (y="% Infant viruses shared", x="Type of data") + 
  geom_sina(aes(fill=source), size=0.6,alpha=0.5) +
  theme_bw()+
  theme(axis.text=element_text(size=12), 
        axis.title=element_text(size=16,face="bold"),
        strip.text.x = element_text(size = 12),
        legend.text = element_text(size=10),
        legend.title = element_text(size=12, face="bold"),
        legend.position = 'none') + 
  annotate(geom="text", x=1.5, y=0.83, 
           label="p-value = 7.674e-14\nmedian difference 0.007")
dev.off()

######### INFANT VLP + prophages TO MATERNAL VLP + prophages (ENRICHNED FOR MGS-DETECTED PROPHAGES STAT) #########
RPKM_counts_VLP_plus_prophages_all <- RPKM_counts_VLP


RPKM_counts_MGS_prophages_only <- RPKM_counts_MGS[row.names(RPKM_counts_MGS) %in% contigs_metadata[contigs_metadata$temperate==1,]$V1,]


for (i in VLP_metadata$Universal_fecal_ID) {
  
  if ( any(MGS_metadata$Short_sample_ID == paste0(i,'M')) ) {
    prophages_detected_in_MGS <- row.names(RPKM_counts_MGS_prophages_only[ RPKM_counts_MGS_prophages_only[,paste0(i,'M')]!=0 , ])
    
    RPKM_counts_VLP_plus_prophages_all[prophages_detected_in_MGS,paste0(i,'V')] <- 1
  }
  
}

VLP_metadata_pro <- VLP_metadata
RPKM_counts_VLP_plus_prophages_all <- RPKM_counts_VLP_plus_prophages_all[,VLP_metadata_pro$Short_sample_ID]
identical(VLP_metadata_pro$Short_sample_ID, colnames(RPKM_counts_VLP_plus_prophages_all) )
VLP_metadata_pro$viral_richness <- colSums(RPKM_counts_VLP_plus_prophages_all>0)
RPKM_counts_VLP_plus_prophages_all$MOCK <- 0

sharedness_timepoints_iVM_mVM <- calculate_sharedness(RPKM_counts_VLP_plus_prophages_all, VLP_metadata_pro, c('P7', 'B'), c('M1', 'M3'))

sharedness_timepoints_iVM_mVM$ID <- row.names(sharedness_timepoints_iVM_mVM)
sharedness_timepoints_iVM_mVM_melt <- melt(sharedness_timepoints_iVM_mVM)
sharedness_timepoints_iVM_mVM_melt$Timepoint <- gsub('.*B_', '', sharedness_timepoints_iVM_mVM_melt$variable)
sharedness_timepoints_iVM_mVM_melt$Timepoint <- factor(sharedness_timepoints_iVM_mVM_melt$Timepoint, levels=c('M1', 'M2', 'M3', 'M6',  "M12"), ordered = T)
sharedness_timepoints_iVM_mVM_melt$variable <- sub('B_.*', '', sharedness_timepoints_iVM_mVM_melt$variable)
sharedness_timepoints_iVM_mVM_melt$variable <- factor(sharedness_timepoints_iVM_mVM_melt$variable, levels=c('to_pre', 'to_post'), ordered=T)


pdf('./04.PLOTS/Infant_Viruses_shared_iVLP_PRO_mVLP_PRO.pdf', width=12/2.54, height=9/2.54)
ggplot(sharedness_timepoints_iVM_mVM_melt, aes(Timepoint, value, fill=variable)) +
  geom_boxplot(outlier.shape = NA) +
  labs (y="% Infant viruses shared", x="Infant Timepoint") + 
  geom_sina(aes(fill=variable), size=0.6,alpha=0.5) +
  theme_bw()+
  theme(axis.text=element_text(size=12), 
        axis.title=element_text(size=16,face="bold"),
        strip.text.x = element_text(size = 12),
        legend.text = element_text(size=10),
        legend.title = element_text(size=12, face="bold")) +
  scale_fill_manual(name = "Maternal\nTimepoints", 
                    labels=c('Pre birth', 'Post birth'),
                    values=c("#FAE3D2", "#BBDED2"))
dev.off()

##############################
# OUTPUT
##############################



