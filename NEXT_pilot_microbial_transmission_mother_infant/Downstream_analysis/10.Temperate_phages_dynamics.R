setwd('~/Desktop/Projects_2022/NEXT_pilot_FUP/')

#############################################################
# Here we explore temperate phages of in gut microbiome and 
# virome
#############################################################

##############################
# Functions
##############################


##############################
# Loading libraries
##############################
library(lme4)
library(RLRsim)
library(lmerTest)
library(reshape2)
library(ggplot2)
library(ggforce)
##############################
# Input data
##############################

VLP_metadata <- read.table("02.CLEAN_DATA/VLP_metadata_with_phenos.txt", sep='\t', header=T, row.names = "Short_sample_ID")
VLP_metadata$Timepoint <- factor(VLP_metadata$Timepoint, levels = c("P7", "B",
                                                                    "M1", "M2", "M3",
                                                                    "M6", "M12"), ordered = T)
VLP_metadata$Type <- as.factor(VLP_metadata$Type)
VLP_metadata$Short_sample_ID <- row.names(VLP_metadata)
VLP_metadata$temperate_perc <- VLP_metadata$temperate_richness/VLP_metadata$viral_richness

MGS_metadata <- read.table("02.CLEAN_DATA/MGS_metadata_with_phenos.txt", sep='\t', header=T, row.names = "Short_sample_ID")
MGS_metadata$Timepoint <- factor(MGS_metadata$Timepoint, levels = c("P3", "P7", "B",
                                                                    "M1", "M2", "M3",
                                                                    "M6", "M9", "M12"), ordered = T)
MGS_metadata$Type <- as.factor(MGS_metadata$Type)
MGS_metadata$Short_sample_ID <- row.names(MGS_metadata)


##############################
# ANALYSIS
##############################
# is there a difference in richness between mothers and infants
type_mod0 <- lm(temperate_richness ~ Type + DNA_CONC + Clean_reads, data = VLP_metadata)
type_mod1  = lmer(temperate_richness ~ Type + DNA_CONC + Clean_reads + (1|NEXT_ID), REML = F, data = VLP_metadata)
BIC(type_mod0, type_mod1)
exactLRT(type_mod1,type_mod0)
summary(type_mod1) #Type 3.777e+02  5.980e+01  6.565e+01   6.315 2.67e-08 ***

# dependant on time in babies? (switching to the exact ages for a higher precision)
btmod0 <- lm(temperate_richness ~ Age_days + DNA_CONC + Clean_reads, data = VLP_metadata[VLP_metadata$Type=="Infant",])
btmod1  = lmer(temperate_richness ~ Age_days + DNA_CONC + Clean_reads + (1|NEXT_ID), REML = F, data = VLP_metadata[VLP_metadata$Type=="Infant",])
BIC(btmod0, btmod1)
exactLRT(btmod1,btmod0)
summary(btmod0) 
summary(btmod1) #Age_days     2.573e-01  4.817e-02  7.948e+01   5.342 8.55e-07 ***

# dependant on time in mothers? NS
mtmod0 <- lm(temperate_richness ~ Timepoint_continuous+ DNA_CONC + Clean_reads, data = VLP_metadata[VLP_metadata$Type=="Mother",])
mtmod1  = lmer(temperate_richness ~ Timepoint_continuous+ DNA_CONC + Clean_reads + (1|NEXT_ID), REML = F, data = VLP_metadata[VLP_metadata$Type=="Mother",])
BIC(mtmod0, mtmod1)
exactLRT(mtmod1,mtmod0)
summary(mtmod1) #Timepoint_continuous -3.816e-01  1.255e+01  9.017e+01  -0.030  0.97580 

# Richness of active temperate phages
pdf('./04.PLOTS/Figure3A_temperate_richness_plot.pdf', width=12/2.54, height=9/2.54)
ggplot(VLP_metadata, aes(Timepoint, temperate_richness/viral_richness, fill=Type)) + 
  labs (y="Richness of active\ntemperate phages (log10)", x="") + 
  geom_violin() + 
  geom_boxplot(width = 0.1, outlier.shape = NA) + 
  geom_sina(colour='#303841', size=0.6) +
  theme_bw()+
  theme(legend.position="none") + 
  facet_grid(~Type, scales="free", space = "free") +
  scale_y_continuous(trans='log10') +
  theme(axis.text=element_text(size=12), 
        axis.title=element_text(size=16,face="bold"),
        strip.text.x = element_text(size = 12),
        legend.text = element_text(size=10),
        legend.title = element_text(size=12, face="bold")) +
  labs(colour = "Type")
dev.off()


# dependent on phenotypes in infants? 
btmod1_feeding  = lmer(temperate_richness ~ Age_days + infant_ever_never_breastfed + viral_richness + DNA_CONC + Clean_reads + (1|NEXT_ID), REML = F, data = VLP_metadata[VLP_metadata$Type=="Infant",])
summary(btmod1_feeding) #infant_ever_never_breastfed  4.674e+01  1.271e+01  7.300e+01   3.678 0.000446 ***

btmod1_feeding  = lmer(temperate_richness/viral_richness ~ Age_days + DNA_CONC + Clean_reads + (1|NEXT_ID), REML = F, data = VLP_metadata[VLP_metadata$Type=="Infant",])
summary(btmod1_feeding) #infant_ever_never_breastfed  4.674e+01  1.271e+01  7.300e+01   3.678 0.000446 ***


btmod1_feeding  = lmer(temperate_RA ~ Age_days + temperate_perc + DNA_CONC + Clean_reads + (1|NEXT_ID), REML = F, data = VLP_metadata[VLP_metadata$Type=="Infant",])
summary(btmod1_feeding)

summary(lm(data = VLP_metadata[VLP_metadata$Type=="Infant",], temperate_RA ~ temperate_perc))
cor.test(VLP_metadata[VLP_metadata$Type=="Infant",]$temperate_RA, VLP_metadata[VLP_metadata$Type=="Infant",]$temperate_perc, method='pearson')

tmp <- melt(VLP_metadata[VLP_metadata$Type=='Infant' & !is.na(VLP_metadata$infant_ever_never_breastfed),c("Timepoint","temperate_richness", "infant_ever_never_breastfed", "Short_sample_ID")])
tmp$infant_ever_never_breastfed <- as.factor(tmp$infant_ever_never_breastfed)

ggplot(tmp, aes(Timepoint, value, fill=infant_ever_never_breastfed)) + 
  labs (y="Richness of active\ntemperate phages (log10)", x="") + 
  geom_boxplot(outlier.shape = NA) + 
  geom_sina(colour='#303841', size=0.6) +
  theme_bw()+
  theme(axis.text=element_text(size=12), 
        axis.title=element_text(size=16,face="bold"),
        strip.text.x = element_text(size = 12),
        legend.text = element_text(size=10),
        legend.title = element_text(size=12, face="bold")) +
  scale_fill_manual(name = "Feeding mode", 
                    labels=c('Ever breastfed', 'Never breastfed'),
                    values=c("#FAE3D9", "#BBDED6"))


# is there a difference in the relative abundance of active temperate phages between mothers and infants?
type_mod0_RA <- lm(temperate_RA ~ Type + DNA_CONC + Clean_reads, data = VLP_metadata)
type_mod1_RA  = lmer(temperate_RA ~ Type + DNA_CONC + Clean_reads + (1|NEXT_ID), REML = F, data = VLP_metadata)
BIC(type_mod0_RA, type_mod1_RA)
exactLRT(type_mod1_RA,type_mod0_RA)
summary(type_mod1_RA) #Type -2.220e+01  2.923e+00  4.694e+01  -7.595 1.04e-09 ***

# dependant on time in babies? (switching to the exact ages for a higher precision)
btmod0_RA <- lm(temperate_RA ~ Age_days + DNA_CONC +  Clean_reads, data = VLP_metadata[VLP_metadata$Type=="Infant",])
btmod1_RA  = lmer(temperate_RA ~ Age_days + DNA_CONC + bacterial_alpha_diversity  + Clean_reads + (1|NEXT_ID), REML = F, data = VLP_metadata[VLP_metadata$Type=="Infant",])
BIC(btmod0_RA, btmod1_RA)
exactLRT(btmod1_RA,btmod0_RA)
summary(btmod0_RA) 
summary(btmod1_RA) #Age_days     2.573e-01  4.817e-02  7.948e+01   5.342 8.55e-07 ***

# dependant on time in mothers? NS
mtmod0 <- lm(temperate_richness ~ Timepoint_continuous+ DNA_CONC + Clean_reads, data = VLP_metadata[VLP_metadata$Type=="Mother",])
mtmod1  = lmer(temperate_richness ~ Timepoint_continuous+ DNA_CONC + Clean_reads + (1|NEXT_ID), REML = F, data = VLP_metadata[VLP_metadata$Type=="Mother",])
BIC(mtmod0, mtmod1)
exactLRT(mtmod1,mtmod0)
summary(mtmod1) #Timepoint_continuous -3.816e-01  1.255e+01  9.017e+01  -0.030  0.97580 





tmp1 <- VLP_metadata
tmp1$proportion_temperate <- tmp1$temperate_richness/tmp1$viral_richness
tmp2 <- melt(tmp1[tmp1$Type=='Infant' & !is.na(tmp1$infant_ever_never_breastfed),c("Timepoint", "proportion_temperate","infant_ever_never_breastfed", "Short_sample_ID")])
tmp2$infant_ever_never_breastfed <- as.factor(tmp$infant_ever_never_breastfed)

ggplot(tmp2, aes(Timepoint, value, fill=infant_ever_never_breastfed)) + 
  labs (y="Relative abundance of active\ntemperate phages", x="") + 
  geom_boxplot(outlier.shape = NA) + 
  geom_sina(colour='#303841', size=0.6) +
  theme_bw()+
  theme(axis.text=element_text(size=12), 
        axis.title=element_text(size=16,face="bold"),
        strip.text.x = element_text(size = 12),
        legend.text = element_text(size=10),
        legend.title = element_text(size=12, face="bold")) +
  scale_fill_manual(name = "Feeding mode", 
                    labels=c('Ever breastfed', 'Never breastfed'),
                    values=c("#FAE3D9", "#BBDED6"))





summary(lm(data=VLP_metadata[VLP_metadata$Type=='Infant',], temperate_RA ~ Timepoint_continuous + viral_alpha_diversity + infant_feeding_mode_imputed_B_to_M3))

ggplot(MGS_metadata, aes(Timepoint, temperate_richness_MGS, fill=Type)) + 
  labs (y="Richness of active\ntemperate phages (log10)", x="") + 
  geom_violin() + 
  geom_boxplot(width = 0.1, outlier.shape = NA) + 
  geom_sina(colour='#303841', size=0.6) +
  theme_bw()+
  theme(legend.position="none") + 
  facet_grid(~Type, scales="free", space = "free") +
  scale_y_continuous(trans='log10') +
  theme(axis.text=element_text(size=12), 
        axis.title=element_text(size=16,face="bold"),
        strip.text.x = element_text(size = 12),
        legend.text = element_text(size=10),
        legend.title = element_text(size=12, face="bold")) +
  labs(colour = "Type")

ggplot(MGS_metadata, aes(Timepoint, temperate_RA_MGS, fill=Type)) + 
  labs (y="Richness of active\ntemperate phages (log10)", x="") + 
  geom_violin() + 
  geom_boxplot(width = 0.1, outlier.shape = NA) + 
  geom_sina(colour='#303841', size=0.6) +
  theme_bw()+
  theme(legend.position="none") + 
  facet_grid(~Type, scales="free", space = "free") +
  theme(axis.text=element_text(size=12), 
        axis.title=element_text(size=16,face="bold"),
        strip.text.x = element_text(size = 12),
        legend.text = element_text(size=10),
        legend.title = element_text(size=12, face="bold")) +
  labs(colour = "Type")

boxplot(MGS_metadata$temperate_RA_MGS ~ MGS_metadata$Type + MGS_metadata$Timepoint)
boxplot(MGS_metadata$temperate_richness ~ MGS_metadata$Type + MGS_metadata$Timepoint)

boxplot(VLP_metadata$temperate_richness ~ VLP_metadata$Type + VLP_metadata$Timepoint)
summary(lm(data = VLP_metadata[VLP_metadata$Type=='Infant',], temperate_richness ~ Age_days))
summary(lm(data = VLP_metadata[VLP_metadata$Type=='Mother',], temperate_richness ~ Age_days))

btmod1_feeding_MGS  = lmer(temperate_RA_MGS ~ Type +  DNA_CONC + Clean_reads + (1|NEXT_ID), REML = F, data = MGS_metadata)
summary(btmod1_feeding_MGS) #Type  -2.117e+00  1.032e+00  8.003e+01  -2.052   0.0434 * 



