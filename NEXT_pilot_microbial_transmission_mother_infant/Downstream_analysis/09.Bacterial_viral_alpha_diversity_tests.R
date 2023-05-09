setwd('~/Desktop/Projects_2022/NEXT_pilot_FUP/')

#############################################################
# Here we explore changes of gut bacteriome and virome
# alpha-diversity
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
library(ggsignif)
##############################
# Input data
##############################

VLP_metadata <- read.table("02.CLEAN_DATA/VLP_metadata_with_phenos.txt", sep='\t', header=T, row.names = "Short_sample_ID")
VLP_metadata$Timepoint <- factor(VLP_metadata$Timepoint, levels = c("P7", "B",
                                                                    "M1", "M2", "M3",
                                                                    "M6", "M12"), ordered = T)
VLP_metadata$Type <- as.factor(VLP_metadata$Type)
VLP_metadata$Short_sample_ID <- row.names(VLP_metadata)

MGS_metadata <- read.table("02.CLEAN_DATA/MGS_metadata_with_phenos.txt", sep='\t', header=T, row.names = "Short_sample_ID")
MGS_metadata$Timepoint <- factor(MGS_metadata$Timepoint, levels = c("P3", "P7", "B",
                                                                    "M1", "M2", "M3",
                                                                    "M6", "M9", "M12"), ordered = T)
MGS_metadata$Type <- as.factor(MGS_metadata$Type)
MGS_metadata$Short_sample_ID <- row.names(MGS_metadata)

# alpha-diversity in virome

# dependent on type?
type_mod0 <- lm(viral_alpha_diversity ~ Type + DNA_CONC + Clean_reads, data = VLP_metadata)
type_mod1  = lmer(viral_alpha_diversity ~ Type + DNA_CONC + Clean_reads + (1|NEXT_ID), REML = F, data = VLP_metadata)
BIC(type_mod0, type_mod1)
exactLRT(type_mod1,type_mod0)
summary(type_mod1) #Type 2.525e+00  3.517e-01 6.182e+01   7.178 1.07e-09 ***


# dependent on type & time?
mod0 <- lm(viral_alpha_diversity ~ Timepoint_continuous + DNA_CONC + Clean_reads, data = VLP_metadata)
mod1  = lmer(viral_alpha_diversity ~ Type + Timepoint_continuous + DNA_CONC + Clean_reads + (1|NEXT_ID), REML = F, data = VLP_metadata)
BIC(mod0, mod1)
exactLRT(mod1,mod0)
summary(mod1) #Timepoint_continuous 1.132e-01  3.393e-02 1.768e+02   3.337 0.001033 **  
as.data.frame(summary(mod1)$coefficients)[,1:5]


# dependant on time in babies? (switching to the exact ages for a higher precision)
btmod0 <- lm(viral_alpha_diversity ~ Age_days + DNA_CONC + Clean_reads, data = VLP_metadata[VLP_metadata$Type=="Infant",])
btmod1  = lmer(viral_alpha_diversity ~ Age_days + DNA_CONC + Clean_reads + (1|NEXT_ID), REML = F, data = VLP_metadata[VLP_metadata$Type=="Infant",])
BIC(btmod0, btmod1)
exactLRT(btmod1,btmod0)
summary(btmod0) 
summary(btmod1) #Age_days     4.733e-03  1.284e-03  8.112e+01   3.687 0.000408 ***


# dependant on time in mothers? NS
mtmod0 <- lm(viral_alpha_diversity ~ Timepoint_continuous+ DNA_CONC + Clean_reads, data = VLP_metadata[VLP_metadata$Type=="Mother",])
mtmod1  = lmer(viral_alpha_diversity ~ Timepoint_continuous+ DNA_CONC + Clean_reads + (1|NEXT_ID), REML = F, data = VLP_metadata[VLP_metadata$Type=="Mother",])
BIC(mtmod0, mtmod1)
exactLRT(mtmod1,mtmod0)
summary(mtmod1) #Timepoint_continuous 4.017e-04  6.795e-02 9.027e+01   0.006   0.9953   


# alpha-diversity plot virome
pdf('./04.PLOTS/Virome_diversity_plot.pdf', width=12/2.54, height=9/2.54)
ggplot(VLP_metadata, aes(Timepoint, viral_alpha_diversity, fill=Type)) + 
  labs (y="Shannon Diversity Index", x="") + 
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
dev.off()

# similar to maternal virome diversity at the later timepoints?
wilcox.test(VLP_metadata[VLP_metadata$Timepoint=='M12',]$viral_alpha_diversity, VLP_metadata[VLP_metadata$Type=='Mother',]$viral_alpha_diversity)

n_m <- length(VLP_metadata[VLP_metadata$Type=='Mother',]$viral_alpha_diversity)
n_i <- length(VLP_metadata[VLP_metadata$Timepoint=='M12',]$viral_alpha_diversity)
s_m <- sd(VLP_metadata[VLP_metadata$Type=='Mother',]$viral_alpha_diversity)
s_i <- sd(VLP_metadata[VLP_metadata$Timepoint=='M12',]$viral_alpha_diversity)
mean_m <- mean(VLP_metadata[VLP_metadata$Type=='Mother',]$viral_alpha_diversity)
mean_i <- mean(VLP_metadata[VLP_metadata$Timepoint=='M12',]$viral_alpha_diversity)
s_p = sqrt(((n_m - 1) * s_m^2 + (n_i - 1) * s_i^2) / (n_m + n_i - 2))
Cohens_d = (mean_m - mean_i) / s_p

# dependent on phenotypes in infants? 
btmod1_feeding  = lmer(viral_alpha_diversity ~ Age_days + infant_ever_never_breastfed + DNA_CONC + Clean_reads + (1|NEXT_ID), REML = F, data = VLP_metadata[VLP_metadata$Type=="Infant",])
summary(btmod1_feeding) #infant_ever_never_breastfed  8.786e-01  3.422e-01  1.941e+01   2.567   0.0187 *


pdf('./04.PLOTS/Virome_diversity_feeding_plot.pdf', width=12/2.54, height=9/2.54)

virome_feeding_plot <- melt(VLP_metadata[VLP_metadata$Type=='Infant' & !is.na(VLP_metadata$infant_ever_never_breastfed),c("Timepoint","viral_alpha_diversity", "infant_ever_never_breastfed", "Short_sample_ID")])
virome_feeding_plot$infant_ever_never_breastfed <- as.factor(virome_feeding_plot$infant_ever_never_breastfed)

ggplot(virome_feeding_plot, aes(Timepoint, value, fill=infant_ever_never_breastfed)) + 
  labs (y="Shannon Diversity Index", x="") + 
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
dev.off()

# alpha-diversity in microbiome

# dependent on type?
type_mod0_bac <- lm(bacterial_alpha_diversity ~ Type + DNA_CONC + Clean_reads, data = MGS_metadata)
type_mod1_bac  = lmer(bacterial_alpha_diversity ~ Type + DNA_CONC + Clean_reads + (1|NEXT_ID), REML = F, data = MGS_metadata)
BIC(type_mod0_bac, type_mod1_bac)
exactLRT(type_mod1_bac,type_mod0_bac)
summary(type_mod1_bac) # Type 1.767321e+00 9.933489e-02  76.76975 17.791543 5.861632e-29
as.data.frame(summary(type_mod1_bac)$coefficients)[,1:5]


# dependent on type & time?
mod0_bac <- lm(bacterial_alpha_diversity ~ Timepoint_continuous + DNA_CONC + Clean_reads, data = MGS_metadata)
mod1_bac  = lmer(bacterial_alpha_diversity ~ Type + Timepoint_continuous + DNA_CONC + Clean_reads + (1|NEXT_ID), REML = F, data = MGS_metadata)
BIC(mod0_bac, mod1_bac)
exactLRT(mod1_bac,mod0_bac)
summary(mod1_bac) #Timepoint_continuous 5.719e-02  7.041e-03  2.759e+02   8.122 1.54e-14 ***



# dependant on time in babies?
btmod0_bac <- lm(bacterial_alpha_diversity ~ Age_days + DNA_CONC + Clean_reads, data = MGS_metadata[MGS_metadata$Type=="Infant",])
btmod1_bac  = lmer(bacterial_alpha_diversity ~ Age_days + DNA_CONC + Clean_reads + (1|NEXT_ID), REML = F, data = MGS_metadata[MGS_metadata$Type=="Infant",])
BIC(btmod0_bac, btmod1_bac)
exactLRT(btmod1_bac,btmod0_bac)
summary(btmod1_bac) #Age_days     2.707e-03  3.521e-04  1.552e+02   7.689 1.59e-12 ***


# dependant on time in mothers? NS
mtmod0_bac <- lm(bacterial_alpha_diversity ~ Timepoint_continuous+ DNA_CONC + Clean_reads, data = MGS_metadata[MGS_metadata$Type=="Mother",])
mtmod1_bac  = lmer(bacterial_alpha_diversity ~ Timepoint_continuous+ DNA_CONC + Clean_reads + (1|NEXT_ID), REML = F, data = MGS_metadata[MGS_metadata$Type=="Mother",])
BIC(mtmod0_bac, mtmod1_bac)
exactLRT(mtmod1_bac,mtmod0_bac)
summary(mtmod1_bac) #Timepoint_continuous 1.125e-02  1.001e-02  1.204e+02   1.124   0.2634   


# alpha-diversity plot bacteriome
pdf('./04.PLOTS/Bacteriome_diversity_plot.pdf', width=12/2.54, height=9/2.54)
ggplot(MGS_metadata, aes(Timepoint, bacterial_alpha_diversity, fill=Type)) + 
  labs (y="Shannon Diversity Index", x="") + 
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
dev.off()


# dependent on phenotypes in infants? NS
btmod1_bac_a  = lmer(bacterial_alpha_diversity ~ Age_days + infant_ever_never_breastfed + DNA_CONC + Clean_reads + (1|NEXT_ID), REML = F, data = MGS_metadata[MGS_metadata$Type=="Infant",])
summary(btmod1_bac_a) #infant_ever_never_breastfed  1.569e-01  1.735e-01  2.834e+01   0.904    0.373  

pdf('./04.PLOTS/Bacteriome_diversity_feeding_plot.pdf', width=12/2.54, height=9/2.54)

bacteriome_feeding_plot <- melt(MGS_metadata[MGS_metadata$Type=='Infant' & !is.na(MGS_metadata$infant_ever_never_breastfed),c("Timepoint","bacterial_alpha_diversity", "infant_ever_never_breastfed", "Short_sample_ID")])
bacteriome_feeding_plot$infant_ever_never_breastfed <- as.factor(bacteriome_feeding_plot$infant_ever_never_breastfed)

ggplot(bacteriome_feeding_plot, aes(Timepoint, value, fill=infant_ever_never_breastfed)) + 
  labs (y="Shannon Diversity Index", x="") + 
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
dev.off()

# plot for virome and bacteriome together:
together_plot <- rbind(data.frame(unname(VLP_metadata[,c("Type", "viral_alpha_diversity")])),
                       data.frame(unname(MGS_metadata[,c("Type", "bacterial_alpha_diversity")])))
together_plot$Source <- substr(row.names(together_plot), 8,8)
together_plot[together_plot$Source=='V',]$Source <- 'Virome'
together_plot[together_plot$Source=='M',]$Source <- 'Bacteriome'
together_plot$Source <- factor(together_plot$Source, levels = c('Virome', 'Bacteriome'), ordered = T)
together_plot$significance_level <- '***'

pdf('./04.PLOTS/Virome_bacteriome_diversity_together.pdf', width=12/2.54, height=12/2.54)

ggplot(together_plot, aes(X1, X2, fill=X1)) + 
  labs (y="Shannon Diversity Index", x="") + 
  geom_violin(alpha=0.7) + 
  geom_boxplot(width = 0.1, outlier.shape = NA, alpha=0.9) + 
  geom_sina(colour='#303841', size=0.5) +
  theme_bw()+
  theme(axis.text=element_text(size=12), 
        axis.title=element_text(size=16,face="bold"),
        strip.text.x = element_text(size = 12),
        legend.position = "none",
        legend.title = element_text(size=12, face="bold")) +
  geom_signif(comparisons = list(c("Infant", "Mother")), 
              map_signif_level=TRUE) +
  facet_grid(~Source)
dev.off()

