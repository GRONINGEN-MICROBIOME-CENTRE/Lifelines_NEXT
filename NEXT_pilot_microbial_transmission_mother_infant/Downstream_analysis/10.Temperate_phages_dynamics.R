setwd('~/Desktop/Projects_2022/NEXT_pilot_FUP/')

#############################################################
# Here we explore temperate phages of in gut microbiome and 
# virome
#############################################################

##############################
# Functions
##############################
mixed_models_taxa <- function(metadata, ID, CLR_transformed_data, pheno_list) {
  df <- metadata
  row.names(df) <- df[,ID]
  df<-merge(df, CLR_transformed_data, by='row.names')
  row.names(df) <- df$Row.names
  df$Row.names <- NULL
  
  Prevalent= c(colnames(CLR_transformed_data))
  #pheno_list= phenotypes
  
  Overall_result_phenos =tibble() 
  
  for (Bug in Prevalent){
    if (! Bug %in% colnames(df)){ next }
    #Prevalence = sum(as.numeric(as_vector(select(df, Bug)) > 0)) / dim(df)[1]
    # print (c(Bug, Prevalence))
    Bug2 = paste(c("`",Bug, "`"), collapse="")
    for ( pheno in pheno_list){
      pheno2 = paste(c("`",pheno, "`"), collapse="")
      df[is.na(df[colnames(df) == pheno]) == F, ID] -> To_keep
      df_pheno = filter(df, !!sym(ID) %in% To_keep )
      Model0 = as.formula(paste( c(Bug2,  " ~ DNA_CONC + Clean_reads + (1|Individual_ID)"), collapse="" )) 
      lmer(Model0, df_pheno) -> resultmodel0
      base_model=resultmodel0
      Model2 = as.formula(paste( c(Bug2,  " ~ DNA_CONC + Clean_reads +",pheno2, "+ (1|Individual_ID)"), collapse="" ))
      lmer(Model2, df_pheno, REML = F) -> resultmodel2
      M = "Mixed"
      as.data.frame(anova(resultmodel2, base_model))['resultmodel2','Pr(>Chisq)']->p_simp
      as.data.frame(summary(resultmodel2)$coefficients)[grep(pheno, row.names(as.data.frame(summary(resultmodel2)$coefficients))),] -> Summ_simple
      Summ_simple %>% rownames_to_column("Feature") %>% as_tibble() %>% mutate(P = p_simp, Model_choice = M, Bug =Bug, Pheno=pheno, Model="simple") -> temp_output
      rbind(Overall_result_phenos, temp_output) -> Overall_result_phenos
    }
  }
  
  p=as.data.frame(Overall_result_phenos)
  p$FDR<-p.adjust(p$P, method = "BH")
  
  return(p)
  
}

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

library(dplyr)
library(tibble)

library(EnhancedVolcano)
##############################
# Input data
##############################

VLP_metadata <- read.table("02.CLEAN_DATA/VLP_metadata_final_10_05_2023.txt", sep='\t', header=T, row.names = "Short_sample_ID")
VLP_metadata$Timepoint <- factor(VLP_metadata$Timepoint, levels = c("P7", "B",
                                                                    "M1", "M2", "M3",
                                                                    "M6", "M12"), ordered = T)
VLP_metadata$Type <- as.factor(VLP_metadata$Type)
VLP_metadata$Short_sample_ID <- row.names(VLP_metadata)
VLP_metadata$temperate_perc <- VLP_metadata$temperate_richness/VLP_metadata$viral_richness

MGS_metadata <- read.table("02.CLEAN_DATA/MGS_metadata_final_10_05_2023.txt", sep='\t', header=T, row.names = "Short_sample_ID")
MGS_metadata$Timepoint <- factor(MGS_metadata$Timepoint, levels = c("P3", "P7", "B",
                                                                    "M1", "M2", "M3",
                                                                    "M6", "M9", "M12"), ordered = T)
MGS_metadata$Type <- as.factor(MGS_metadata$Type)
MGS_metadata$Short_sample_ID <- row.names(MGS_metadata)


host_assignment_species <- read.csv('01.RAW_DATA/Host_assignment/MERGED_Host_prediction_to_genome_m90.csv', header=T)
# pre-cleaning:
host_assignment_species$Genus <- gsub('.*g__','', sapply(strsplit(host_assignment_species$Host.taxonomy, '\\;'), "[", 6))
host_assignment_species$Species <- gsub('.*s__','', sapply(strsplit(host_assignment_species$Host.taxonomy, '\\;'), "[", 7))
host_assignment_species <- host_assignment_species[host_assignment_species$Species!='',]
host_assignment_species <- host_assignment_species[order(host_assignment_species$Virus, host_assignment_species$Confidence.score ), ]
host_assignment_species_unique <- host_assignment_species[ !duplicated(host_assignment_species$Virus), ]
host_assignment_species_unique$Species_new <- gsub(' ', '_', host_assignment_species_unique$Species)
#### refining host assignment taxa based on most prevalent species and microbiome notation
host_assignment_species_unique$Species_new <- gsub('_[[:upper:]]', '',host_assignment_species_unique$Species_new)
# there is a typo in Metaphlan:
host_assignment_species_unique$Species_new <- sub('Hydrogeniiclostridium', 'Hydrogeniiclostidium', host_assignment_species_unique$Species_new)

contigs_metadata <- read.table('02.CLEAN_DATA/VLP_viral_contigs_metadata.txt', sep='\t', header=T)

RPKM_counts_VLP <- read.table('02.CLEAN_DATA/RPKM_counts_VLP.txt', sep='\t', header=T)
RPKM_counts_VLP <- RPKM_counts_VLP[row.names(RPKM_counts_VLP) %in% contigs_metadata[contigs_metadata$temperate==1,]$V1,]
RPKM_counts_VLP[RPKM_counts_VLP!=0] <- 1
temperate_binary <- RPKM_counts_VLP

host_assignment_species_temperate <- host_assignment_species_unique[host_assignment_species_unique$Virus %in% row.names(temperate_binary),]

### aggregating temperate phages by host:
temperate_by_BacSp <- data.frame(matrix(ncol = ncol(temperate_binary), nrow=length(unique(host_assignment_species_temperate$Species_new)) ) )
colnames(temperate_by_BacSp) <- colnames(temperate_binary)
row.names(temperate_by_BacSp) <- unique(host_assignment_species_temperate$Species_new)

for (i in colnames(temperate_by_BacSp) ) {
  
  for (j in row.names(temperate_by_BacSp) ) {
    
    temperate_by_BacSp[row.names(temperate_by_BacSp)==j,i] <- sum( temperate_binary[ c( host_assignment_species_temperate[ host_assignment_species_temperate$Species_new==j, ]$Virus ) , i] )   
    
  }
  
}

#considered_infants_feed has all infants that are considered for testing the differentially prevalent temperate phages
considered_infants_feed <- VLP_metadata[VLP_metadata$Type=='Infant' & !is.na(VLP_metadata$infant_ever_never_breastfed),]$Short_sample_ID
temperate_by_BacSp_infant <- temperate_by_BacSp[,considered_infants_feed]
temperate_by_BacSp_infant <- temperate_by_BacSp_infant[rowSums(temperate_by_BacSp_infant!=0)>round(0.05*ncol(temperate_by_BacSp_infant)),]
temperate_by_BacSp_infant <- as.data.frame(t(temperate_by_BacSp_infant ))

## microbiome
microbiome <- read.table('02.CLEAN_DATA/Microbiome_species_unfiltred.txt', sep='\t', header=T)
microbiome <- microbiome[,MGS_metadata[MGS_metadata$Type=='Infant' & !is.na(MGS_metadata$infant_ever_never_breastfed),]$Short_sample_ID_bact]
microbiome <- microbiome[rowSums(microbiome)!=0,]
microbiome <- as.data.frame(t(microbiome))
microbiome$Short_sample_ID_bact <- row.names(microbiome)
colnames(microbiome) <- gsub('.*s__', '', colnames(microbiome))

RPKM_counts_MGS <- read.table('02.CLEAN_DATA/RPKM_counts_MGS.txt', sep='\t', header = T)
RPKM_counts_MGS <- RPKM_counts_MGS[row.names(RPKM_counts_MGS) %in% contigs_metadata[contigs_metadata$temperate==1,]$V1,]
RPKM_counts_MGS[RPKM_counts_MGS!=0] <- 1
temperate_binary_MGS <- RPKM_counts_MGS

host_assignment_species_temperate_MGS <- host_assignment_species_unique[host_assignment_species_unique$Virus %in% row.names(temperate_binary_MGS),]
### aggregating temperate phages by host:
temperate_by_BacSp_MGS <- data.frame(matrix(ncol = ncol(temperate_binary_MGS), nrow=length(unique(host_assignment_species_temperate_MGS$Species_new)) ) )
colnames(temperate_by_BacSp_MGS) <- colnames(temperate_binary_MGS)
row.names(temperate_by_BacSp_MGS) <- unique(host_assignment_species_temperate_MGS$Species_new)

for (i in colnames(temperate_by_BacSp_MGS) ) {
  
  for (j in row.names(temperate_by_BacSp_MGS) ) {
    
    temperate_by_BacSp_MGS[row.names(temperate_by_BacSp_MGS)==j,i] <- sum( temperate_binary_MGS[ c( host_assignment_species_temperate_MGS[ host_assignment_species_temperate_MGS$Species_new==j, ]$Virus ) , i] )   
    
  }
  
}
considered_infants_feed_MGS <- MGS_metadata[MGS_metadata$Type=='Infant' & !is.na(MGS_metadata$infant_ever_never_breastfed),]$Short_sample_ID
temperate_by_BacSp_infant_MGS <- temperate_by_BacSp_MGS[,considered_infants_feed_MGS]
temperate_by_BacSp_infant_MGS <- temperate_by_BacSp_infant_MGS[rowSums(temperate_by_BacSp_infant_MGS!=0)>round(0.05*ncol(temperate_by_BacSp_infant_MGS)),]
temperate_by_BacSp_infant_MGS <- as.data.frame(t(temperate_by_BacSp_infant_MGS ))

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


# Differential Richness of active temperate phages
pdf('./04.PLOTS/Differential_temperate_richness.pdf', width=7/2.54, height=8.5/2.54)
ggplot(VLP_metadata, aes(Type, temperate_richness, fill=Type)) + 
  labs (y="Richness of active\ntemperate phages (log10)", x="") + 
  geom_violin() + 
  geom_boxplot(width = 0.1, outlier.shape = NA) + 
  geom_sina(colour='#303841', size=0.6) +
  geom_signif(comparisons = list(c("Infant", "Mother")), 
              map_signif_level=TRUE) +
  theme_bw()+
  theme(legend.position="none") + 
  scale_y_continuous(trans='log10') +
  theme(axis.text=element_text(size=10), 
        axis.title=element_text(size=12,face="bold"),
        strip.text.x = element_text(size = 12),
        legend.text = element_text(size=10),
        legend.title = element_text(size=12, face="bold")) +
  labs(colour = "Type")
dev.off()

# Richness of active temperate phages over time
pdf('./04.PLOTS/Dynamics_temperate_richness_infant.pdf', width=16/2.54, height=9/2.54)
ggplot(VLP_metadata, aes(Timepoint, temperate_richness, fill=Type)) + 
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
btmod1_feeding  = lmer(temperate_richness ~ Age_days + infant_ever_never_breastfed + DNA_CONC + Clean_reads + (1|NEXT_ID), REML = F, data = VLP_metadata[VLP_metadata$Type=="Infant",])
summary(btmod1_feeding) #infant_ever_never_breastfed  4.602e+01  1.213e+01  7.700e+01   3.792 0.000295 ***

btmod1_feeding  = lmer(temperate_RA ~ Age_days + infant_ever_never_breastfed + DNA_CONC + Clean_reads + (1|NEXT_ID), REML = F, data = VLP_metadata[VLP_metadata$Type=="Infant",])
summary(btmod1_feeding) #NS

Richness_Feeding <- melt(VLP_metadata[VLP_metadata$Type=='Infant' & !is.na(VLP_metadata$infant_ever_never_breastfed),c("Timepoint","temperate_richness", "infant_ever_never_breastfed", "Short_sample_ID")])
Richness_Feeding$infant_ever_never_breastfed <- as.factor(Richness_Feeding$infant_ever_never_breastfed)

pdf('./04.PLOTS/Temperate_richness_infant_feeding_mode.pdf', width=13/2.54, height=9/2.54)
ggplot(Richness_Feeding, aes(Timepoint, value, fill=infant_ever_never_breastfed)) + 
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
dev.off()


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
mtmod0 <- lm(temperate_RA ~ Timepoint_continuous+ DNA_CONC + Clean_reads, data = VLP_metadata[VLP_metadata$Type=="Mother",])
mtmod1  = lmer(temperate_RA ~ Timepoint_continuous+ DNA_CONC + Clean_reads + (1|NEXT_ID), REML = F, data = VLP_metadata[VLP_metadata$Type=="Mother",])
BIC(mtmod0, mtmod1)
exactLRT(mtmod1,mtmod0)
summary(mtmod1) #Timepoint_continuous -5.295e-01  3.108e-01  9.399e+01  -1.704   0.0917 . 

# still different at M12?
wilcox.test(VLP_metadata[VLP_metadata$Timepoint=='M12',]$temperate_RA, 
            VLP_metadata[VLP_metadata$Type=='Mother',]$temperate_RA)


# Differential relative abundance of active temperate phages
pdf('./04.PLOTS/Differential_temperate_RelAb.pdf', width=7/2.54, height=8.5/2.54)
ggplot(VLP_metadata, aes(Type, temperate_RA, fill=Type)) + 
  labs (y="Relative abundance of\nactive temperate phages", x="") + 
  geom_violin() + 
  geom_boxplot(width = 0.05, outlier.shape = NA) + 
  geom_sina(colour='#303841', size=0.6) +
  geom_signif(comparisons = list(c("Infant", "Mother")), 
              map_signif_level=TRUE) +
  theme_bw()+
  theme(legend.position="none") + 
  theme(axis.text=element_text(size=10), 
        axis.title=element_text(size=12,face="bold"),
        strip.text.x = element_text(size = 12),
        legend.text = element_text(size=10),
        legend.title = element_text(size=12, face="bold")) +
  labs(colour = "Type")
dev.off()

# Growth of relative abundance of active temperate phages
pdf('./04.PLOTS/Dynamics_of_temperate_RelAb.pdf', width=12/2.54, height=8.5/2.54)
ggplot(VLP_metadata, aes(Timepoint, temperate_RA, fill=Type)) + 
  labs (y="Relative abundance of\nactive temperate phages", x="") + 
  #geom_violin() + 
  geom_boxplot(width=0.5,outlier.shape = NA, alpha=0.5) + 
  geom_sina(aes(color=Type), size=0.6) +
  facet_grid(~Type, scales = 'free_x') +
  theme_bw()+
  theme(legend.position="none") + 
  theme(axis.text=element_text(size=10), 
        axis.title=element_text(size=12,face="bold"),
        strip.text.x = element_text(size = 12),
        legend.text = element_text(size=10),
        legend.title = element_text(size=12, face="bold")) +
  labs(colour = "Type")
dev.off()

#### richness of temperate phages depending on the feeding mode
diff_ab_temperate <- mixed_models_taxa(VLP_metadata[considered_infants_feed,], "Short_sample_ID", temperate_by_BacSp_infant, c('infant_ever_never_breastfed'))
diff_ab_temperate$log2FC <- NA

# necessary to calculate log2FC in cases when a bug is not present in one of the groups:
considered_infants_metadata <- VLP_metadata[considered_infants_feed,]
ever_fed <- considered_infants_metadata[considered_infants_metadata$infant_ever_never_breastfed=='ever_BF',]$Short_sample_ID
never_fed <- considered_infants_metadata[considered_infants_metadata$infant_ever_never_breastfed=='never_BF',]$Short_sample_ID

min_ever <- min(colMeans(temperate_by_BacSp_infant[ever_fed,])[colMeans(temperate_by_BacSp_infant[ever_fed,]) > 0])
min_never <- min(colMeans(temperate_by_BacSp_infant[never_fed,])[colMeans(temperate_by_BacSp_infant[never_fed,]) > 0])

for (i in diff_ab_temperate$Bug) {
  
  A <- mean( temperate_by_BacSp_infant[ ever_fed,  i ] ) 
  B <- mean( temperate_by_BacSp_infant[ never_fed,  i ] ) 
  
  if ( abs((log2(A) - log2(B)))==Inf ) {
    if (A==0) {
      A <- min_ever
    } else {
      B <- min_never
    }
    
  }
  diff_ab_temperate[diff_ab_temperate$Bug==i,]$log2FC <-  log2(B) - log2(A) 
}

row.names(diff_ab_temperate) <- diff_ab_temperate$Bug

pdf('./04.PLOTS/Differentially_prevalent_temperate_feeding.pdf', width=20/2.54, height=15/2.54)
EnhancedVolcano(diff_ab_temperate,
                lab = rownames(diff_ab_temperate),
                x = 'log2FC',
                y = 'FDR', 
                ylim = c(0, -log10(1e-2)),
                xlim = c(-6,6),
                title = 'Differentially prevalent temperate vOTUs\naggregated by host in infants\nexclusively formula fed and mixed fed',
                pCutoff = 0.05,
                FCcutoff = 1.5,
                drawConnectors = TRUE,
                widthConnectors = 0.75,
                boxedLabels = TRUE) + 
  ggplot2::labs(y=expression(-Log["10"]~FDR))

dev.off()

# is it mirrored by host abundance?
Bacteroides_fragilis <- merge(MGS_metadata, microbiome[,c("Bacteroides_fragilis", "Short_sample_ID_bact")], by='Short_sample_ID_bact')

B_frag_mod0 <- lm(Bacteroides_fragilis ~ infant_ever_never_breastfed + DNA_CONC + Clean_reads, data = Bacteroides_fragilis)
B_frag_mod1  = lmer(Bacteroides_fragilis ~ infant_ever_never_breastfed + DNA_CONC + Clean_reads + (1|NEXT_ID), REML = F, data = Bacteroides_fragilis)
BIC(B_frag_mod0, B_frag_mod1)
exactLRT(B_frag_mod1,B_frag_mod0)
summary(B_frag_mod1) #infant_ever_never_breastfed 1.032e+00  1.761e+00 2.832e+01   0.586    0.562

### is richness of temperates is still different when prophages are taken into account?
btmod1_mgs_feeding  = lmer(temperate_richness_MGS ~ Age_days + infant_ever_never_breastfed + DNA_CONC + Clean_reads + (1|NEXT_ID), REML = F, data = MGS_metadata[MGS_metadata$Type=="Infant",])
summary(btmod1_mgs_feeding) #infant_ever_never_breastfed  5.414e+01  2.454e+01 2.925e+01   2.206   0.0354 *

pdf('./04.PLOTS/Temperate_richness_MGS_infant_feeding_mode.pdf', width=13/2.54, height=9/2.54)
ggplot(MGS_metadata[considered_infants_feed_MGS,], aes(Timepoint, temperate_richness_MGS, fill=infant_ever_never_breastfed)) + 
  labs (y="Richness of temperate phages (log10)", x="") + 
  geom_boxplot(outlier.shape = NA, alpha=0.5) + 
  geom_sina(colour='#303841', size=0.6) +
  theme_bw()+
  scale_y_continuous(trans='log10') +
  theme(axis.text=element_text(size=12), 
        axis.title=element_text(size=12,face="bold"),
        strip.text.x = element_text(size = 12),
        legend.text = element_text(size=10),
        legend.title = element_text(size=12, face="bold")) +
  scale_fill_manual(name = "Feeding mode", 
                    labels=c('Ever breastfed', 'Never breastfed'),
                    values=c("#FAE3D9", "#BBDED6"))
dev.off()


#### richness of temperate phages depending on the feeding mode in MGS
diff_ab_temperate_MGS <- mixed_models_taxa(MGS_metadata[considered_infants_feed_MGS,], "Short_sample_ID", temperate_by_BacSp_infant_MGS, c('infant_ever_never_breastfed'))

###### OUTPUT #####
write.table(diff_ab_temperate, '03a.RESULTS/Differentially_prevalent_aggregated_temperates_Mixed_Models_Feeding.txt', sep='\t', quote=F, row.names=F)
write.table(diff_ab_temperate_MGS, '03a.RESULTS/Differentially_prevalent_aggregated_temperates_MGS_Mixed_Models_Feeding.txt', sep='\t', quote=F, row.names=F)



