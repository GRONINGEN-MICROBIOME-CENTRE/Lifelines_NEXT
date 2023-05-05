setwd('~/Desktop/Projects_2022/NEXT_pilot_FUP/')

#############################################################
# Here we explore differential abundance of gut bacterial
# species and viral vOTUs between maternal and infant gut 
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
      df$Short_sample_ID_bact[is.na(df[colnames(df) == pheno]) == F] -> To_keep
      df_pheno = filter(df, Short_sample_ID_bact %in% To_keep )
      Model0 = as.formula(paste( c(Bug2,  " ~ DNA_CONC + Clean_reads + (1|Individual_ID)"), collapse="" )) 
      lmer(Model0, df_pheno) -> resultmodel0
      base_model=resultmodel0
      Model2 = as.formula(paste( c(Bug2,  " ~ DNA_CONC + Clean_reads +",pheno2, "+ (1|Individual_ID)"), collapse="" ))
      lmer(Model2, df_pheno, REML = F) -> resultmodel2
      M = "Mixed"
      as.data.frame(anova(resultmodel2, base_model))[2,8]->p_simp
      as.data.frame(summary(resultmodel2)$coefficients)[3,1:5] -> Summ_simple
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
library(dplyr)
library(lme4)
library(lmerTest)
library(tibble)

library(EnhancedVolcano)
##############################
# Input data
##############################
MGS_metadata <- read.table('02.CLEAN_DATA/MGS_metadata_with_phenos.txt', sep='\t', header=T)
MGS_metadata$Type <- factor(MGS_metadata$Type, levels=c('Infant', 'Mother'), ordered = T)

microbiome <- read.table('02.CLEAN_DATA/Bacterial_spp_filtered_CLR_transformed.txt', sep='\t', header=T)
microbiome_filt_raw <- read.table('02.CLEAN_DATA/Bacterial_spp_filtered.txt', sep='\t', header=T)

VLP_metadata <- read.table('02.CLEAN_DATA/VLP_metadata_with_phenos.txt', sep='\t', header=T)
# it will take some time to load vOTUs, because it has vOTUs as columns
virome_vOTUs <- read.table('02.CLEAN_DATA/Viral_vOTUs_filtered_CLR_transformed.txt', sep='\t', header=T)
virome_vOTUs_filt_raw <- read.table('02.CLEAN_DATA/Viral_vOTUs_filtered.txt', sep='\t', header=T)

contigs_metadata <- read.table('02.CLEAN_DATA/VLP_viral_contigs_metadata.txt', sep='\t', header=T)
colnames(contigs_metadata)[1] <- 'Bug'

virome_aggr <- read.table('02.CLEAN_DATA/RPKM_counts_aggregated_bacterial_host_filtered_CLR_transformed.txt', sep='\t', header=T)
virome_aggr_filt_raw <- read.table('02.CLEAN_DATA/RPKM_counts_aggregated_bacterial_host_filt.txt', sep='\t', header=T)
##############################
# ANALYSIS
##############################

diff_ab_bacteria <- mixed_models_taxa(MGS_metadata, "Short_sample_ID_bact", microbiome, c('Type'))
diff_ab_bacteria$log2FC <- NA

# necessary to calculate log2FC in cases when a bug is not present in one of the groups:
min_mother <- min(colMeans(microbiome_filt_raw[MGS_metadata[MGS_metadata$Type=='Mother',]$Short_sample_ID_bact,])[colMeans(microbiome_filt_raw[MGS_metadata[MGS_metadata$Type=='Mother',]$Short_sample_ID_bact,]) > 0])
min_infant <- min(colMeans(microbiome_filt_raw[MGS_metadata[MGS_metadata$Type=='Infant',]$Short_sample_ID_bact,])[colMeans(microbiome_filt_raw[MGS_metadata[MGS_metadata$Type=='Infant',]$Short_sample_ID_bact,]) > 0])
for (i in diff_ab_bacteria$Bug) {
  
  A <- mean( microbiome_filt_raw[ MGS_metadata[MGS_metadata$Type=='Mother',]$Short_sample_ID_bact,  i ] ) 
  B <- mean( microbiome_filt_raw[ MGS_metadata[MGS_metadata$Type=='Infant',]$Short_sample_ID_bact,  i ] ) 
  
  if ( abs((log2(A) - log2(B)))==Inf ) {
    if (A==0) {
      A <- min_mother
    } else {
      B <- min_infant
    }
    
  }
  diff_ab_bacteria[diff_ab_bacteria$Bug==i,]$log2FC <-  log2(B) - log2(A) 
}
diff_ab_bacteria$Species <- gsub('.*s__','',diff_ab_bacteria$Bug)
row.names(diff_ab_bacteria) <- diff_ab_bacteria$Species

pdf('./04.PLOTS/Differentially_ab_BacSp_FC5_10e-10.pdf', width=22/2.54, height=17/2.54)
EnhancedVolcano(diff_ab_bacteria,
                lab = rownames(diff_ab_bacteria),
                x = 'log2FC',
                y = 'FDR', 
                selectLab = c('Clostridium_sp_Marseille_P3244','Clostridium_fessum','Clostridium_sp_C5_48',
                              'Clostridium_neonatale','Clostridium_paraputrificum',
                              'Bifidobacterium_scardovii','Bifidobacterium_breve','Bifidobacterium_dentium',
                              'Cutibacterium_acnes','Escherichia_coli'),
                title = 'Differentially abundant bacterial species in mothers and infants',
                pCutoff = 0.05,
                FCcutoff = 2,
                drawConnectors = TRUE,
                widthConnectors = 0.75,
                boxedLabels = TRUE) + 
  ggplot2::labs(y=expression(-Log["10"]~FDR))
                
dev.off()

# color by phylum?

##### vOTUs
diff_ab_vOTUs <- mixed_models_taxa(VLP_metadata, "Short_sample_ID", virome_vOTUs , c('Type'))
diff_ab_vOTUs$log2FC <- NA

# necessary to calculate log2FC in cases when a bug is not present in one of the groups:
min_mother <- min(colMeans(virome_vOTUs_filt_raw[VLP_metadata[VLP_metadata$Type=='Mother',]$Short_sample_ID,])[colMeans(virome_vOTUs_filt_raw[VLP_metadata[VLP_metadata$Type=='Mother',]$Short_sample_ID,]) > 0])
min_infant <- min(colMeans(virome_vOTUs_filt_raw[VLP_metadata[VLP_metadata$Type=='Infant',]$Short_sample_ID,])[colMeans(virome_vOTUs_filt_raw[VLP_metadata[VLP_metadata$Type=='Infant',]$Short_sample_ID,]) > 0])
for (i in diff_ab_vOTUs$Bug) {
  
  A <- mean( virome_vOTUs_filt_raw[ VLP_metadata[VLP_metadata$Type=='Mother',]$Short_sample_ID,  i ] ) 
  B <- mean( virome_vOTUs_filt_raw[ VLP_metadata[VLP_metadata$Type=='Infant',]$Short_sample_ID,  i ] ) 
  
  if ( abs((log2(A) - log2(B)))==Inf ) {
    if (A==0) {
      A <- min_mother
    } else {
      B <- min_infant
    }
    
  }
  diff_ab_vOTUs[diff_ab_vOTUs$Bug==i,]$log2FC <-  log2(B) - log2(A) 
}

# what are those differentially_abundant_vOTUs?
diff_ab_vOTUs <- merge(diff_ab_vOTUs, contigs_metadata, by='Bug')

# some descriptive stat:
significant_vOTUs <- diff_ab_vOTUs[diff_ab_vOTUs$FDR < 0.05,] #2,098 diff. abundant vOTUs
length(significant_vOTUs[!is.na(significant_vOTUs$class),]$class) # Caudoviricetes
sum(significant_vOTUs$temperate==1) # total number of temperate phaes among differentially abundant vOTUs
sum(significant_vOTUs[significant_vOTUs$log2FC<0,]$temperate==1) # how many enriched in mothers

#### contigID_easy with length, lifestyle and k-mer coverage (latter is added for uniqness)
diff_ab_vOTUs$cov <- round(as.numeric(sapply(strsplit(diff_ab_vOTUs$Bug, '\\_'), "[", 10)), 3)
diff_ab_vOTUs$ContigID_easy <- paste0('L', diff_ab_vOTUs$length, '_LS', diff_ab_vOTUs$temperate, '_cov_', diff_ab_vOTUs$cov)
row.names(diff_ab_vOTUs) <- diff_ab_vOTUs$ContigID_easy



pdf('./04.PLOTS/Differentially_ab_vOTUs_FC5_10e-3.pdf', width=22/2.54, height=17/2.54)
EnhancedVolcano(diff_ab_vOTUs,
                lab = rownames(diff_ab_vOTUs),
                x = 'log2FC',
                y = 'FDR', 
                selectLab = c(diff_ab_vOTUs[diff_ab_vOTUs$checkv_quality %in% c('Complete', 'High-quality') & diff_ab_vOTUs$FDR < 0.05,]$ContigID_easy),
                title = 'Differentially abundant vOTUs in mothers and infants',
                pCutoff = 0.05,
                FCcutoff = 2,
                drawConnectors = TRUE,
                widthConnectors = 0.75,
                boxedLabels = TRUE) + 
  ggplot2::labs(y=expression(-Log["10"]~FDR))

dev.off()

### viruses aggregated at the species-level of the predicted host
diff_ab_aggr_vir <- mixed_models_taxa(VLP_metadata, "Short_sample_ID", virome_aggr, c('Type'))
diff_ab_aggr_vir$log2FC <- NA

# necessary to calculate log2FC in cases when a bug is not present in one of the groups:
min_mother <- min(colMeans(virome_aggr_filt_raw[VLP_metadata[VLP_metadata$Type=='Mother',]$Short_sample_ID,])[colMeans(virome_aggr_filt_raw[VLP_metadata[VLP_metadata$Type=='Mother',]$Short_sample_ID,]) > 0])
min_infant <- min(colMeans(virome_aggr_filt_raw[VLP_metadata[VLP_metadata$Type=='Infant',]$Short_sample_ID,])[colMeans(virome_aggr_filt_raw[VLP_metadata[VLP_metadata$Type=='Infant',]$Short_sample_ID,]) > 0])
for (i in diff_ab_aggr_vir$Bug) {
  
  A <- mean( virome_aggr_filt_raw[ VLP_metadata[VLP_metadata$Type=='Mother',]$Short_sample_ID,  i ] ) 
  B <- mean( virome_aggr_filt_raw[ VLP_metadata[VLP_metadata$Type=='Infant',]$Short_sample_ID,  i ] ) 
  
  if ( abs((log2(A) - log2(B)))==Inf ) {
    if (A==0) {
      A <- min_mother
    } else {
      B <- min_infant
    }
    
  }
  diff_ab_aggr_vir[diff_ab_aggr_vir$Bug==i,]$log2FC <-  log2(B) - log2(A) 
}

row.names(diff_ab_aggr_vir) <- diff_ab_aggr_vir$Bug
significant_aggr_vir <- diff_ab_aggr_vir[diff_ab_aggr_vir$FDR<0.05,]

pdf('./04.PLOTS/Differentially_ab_host_agr_viruses_FC2_FDR.pdf', width=22/2.54, height=17/2.54)
EnhancedVolcano(diff_ab_aggr_vir,
                lab = rownames(diff_ab_aggr_vir),
                x = 'log2FC',
                y = 'FDR', 
                selectLab = c("Veillonella_rogosae","Veillonella_ratti","Bifidobacterium_scaligerum",
                              "Bifidobacterium_scardovii", "Bifidobacterium_dentium","Bifidobacterium_breve",
                              "Clostridium_neonatale","Escherichia_coli", "Alistipes_putredinis","Alistipes_ihumii", 
                              "Ruminiclostridium_sp900539195"),
                title = 'Differentially abundant viruses aggregated by host\nin mothers and infants',
                pCutoff = 0.05,
                FCcutoff = 2,
                drawConnectors = TRUE,
                widthConnectors = 0.75,
                boxedLabels = TRUE) + 
  ggplot2::labs(y=expression(-Log["10"]~FDR))

dev.off()


##############################
# OUTPUT
##############################
write.table(diff_ab_bacteria, "03a.RESULTS/Differentially_abundant_bacteria_CLR_Mixed_Models_Type.txt", sep="\t", row.names=F, quote = F)
write.table(diff_ab_vOTUs, "03a.RESULTS/Differentially_vOTUs_CLR_Mixed_Models_Type.txt", sep="\t", row.names=F, quote = F)
write.table(diff_ab_aggr_vir, "03a.RESULTS/Differentially_abundant_aggregated_viruses_CLR_Mixed_Models_Type.txt", sep="\t", row.names=F, quote = F)
