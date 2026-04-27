setwd("~/Desktop/Projects_2021/NEXT_virome/09.DATA_ANALYSIS/")

#############################################################
# Here we explore which viral taxa shape the virome
#############################################################

#############################################################
# 0. Used files source
#############################################################

#############################################################
# 1. Functions
#############################################################
# 1. function is rigid to the name of the column containing sample ID:
# it has to be Sequencing_ID
# 2. check how it pulls the beta and p-value IF you decide to correct for
# anything else besides time
running_LMM <- function(phenos_list, # vector of phenotypes to test
                        metadata_w_phenos, # metadata containing these phenos as columns
                        RPKM_table, # filtered by prevalence table
                        pseudocount # calculated pseudocount, no default!
                        ) { 
  
  map_dfr(phenos_list, function(pheno){
    
    pheno_df <- metadata_w_phenos %>%
      filter(!is.na(!!sym(pheno)))
    
    smaller <- RPKM_table[,colnames(RPKM_table) %in% pheno_df$Sequencing_ID]
    smaller <- smaller[rowSums(smaller) > 0,]
    
    feature_set <- row.names(smaller)
    
    smaller <- smaller %>%
      rownames_to_column(var = "New_CID") %>%
      pivot_longer(!New_CID) %>%
      pivot_wider(names_from = New_CID) %>%
      left_join(pheno_df, by = c("name" = "Sequencing_ID"))
    
    map_dfr(feature_set, function(feature) {
      
      #F1 <- as.formula(paste0("log(", feature, "+", pseudocount/2, ")", " ~ ", pheno, " + Timepoint_new + (1|NEXT_ID)"))
      
      F1 <- as.formula(paste0("log(", feature, "+", pseudocount/2, ")", " ~ ", pheno, " + perc_aligned_cf + reads_lost_QC + sequencing_batch + sc_enrichment + Timepoint_new + (1|NEXT_ID)"))
      
      model <- lmer(F1, 
                    data = smaller,
                    REML = F)  
      
      summary(model)$coefficients %>%
        as.data.frame() %>%
        rownames_to_column(var = "rowname") %>%
        filter(!grepl("Timepoint", rowname) & # I
                 !grepl("Intercept", rowname) & # don't
                 !grepl("perc_aligned_cf", rowname) & # know
                 !grepl("reads_lost_QC", rowname) & # why
                 !grepl("sequencing_batch", rowname) & # lmer adds first factor value to the output........
                 !grepl("sc_enrichment", rowname)
                 ) %>%
        mutate(feature = feature,
               n_nzero_feature = sum(smaller %>% 
                                       filter(!!sym(feature) != 0) %>%
                                       pull(!!sym(feature)) > 0),
               n_non_NA_pheno = nrow(smaller)  )
      
    })
  })
}

#############################################################
# 1. Loading libraries
#############################################################
library(dplyr)
library(tidyverse)
library(lme4)
library(lmerTest)
library(MetBrewer)
library(vegan)
library(skimr)
library(ggplot2)
library(ggrepel)
library(patchwork)
#############################################################
# 2. Load Input Data
#############################################################
# phenotype selection:
phenos_selected <- readxl::read_xlsx("07.RESULTS/Adonis2_betadisper_VLP_virome.xlsx")

# metadatas:
smeta_w_phenos <- read.delim('06.CLEAN_DATA/02.FINAL/Chiliadal_meta_VLP_matched_v05_suppl_w_phenotypes.txt', sep='\t', header=T)
smeta_w_phenos <- smeta_w_phenos %>%
  mutate(Timepoint_new = factor(Timepoint_new, levels = c("M1", "M3", "M6", "M12"), ordered = T))

# abundance tables etc:
# VLP
VLP <- read.table("06.CLEAN_DATA/02.FINAL/VLP_only_RPKM_table_VLP_MGS_dec99ANI_ALL_CS_ab3kbp_1110_samples.txt", sep='\t', header=T)
VLP <- VLP[,colnames(VLP) %in% smeta_w_phenos$Sequencing_ID]
VLP <- VLP[rowSums(VLP) > 0,]

# transform:
VLP_RAB <- as.data.frame(t(as.data.frame(t(VLP)/colSums(VLP))))

# HLV
holovirome <- read.table("06.CLEAN_DATA/Intermediate/Holovirome_RPKM_1110samples_120997vOTUs.txt", sep='\t', header=T)
holovirome <- holovirome[,colnames(holovirome) %in% smeta_w_phenos$Universal_ID]
holovirome <- holovirome[rowSums(holovirome) > 0,]

# virus metadata
ETOF_vOTUr <- read.table('06.CLEAN_DATA/02.FINAL/Working_ETOF_120997vOTUr_ab3kbp_in_2200_VLP_MGS.txt', sep='\t', header=T)
ETOF_vOTUr <- ETOF_vOTUr %>%
  #filter(New_CID %in% row.names(VLP)) %>%
  separate(
    col    = tax_ictv_aai,
    into = c('Life', 'Realm', 'Kingdom', 'Phylum',  'Class', 'Order', 'Family',  'Genus',  'Species', 'Strain'),
    sep    = ";",
    fill   = "right",               
    remove = FALSE                  
  ) %>%
  separate(
    col    = Host_taxonomy,
    into = paste0('Host_', c('Domain', 'Phylum',  'Class', 'Order', 'Family',  'Genus')),
    sep    = ";",
    fill   = "right",               
    remove = FALSE                  
  )

# temperate:
temperates <- ETOF_vOTUr$New_CID[ETOF_vOTUr$lifestyle == "Temperate"]

# adding to the metadata
smeta_w_phenos$vir_richness_holo <- colSums(holovirome)[match(smeta_w_phenos$Universal_ID, colnames(holovirome))]
smeta_w_phenos$temperate_richness_holo <- colSums(holovirome[row.names(holovirome) %in% temperates,])[match(smeta_w_phenos$Universal_ID, colnames(holovirome))]

set.seed(444)

#############################################################
# 3.1 Analysis: vOTU level
#############################################################
min_rpkm <- VLP %>%
  rownames_to_column(var = "New_CID") %>%
  pivot_longer(!New_CID) %>%
  filter(value != 0) %>%
  pull(value) %>%
  min()

VLP_filtered_10 <- VLP[rowSums(VLP > 0) >= round(0.1 * ncol(VLP), 0), ]
VLP_filtered_10 <- VLP_filtered_10[,colSums(VLP_filtered_10) > 0]


vOTU_vs_pehnos <- running_LMM(phenos_list = phenos_selected %>%
                       filter(analysis == "shaping" & FDR < 0.05 & new_name!="mother_birthcardhealth_gravidity") %>%
                       pull(new_name),
                     metadata_w_phenos = smeta_w_phenos,
                     RPKM_table = VLP_filtered_10,
                     pseudocount = min_rpkm)



vOTU_vs_pehnos$FDR <- p.adjust(vOTU_vs_pehnos$`Pr(>|t|)`, 'BH')

writexl::write_xlsx(vOTU_vs_pehnos, '07.RESULTS/LMM_vOTUs_shaping_phenos.xlsx')

#############################################################
# 3.1.1 Analysis: visualization: feeding mode vOTU
#############################################################
feeding_vOTU_df <- vOTU_vs_pehnos %>%
  filter(rowname == "infant_ffq_feeding_mode_simpleexcl_FF") %>%
  left_join(ETOF_vOTUr %>% 
              select(New_CID, Host_Genus) %>%
              mutate(Host_Genus = gsub("g__", "", Host_Genus)), by = c("feature" = "New_CID")) %>%
  mutate(
    neg_log10_p = -log10(`Pr(>|t|)`),
    significant = FDR < 0.05,
    label = case_when(
      FDR < 0.05 ~ Host_Genus,
      TRUE ~ NA_character_
    )
  ) %>%
  group_by(Host_Genus) %>%
  mutate(
    babel = case_when(
      FDR < 0.05 & `Pr(>|t|)` == min(`Pr(>|t|)`[FDR < 0.05], na.rm = TRUE) ~ Host_Genus,
      TRUE ~ NA_character_
    )
  ) %>%
  ungroup()

feeding_vulcano <- ggplot(feeding_vOTU_df, aes(x = Estimate, y = neg_log10_p, color = label)) +
  ggrastr::rasterise(geom_point(alpha = 0.7, size = 1), dpi = 600) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey50") +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_text_repel(aes(label = babel), size = 2.5, max.overlaps = 10) +
  labs(x = "Effect size (log RPKM)", 
       y = expression(-log[10]*"(p-value)"),
       color = "Host genus") +
  scale_color_manual(values = c(met.brewer("Renoir")[sample(1:12)]),
                     na.value = "grey90") +
  annotate(geom = "text", x = 0.8, y = 13, label = "Exclusively\nformula-fed", size = 3) + 
  annotate(geom = "text", x = -0.8, y = 13, label = "Exclusively\nbreastfed", size = 3) +
  theme_minimal() + 
  theme(axis.text = element_text(size = 8),
        axis.title = element_text(size = 9),
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 8),
        legend.spacing.y = unit(0, 'cm'),
        legend.key.height = unit(0.5, "cm")) +
  guides(color = guide_legend(byrow = TRUE))

ggsave("05.PLOTS/07.Health_outcomes/VLP_feeding_vOTU_vulcano.pdf",
       feeding_vulcano, "pdf", width=10, height=8, units="cm", dpi = 300)


# stat for text:
feeding_vOTU_genus <- vOTU_vs_pehnos %>%
  filter(rowname == "infant_ffq_feeding_mode_simpleexcl_FF" & FDR < 0.05) %>%
  left_join(ETOF_vOTUr %>% 
              select(New_CID, Genus), by = c("feature" = "New_CID")) %>%
  group_by(Genus) %>%
  summarise(n = n())
#############################################################
# 3.1.2 Analysis: visualization: delivery mode vOTU
#############################################################
delivery_vOTU_df <- vOTU_vs_pehnos %>%
  filter(rowname == "birth_deliverybirthcard_mode_binaryVG") %>%
  left_join(ETOF_vOTUr %>% 
              select(New_CID, Host_Genus) %>%
              mutate(Host_Genus = gsub("g__", "", Host_Genus)), by = c("feature" = "New_CID")) %>%
  mutate(
    neg_log10_p = -log10(`Pr(>|t|)`),
    significant = FDR < 0.05,
    label = case_when(
      FDR < 0.05 ~ Host_Genus,
      TRUE ~ NA_character_
    )
  ) %>%
  group_by(Host_Genus) %>%
  mutate(
    babel = case_when(
      FDR < 0.05 & `Pr(>|t|)` == min(`Pr(>|t|)`[FDR < 0.05], na.rm = TRUE) ~ Host_Genus,
      TRUE ~ NA_character_
    )
  ) %>%
  ungroup()

delivery_vulcano <- ggplot(delivery_vOTU_df, aes(x = Estimate, y = neg_log10_p, color = label)) +
  ggrastr::rasterise(geom_point(alpha = 0.7, size = 1), dpi = 600) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey50") +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_text_repel(aes(label = babel), size = 2.5, max.overlaps = 10) +
  labs(x = "Effect size (log RPKM)", 
       y = expression(-log[10]*"(p-value)"),
       color = "Host genus") +
  scale_color_manual(values = c(met.brewer("Renoir")[sample(1:12)]),
                     na.value = "grey90") +
  annotate(geom = "text", x = 1.5, y = 11, label = "Vaginal delivery", size = 3) + 
  annotate(geom = "text", x = -1.5, y = 11, label = "C-section", size = 3) +
  theme_minimal() + 
  theme(axis.text = element_text(size = 8),
        axis.title = element_text(size = 9),
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 8),
        legend.spacing.y = unit(0, 'cm'),
        legend.key.height = unit(0.5, "cm")) +
  guides(color = guide_legend(byrow = TRUE))

ggsave("05.PLOTS/07.Health_outcomes/VLP_delivery_vOTU_vulcano.pdf",
       delivery_vulcano, "pdf", width=10, height=8, units="cm", dpi = 300)

# stat for text:
delivery_vOTU_host_genus <- vOTU_vs_pehnos %>%
  filter(rowname == "birth_deliverybirthcard_mode_binaryVG" & FDR < 0.05) %>%
  left_join(ETOF_vOTUr %>% 
              select(New_CID, Host_Genus), by = c("feature" = "New_CID")) %>%
  group_by(Host_Genus) %>%
  summarise(n = n()) %>%
  arrange(desc(n))
#############################################################
# 3.2 Analysis: defining viral taxa level to analyze
#############################################################

UNC_stat <- ETOF_vOTUr %>%
  filter(New_CID %in% row.names(VLP)) %>%
  select(New_CID, Species, Genus, Family) %>%
  pivot_longer(-New_CID, names_to = "level", values_to = "taxon") %>%
  group_by(level) %>%
  summarise(unclassified_frac = sum(taxon %in% c("Unclassified", "Unassigned")) / nrow(VLP)) %>%
  mutate(classified_frac = 1 - unclassified_frac)

#############################################################
# 3.3 Analysis: trying family again
#############################################################
by_family <- VLP %>%
  rownames_to_column(var = "New_CID") %>%
  left_join(ETOF_vOTUr) %>%
  group_by(Family) %>%
  summarise(across(all_of(colnames(VLP)), ~ sum(.x))) %>%
  filter( rowSums( across(all_of(colnames(VLP))) > 0) >= 0.1*ncol(VLP)) %>%
  filter(!Family %in% c("Unclassified", "Unassigned")) %>%
  column_to_rownames(var = "Family")

family_vs_phenos <- running_LMM(phenos_list = phenos_selected %>%
                                  filter(analysis == "shaping" & FDR < 0.05 & new_name!="mother_birthcardhealth_gravidity") %>%
                                  pull(new_name),
                                metadata_w_phenos = smeta_w_phenos,
                                RPKM_table = by_family,
                                pseudocount = min_rpkm)

family_vs_phenos$FDR <- p.adjust(family_vs_phenos$`Pr(>|t|)`, 'BH')

writexl::write_xlsx(family_vs_phenos, '07.RESULTS/LMM_virFamily_shaping_phenos.xlsx')

ETOF_vOTUr %>%
  filter(New_CID %in% row.names(VLP) & Family == "Suoliviridae") %>%
  group_by(Host_Genus) %>%
  summarise(n = n()) %>%
  arrange(desc(n))

#############################################################
# 3.3.1 Analysis: Visualization of Family association
#############################################################
plot_Suoliviridae_data <- by_family %>%
  rownames_to_column(var = "Family") %>%
  pivot_longer(!Family) %>%
  pivot_wider(names_from = Family) %>%
  select(name, all_of(family_vs_phenos$feature[family_vs_phenos$FDR < 0.05])) %>%
  left_join(smeta_w_phenos, by = c("name" = "Sequencing_ID")) %>%
  filter(!is.na(birth_deliverybirthcard_mode_binary))

# detection panel
p1 <- plot_Suoliviridae_data %>%
  group_by(Timepoint_new, birth_deliverybirthcard_mode_binary) %>%
  summarise(n_total = n(),
    n_present = sum(Suoliviridae > 0),
    .groups = "drop") %>%
  ggplot(aes(x = Timepoint_new, y = n_present, fill = birth_deliverybirthcard_mode_binary)) +
  geom_col(position = position_dodge(width = 0.9), alpha = 0.5, color = "black") +
  geom_text(aes(label = paste0(n_present, "/", n_total)),
    position = position_dodge(width = 0.9),
    vjust = -0.4,
    size = 2) +
  labs(y = "Presence", x = NULL, fill = "Delivery mode") +
  scale_fill_manual(values = met.brewer("Morgenstern")[c(2,5)]) +
  ylim(0, 180) +
  theme_minimal() +
  theme(axis.title = element_text(size = 9),
        axis.text = element_text(size = 8),
        legend.title = element_text(size = 9),
        legend.text = element_text(size=8),
        legend.key.size = unit(0.5, units = "cm"))

# abundance
p2 <- plot_Suoliviridae_data %>%
  filter(Suoliviridae > 0) %>%
  ggplot(aes(x = Timepoint_new, y = log(Suoliviridae), 
             fill = birth_deliverybirthcard_mode_binary,
             color = birth_deliverybirthcard_mode_binary)) +
  ggrastr::rasterise(geom_jitter(aes(color = birth_deliverybirthcard_mode_binary), 
                                 position = position_jitterdodge(), 
                                 size=0.1, alpha=1), dpi = 600) +
  geom_boxplot(alpha=0.4, color = "black") +
  labs(y = "Suoliviridae, log(RPKM > 0)", x = "Timepoint", fill = "Delivery mode", color = "Delivery mode") +
  scale_fill_manual(values = met.brewer("Morgenstern")[c(2,5)]) +
  scale_color_manual(values = met.brewer("Morgenstern")[c(2,5)]) +
  theme_minimal() +
  theme(axis.title = element_text(size = 9),
        axis.text = element_text(size = 8),
        legend.title = element_text(size = 9),
        legend.text = element_text(size=8),
        legend.key.size = unit(0.5, units = "cm"))

plot_Suoliviridae <- (p1 / p2) + plot_layout(guides = "collect", heights = c(3, 4)) & theme(legend.position = 'bottom')

ggsave("05.PLOTS/07.Health_outcomes/VLP_delivery_Suoliviridae.pdf",
       plot_Suoliviridae, "pdf", width=7, height=8, units="cm", dpi = 300)
#############################################################
# 3.4 Analysis: genome composition
#############################################################
by_genome_RPKM <- VLP %>%
  rownames_to_column(var = "New_CID") %>%
  left_join(ETOF_vOTUr %>% select(genome, New_CID)) %>%
  group_by(genome) %>%
  summarise(across(all_of(colnames(VLP)), ~ sum(.x))) %>%
  filter( rowSums( across(all_of(colnames(VLP))) > 0) >= 0.1*ncol(VLP)) %>%
  filter(genome != "Unclassified") %>%
  mutate(genome = paste0(genome, "_RPKM")) %>%
  column_to_rownames(var = "genome")

genome_vs_phenos <- running_LMM(phenos_list = phenos_selected %>%
                                  filter(analysis == "shaping" & FDR < 0.05 & new_name!="mother_birthcardhealth_gravidity") %>%
                                  pull(new_name),
                                metadata_w_phenos = smeta_w_phenos,
                                RPKM_table = by_genome_RPKM,
                                pseudocount = min_rpkm)


genome_vs_phenos$FDR <- p.adjust(genome_vs_phenos$`Pr(>|t|)`, 'BH')

writexl::write_xlsx(genome_vs_phenos, '07.RESULTS/LMM_virGenome_shaping_phenos.xlsx')
#############################################################
# 3.5 Analysis: host composition
#############################################################
by_host_RPKM <- VLP %>%
  rownames_to_column(var = "New_CID") %>%
  left_join(ETOF_vOTUr %>% select(New_CID, Host) %>% mutate(Host = ifelse(Host == "invertebrates, vertebrates", "inver_ver", Host))) %>%
  group_by(Host) %>%
  summarise(across(all_of(colnames(VLP)), ~ sum(.x))) %>%
  filter( rowSums(across(all_of(colnames(VLP))) > 0) >= 0.1*ncol(VLP) ) %>%
  filter(!Host %in% c("Unknown")) %>%
  mutate(Host = paste0(Host, "_RPKM")) %>%
  column_to_rownames(var = "Host")

host_vs_phenos <- running_LMM(phenos_list = phenos_selected %>%
                                filter(analysis == "shaping" & FDR < 0.05 & new_name!="mother_birthcardhealth_gravidity") %>%
                                pull(new_name),
                              metadata_w_phenos = smeta_w_phenos,
                              RPKM_table = by_host_RPKM,
                              pseudocount = min_rpkm)

host_vs_phenos$FDR <- p.adjust(host_vs_phenos$`Pr(>|t|)`, 'BH') 

writexl::write_xlsx(host_vs_phenos, '07.RESULTS/LMM_virHost_shaping_phenos.xlsx')
#############################################################
# 3.6 Analysis: viral sample metrics
#############################################################
virvars <- c("vir_richness_cf", "vir_diversity", "temperate_RAb", "lytic_RAb",
             "temperate_richness", "lytic_richness")

virvar_pseudocounts <- smeta_w_phenos %>%
  summarise(across(all_of(virvars), 
                   ~ min(.x[.x > 0], na.rm = TRUE))) %>%
  pivot_longer(cols = all_of(virvars))

virvars_vs_phenos <- map_dfr(phenos_selected %>%
                   filter(analysis == "shaping" & FDR < 0.05 & new_name!="mother_birthcardhealth_gravidity") %>%
                   pull(new_name), 
                 function(pheno){
                   
                   pheno_df <- smeta_w_phenos %>%
                     filter(!is.na(!!sym(pheno)))
                   
                   map_dfr(virvars, function(virvar) {
                     
                     pseudocount <- virvar_pseudocounts$value[virvar_pseudocounts$name == virvar]
                     
                     if (virvar == "vir_diversity") {
                       F1 <- as.formula(paste0(virvar, " ~ ", pheno, " + Timepoint_new + (1|NEXT_ID)"))
                     } else {
                       F1 <- as.formula(paste0("log(", virvar, " + ", pseudocount/2, ") ~ ", pheno,  " + perc_aligned_cf + reads_lost_QC + sequencing_batch + sc_enrichment +  Timepoint_new + (1|NEXT_ID)"))
                     }
                     
                     model <- lmer(F1, 
                                   data = pheno_df,
                                   REML = F)  
                     
                     summary(model)$coefficients %>%
                       as.data.frame() %>%
                       rownames_to_column(var = "rowname") %>%
                       filter(!grepl("Timepoint", rowname) & # I
                                !grepl("Intercept", rowname) & # don't
                                !grepl("perc_aligned_cf", rowname) & # know
                                !grepl("reads_lost_QC", rowname) & # why
                                !grepl("sequencing_batch", rowname) & # lmer adds first factor value to the output........
                       !grepl("sc_enrichment", rowname)
                                ) %>%
                       mutate(virvar = virvar,
                              n_nzero_virvar = sum(pheno_df %>% 
                                                      filter(!!sym(virvar) != 0) %>%
                                                      pull(!!sym(virvar)) > 0),
                              n_non_NA_pheno = nrow(pheno_df) )
                   })
                 })

virvars_vs_phenos <- virvars_vs_phenos %>%
  mutate(cat_virvar = ifelse(grepl("vir_", virvar), "Richness", "Lifestyle")) %>%
  group_by(cat_virvar) %>%
  mutate(FDR = p.adjust(`Pr(>|t|)`, 'BH')) %>%
  ungroup()

writexl::write_xlsx(virvars_vs_phenos, '07.RESULTS/LMM_virVars_shaping_phenos.xlsx')  

# sub-analysis: correcting for skunaviruses in richness of virulent phages:
virulent_no_skuna <- ETOF_vOTUr %>% 
  filter(Host_Genus != "g__Lactococcus" & lifestyle == "Virulent") %>% 
  pull(New_CID)

smeta_w_phenos$virulent_richness_no_skuna <- colSums(VLP[row.names(VLP) %in% virulent_no_skuna,] > 0)[match(smeta_w_phenos$Sequencing_ID, colnames(VLP))]
smeta_w_phenos$virulent_RAB <- colSums(VLP_RAB[row.names(VLP_RAB) %in% virulent_no_skuna,])[match(smeta_w_phenos$Sequencing_ID, colnames(VLP))]*100

sum(!is.na(smeta_w_phenos$infant_ffq_feeding_mode_simple))
summary(lmer(log(virulent_richness_no_skuna + 1/2) ~ infant_ffq_feeding_mode_simple + Timepoint_new + perc_aligned_cf + reads_lost_QC + sequencing_batch + (1|NEXT_ID), data = smeta_w_phenos, REML = F))

sum(smeta_w_phenos$virulent_RAB!=0 & !is.na(smeta_w_phenos$infant_ffq_feeding_mode_simple))
summary(lmer(log(virulent_RAB + 1/2) ~ infant_ffq_feeding_mode_simple + Timepoint_new + perc_aligned_cf + reads_lost_QC + sequencing_batch + (1|NEXT_ID), data = smeta_w_phenos, REML = F))
#############################################################
# 3.6.1 Analysis: visualization: temperate to lytic ratio & feeding
#############################################################

feeding_lytic_richness <- smeta_w_phenos %>%
  filter(!is.na(infant_ffq_feeding_mode_simple)) %>%
  filter(Timepoint_new != "M12") %>%
ggplot(aes(Timepoint_new, lytic_richness, 
           fill = infant_ffq_feeding_mode_simple)) +
  ggrastr::rasterise(geom_jitter(aes(color = infant_ffq_feeding_mode_simple), position = position_jitterdodge(), size=0.1, alpha=1), dpi = 300)+
  geom_boxplot(alpha=0.6, outlier.shape = NA) +
  theme_minimal() +
  theme(strip.background = element_rect(NA),
        legend.position = "bottom",
        axis.text = element_text(size=8),
        axis.title = element_text(size=9),
        legend.title = element_text(size=9),
        legend.text = element_text(size=8),
        legend.spacing.y = unit(0, 'cm'),
        legend.key.height = unit(0.4, "cm")) + 
  labs(y = "Virulent phages richness", x="Timepoint", fill = "Feeding mode", color = "Feeding mode") +
  scale_fill_manual(values = c("#6F0000", "#FF5200"), labels = c("Exclusive breastfeeding", "Exclusive formula-feeding")) +
  scale_color_manual(values = c("#6F0000", "#FF5200"), labels = c("Exclusive breastfeeding", "Exclusive formula-feeding")) +
  guides(fill = guide_legend(theme = theme(legend.title.position = "top", 
                                           legend.title = element_text(hjust = 0.5)), 
                             nrow = 2,
                             byrow = TRUE)) 

ggsave("05.PLOTS/07.Health_outcomes/VLP_feeding_lifestyle_richness.pdf",
       feeding_lytic_richness, "pdf", width=6, height=8, units="cm", dpi = 300)

#############################################################
# 3.6.2 Analysis: visualization: viral richness, feeding mode
#############################################################

VLP_richness <- smeta_w_phenos %>%
  filter(Timepoint_new != "M12") %>%
  filter(!is.na(infant_ffq_feeding_mode_simple)) %>%
  ggplot(aes(Timepoint_new, vir_richness_cf, 
             fill = infant_ffq_feeding_mode_simple)) +
  ggrastr::rasterise(geom_jitter(aes(color = infant_ffq_feeding_mode_simple), position = position_jitterdodge(), size=0.1, alpha=0.7), dpi = 300)+
  geom_boxplot(alpha=0.5, outlier.shape = NA) +
  labs(x = "Timepoint", y = "Viral richness, VLP", fill = "Feeding mode", color = "Feeding mode") +
  scale_fill_manual(values = c("#6F0000", "#FF5200"), labels = c("Exclusive breastfeeding", "Exclusive formula-feeding")) +
  scale_color_manual(values = c("#6F0000", "#FF5200"), labels = c("Exclusive breastfeeding", "Exclusive formula-feeding")) +
  theme_minimal() +
  theme(legend.position = "bottom",
        axis.text = element_text(size=8),
        axis.title = element_text(size=9),
        legend.title = element_text(size=9),
        legend.text = element_text(size=8),
        legend.spacing.y = unit(0, 'cm'),
        legend.key.height = unit(0.4, "cm")) +
  guides(fill = guide_legend(theme = theme(legend.title.position = "top", 
                                           legend.title = element_text(hjust = 0.5)), 
                             nrow = 2,
                             byrow = TRUE)) 

ggsave("05.PLOTS/07.Health_outcomes/VLP_feeding_richness.pdf",
       VLP_richness, "pdf", width=6, height=8, units="cm", dpi = 300)

#############################################################
# 3.6.3 Analysis: visualization: temperate RAb & delivery
#############################################################

temp_RAb_delivery <- smeta_w_phenos %>%
  filter(!is.na(birth_deliverybirthcard_mode_binary)) %>%
  ggplot(aes(Timepoint_new, temperate_RAb, 
             fill = birth_deliverybirthcard_mode_binary)) +
  ggrastr::rasterise(geom_jitter(aes(color = birth_deliverybirthcard_mode_binary), position = position_jitterdodge(), size=0.1, alpha=1), dpi = 300)+
  geom_boxplot(alpha=0.6, outlier.shape = NA) +
  theme_minimal() +
  theme(strip.background = element_rect(NA),
        legend.position = "bottom",
        axis.text = element_text(size=8),
        axis.title = element_text(size=9),
        legend.title = element_text(size=9),
        legend.text = element_text(size=8),
        legend.spacing.y = unit(0, 'cm'),
        legend.key.height = unit(0.4, "cm")) + 
  labs(y = "Temperate phages abundance", x="Timepoint", fill = "Delivery mode", color = "Delivery mode") +
  scale_fill_manual(values = met.brewer("Morgenstern")[c(2,5)], labels = c("C-section", "Vaginal delivery")) +
  scale_color_manual(values = met.brewer("Morgenstern")[c(2,5)], labels = c("C-section", "Vaginal delivery")) +
  guides(fill = guide_legend(theme = theme(legend.title.position = "top", 
                                           legend.title = element_text(hjust = 0.5)), 
                             nrow = 2,
                             byrow = TRUE)) 

ggsave("05.PLOTS/07.Health_outcomes/VLP_delivery_temperate_RAb.pdf",
       temp_RAb_delivery, "pdf", width=6, height=8, units="cm", dpi = 300)

#############################################################
# 3.6.4 Analysis: visualization: lytic RAb & feeding
#############################################################

lytic_RAb_feeding <- smeta_w_phenos %>%
  filter(!is.na(infant_ffq_feeding_mode_simple)) %>%
  filter(Timepoint_new != "M12") %>%
  ggplot(aes(Timepoint_new, lytic_RAb, 
             fill = infant_ffq_feeding_mode_simple)) +
  ggrastr::rasterise(geom_jitter(aes(color = infant_ffq_feeding_mode_simple), position = position_jitterdodge(), size=0.1, alpha=1), dpi = 300)+
  geom_boxplot(alpha=0.6, outlier.shape = NA) +
  theme_minimal() +
  theme(strip.background = element_rect(NA),
        legend.position = "bottom",
        axis.text = element_text(size=8),
        axis.title = element_text(size=9),
        legend.title = element_text(size=9),
        legend.text = element_text(size=8),
        legend.spacing.y = unit(0, 'cm'),
        legend.key.height = unit(0.4, "cm")) + 
  labs(y = "Virulent phages abundance", x="Timepoint", fill = "Feeding mode", color = "Feeding mode") +
  scale_fill_manual(values = c("#6F0000", "#FF5200"), labels = c("Exclusive breastfeeding", "Exclusive formula-feeding")) +
  scale_color_manual(values = c("#6F0000", "#FF5200"), labels = c("Exclusive breastfeeding", "Exclusive formula-feeding")) +
  guides(fill = guide_legend(theme = theme(legend.title.position = "top", 
                                           legend.title = element_text(hjust = 0.5)), 
                             nrow = 2,
                             byrow = TRUE)) 

ggsave("05.PLOTS/07.Health_outcomes/VLP_feeding_lytic_RAb.pdf",
       lytic_RAb_feeding, "pdf", width=6, height=8, units="cm", dpi = 300)

#############################################################
# 3.7 Analysis: holovirome viral richness
#############################################################
holovirvars <- c("vir_richness_holo", "temperate_richness_holo")

holovirvars_pseudocounts <- smeta_w_phenos %>%
  summarise(across(all_of(holovirvars), 
                   ~ min(.x[.x > 0], na.rm = TRUE))) %>%
  pivot_longer(cols = all_of(holovirvars))

holovirvars_vs_phenos <- map_dfr(phenos_selected %>%
                   filter(analysis == "shaping" & FDR < 0.05 & new_name!="mother_birthcardhealth_gravidity") %>%
                   pull(new_name), 
                 function(pheno){
                   
                   pheno_df <- smeta_w_phenos %>%
                     filter(!is.na(!!sym(pheno)))
                   
                   map_dfr(holovirvars, function(virvar) {
                     
                     pseudocount <- holovirvars_pseudocounts$value[holovirvars_pseudocounts$name == virvar]
                     
                     F1 <- as.formula(paste0("log(", virvar, " + ", pseudocount/2, ") ~ ", pheno, " + Timepoint_new + (1|NEXT_ID)"))
                     
                     model <- lmer(F1, 
                                   data = pheno_df,
                                   REML = F)  
                     
                     summary(model)$coefficients %>%
                       as.data.frame() %>%
                       rownames_to_column(var = "rowname") %>%
                       filter(!grepl("Timepoint", rowname) & !grepl("Intercept", rowname)) %>%
                       mutate(virvar = virvar,
                              n_nzero_virvar = sum(pheno_df %>% 
                                                     filter(!!sym(virvar) != 0) %>%
                                                     pull(!!sym(virvar)) > 0),
                              n_non_NA_pheno = nrow(pheno_df))
                   })
                 })

holovirvars_vs_phenos <- holovirvars_vs_phenos %>%
  mutate(cat_virvar = ifelse(virvar == "vir_richness_holo", "Richness", "Lifestyle")) %>%
  group_by(cat_virvar) %>%
  mutate(FDR = p.adjust(`Pr(>|t|)`, 'BH')) %>%
  ungroup()

writexl::write_xlsx(holovirvars_vs_phenos, '07.RESULTS/LMM_holovirVars_shaping_phenos.xlsx') 
#############################################################
# 3.7.1 Analysis: virusalization: holovirome richness vs feeding mode
#############################################################

feeding_holo_rich <- smeta_w_phenos %>%
  filter(!is.na(infant_ffq_feeding_mode_simple)) %>%
  filter(Timepoint_new != "M12") %>%
  ggplot(aes(Timepoint_new, vir_richness_holo, 
             fill = infant_ffq_feeding_mode_simple)) +
  ggrastr::rasterise(geom_jitter(aes(color = infant_ffq_feeding_mode_simple), position = position_jitterdodge(), size=0.1, alpha=1), dpi = 300)+
  geom_boxplot(alpha=0.6, outlier.shape = NA) +
  theme_minimal() +
  theme(strip.background = element_rect(NA),
        legend.position = "bottom",
        axis.text = element_text(size=8),
        axis.title = element_text(size=9),
        legend.title = element_text(size=9),
        legend.text = element_text(size=8),
        legend.spacing.y = unit(0, 'cm'),
        legend.key.height = unit(0.4, "cm")) + 
  labs(y = "Holovirome richness", x="Timepoint", fill = "Feeding mode", color = "Feeding mode") +
  scale_fill_manual(values = c("#6F0000", "#FF5200"), labels = c("Exclusive breastfeeding", "Exclusive formula-feeding")) +
  scale_color_manual(values = c("#6F0000", "#FF5200"), labels = c("Exclusive breastfeeding", "Exclusive formula-feeding")) +
  guides(fill = guide_legend(theme = theme(legend.title.position = "top", 
                                           legend.title = element_text(hjust = 0.5)), 
                             nrow = 2,
                             byrow = TRUE)) 

ggsave("05.PLOTS/07.Health_outcomes/Holovirome_feeding_richness.pdf",
       feeding_holo_rich, "pdf", width=6, height=8, units="cm", dpi = 300)

#############################################################
# 3.7.2 Analysis: virusalization: holovirome temperate richness vs feeding mode
#############################################################

# highlight the difference in holovirome?

feeding_holo_temp_rich <- smeta_w_phenos %>%
  filter(!is.na(infant_ffq_feeding_mode_simple)) %>%
  filter(Timepoint_new != "M12") %>%
  ggplot(aes(Timepoint_new, temperate_richness_holo, 
             fill = infant_ffq_feeding_mode_simple)) +
  ggrastr::rasterise(geom_jitter(aes(color = infant_ffq_feeding_mode_simple), position = position_jitterdodge(), size=0.1, alpha=1), dpi = 300)+
  geom_boxplot(alpha=0.6, outlier.shape = NA) +
  theme_minimal() +
  theme(strip.background = element_rect(NA),
        legend.position = "bottom",
        axis.text = element_text(size=8),
        axis.title = element_text(size=9),
        legend.title = element_text(size=9),
        legend.text = element_text(size=8),
        legend.spacing.y = unit(0, 'cm'),
        legend.key.height = unit(0.4, "cm")) + 
  labs(y = "Holovirome temperate phage richness", x="Timepoint", fill = "Feeding mode", color = "Feeding mode") +
  scale_fill_manual(values = c("#6F0000", "#FF5200"), labels = c("Exclusive breastfeeding", "Exclusive formula-feeding")) +
  scale_color_manual(values = c("#6F0000", "#FF5200"), labels = c("Exclusive breastfeeding", "Exclusive formula-feeding")) +
  guides(fill = guide_legend(theme = theme(legend.title.position = "top", 
                                           legend.title = element_text(hjust = 0.5)), 
                             nrow = 2,
                             byrow = TRUE)) 

ggsave("05.PLOTS/07.Health_outcomes/Holovirome_feeding_temperate_richness.pdf",
       feeding_holo_temp_rich, "pdf", width=5.5, height=8, units="cm", dpi = 300)

#############################################################
# 3.7.3 Analysis: virusalization: holovirome temperate richness vs delivery mode
#############################################################
delivery_holo_temp_rich <- smeta_w_phenos %>%
  filter(!is.na(birth_deliverybirthcard_mode_binary)) %>%
  ggplot(aes(Timepoint_new, temperate_richness_holo, 
             fill = birth_deliverybirthcard_mode_binary)) +
  ggrastr::rasterise(geom_jitter(aes(color = birth_deliverybirthcard_mode_binary), position = position_jitterdodge(), size=0.08, alpha=1), dpi = 300)+
  geom_boxplot(alpha=0.6, outlier.shape = NA) +
  theme_minimal() +
  theme(strip.background = element_rect(NA),
        legend.position = "bottom",
        axis.text = element_text(size=8),
        axis.title = element_text(size=9),
        legend.title = element_text(size=9),
        legend.text = element_text(size=8),
        legend.spacing.y = unit(0, 'cm'),
        legend.key.height = unit(0.4, "cm")) + 
  labs(y = "Holovirome temperate phage richness", x="Timepoint", fill = "Delivery mode", color = "Delivery mode") +
  scale_fill_manual(values = met.brewer("Morgenstern")[c(2,5)], labels = c("C-section", "Vaginal delivery")) +
  scale_color_manual(values = met.brewer("Morgenstern")[c(2,5)], labels = c("C-section", "Vaginal delivery")) +
  guides(fill = guide_legend(theme = theme(legend.title.position = "top", 
                                           legend.title = element_text(hjust = 0.5)), 
                             nrow = 2,
                             byrow = TRUE)) 

ggsave("05.PLOTS/07.Health_outcomes/Holovirome_delivery_temperate_richness.pdf",
       delivery_holo_temp_rich, "pdf", width=6, height=8, units="cm", dpi = 300)

#############################################################
# 3.8 Analysis: bacterial diversity between feeding modes
#############################################################
sum(!is.na(smeta_w_phenos$infant_ffq_feeding_mode_simple) & smeta_w_phenos$bacShannon > 0)
summary(lmer(bacShannon ~ infant_ffq_feeding_mode_simple + Timepoint_new + (1|NEXT_ID), data = smeta_w_phenos, REML = F))
