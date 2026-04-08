setwd("~/Desktop/Projects_2021/NEXT_virome/09.DATA_ANALYSIS/")

#############################################################
# vOTU sharing
#############################################################

#############################################################
# 0. Used files source
#############################################################

#############################################################
# 1. Functions
#############################################################

#############################################################
# 1. Loading libraries
#############################################################
library(dplyr)
library(tidyverse)
library(ggplot2)
library(lme4)
library(lmerTest)
library(Matrix)
library(mediation)
#############################################################
# 2. Load Input Data
#############################################################
# metadatas:
smeta <- read.delim('06.CLEAN_DATA/02.FINAL/Chiliadal_meta_VLP_matched_v05_suppl_w_virmetrics.txt', sep='\t', header=T)

smeta <- smeta %>%
  mutate(secpreg = grepl("P2", Family_structure)) %>%
  mutate(FAMILYupd = if_else(secpreg, paste0(FAMILY, "_P2"), FAMILY)) %>% # treating 2nd pregnancy as a separate family:
  mutate(Timepoint_new = factor(Timepoint_new, levels=c("M1", "M3", "M6", "M12", "Mother"), ordered = T))

long_smeta <- read.delim("06.CLEAN_DATA/02.FINAL/Chiliadal_meta_VLP_MGS_matched_v05_suppl_w_virmetrics.txt", sep='\t', header=T)
long_smeta <- long_smeta %>%
  filter(seq_type == "MGS")

smeta$corresp_MGS <- long_smeta$Sequencing_ID[match(smeta$Universal_ID, long_smeta$Universal_ID)]

# abundance tables etc:
VLP <- read.table('06.CLEAN_DATA/02.FINAL/VLP_only_RPKM_table_VLP_MGS_dec99ANI_ALL_CS_ab3kbp_1110_samples.txt', sep='\t', header=T)

ETOF_vOTUr <- read.table('06.CLEAN_DATA/02.FINAL/Working_ETOF_120997vOTUr_ab3kbp_in_2200_VLP_MGS.txt', sep='\t', header = T)
ETOF_vOTUr <- ETOF_vOTUr %>%
  mutate(Host = ifelse(Host == "invertebrates, vertebrates", "invertebrates_vertebrates", Host))

metaphlan <- read.table('06.CLEAN_DATA/02.FINAL/MGS_Chiliadal_metaphlan_full_taxonomy_ver_01_07102025.txt', sep='\t', header=T)
metaphlan <- metaphlan[grep('g__', row.names(metaphlan)),]
metaphlan <- metaphlan[!grepl('s__', row.names(metaphlan)),]
metaphlan <- metaphlan[rowSums(metaphlan) > 0,]
#############################################################
# 3.0 Analysis: calculating vOTU sharing
#############################################################
binary_matrix <- as.matrix(VLP > 0) + 0 
binary_matrix_t <- t(binary_matrix)

# matrix multiplication
shared_species_matrix <- binary_matrix_t %*% t(binary_matrix_t)
#rm(binary_matrix, binary_matrix_t)

upper_tri_idx <- which(upper.tri(shared_species_matrix, diag = FALSE), arr.ind = TRUE)

sharing <- data.frame(
  sample_1 = rownames(shared_species_matrix)[upper_tri_idx[, 1]],
  sample_2 = colnames(shared_species_matrix)[upper_tri_idx[, 2]],
  N_shared = shared_species_matrix[upper_tri_idx]) %>%
  left_join(smeta %>% select(Sequencing_ID,
                             Universal_ID,
                             Modified_NEXT_ID_without_preg_number,
                             FAMILY, 
                             FAMILYupd,
                             Type, 
                             seq_type,
                             Timepoint_new), by = c("sample_1" = "Sequencing_ID" )) %>%
  left_join(smeta %>% select(Sequencing_ID,
                             Universal_ID,
                             Modified_NEXT_ID_without_preg_number,
                             FAMILY,
                             FAMILYupd,
                             Type, seq_type,
                             Timepoint_new), by = c("sample_2" = "Sequencing_ID"), suffix = c("_1", "_2")) %>%
  select(sample_1, sample_2, N_shared, sort(colnames(.))) %>%
  # twins
  mutate(twins = ifelse(FAMILYupd_1 == FAMILYupd_2 &
                          Modified_NEXT_ID_without_preg_number_1 != Modified_NEXT_ID_without_preg_number_2 &
                          Type_1 == "K" &
                          Type_2 == "K", T, F)) %>%
  # siblings
  mutate(siblings = ifelse(FAMILY_1 == FAMILY_2 &
                             FAMILYupd_1 != FAMILYupd_2 &
                             Modified_NEXT_ID_without_preg_number_1 != Modified_NEXT_ID_without_preg_number_2 &
                             Type_1 == "K" &
                             Type_2 == "K", T, F)) %>%
  # mother-infant related
  mutate(same_dyad = ifelse(FAMILYupd_1 == FAMILYupd_2 &
                              Modified_NEXT_ID_without_preg_number_1 != Modified_NEXT_ID_without_preg_number_2 &
                              Type_1 != Type_2, T, F)) %>%
  # mother-infant unrelated
  mutate(unrelated_m_b = ifelse(Type_1 != Type_2 &
                                  FAMILY_1 != FAMILY_2, T, F)) %>%
  # M1 preservation
  mutate(longitudinal = ifelse(FAMILYupd_1 == FAMILYupd_2 &
                               Modified_NEXT_ID_without_preg_number_1 == Modified_NEXT_ID_without_preg_number_2, T, F))


#############################################################
# 3.1 Analysis: vOTU sharing to mom
#############################################################
caring <- sharing %>%
  filter(same_dyad == T | unrelated_m_b == T) %>%
  mutate(
    temp_s1 = sample_1,
    sample_1 = ifelse(Timepoint_new_1 != "Mother", sample_2, sample_1),
    sample_2 = ifelse(Timepoint_new_1 != "Mother", temp_s1, sample_2)
  ) %>%
  select(-temp_s1) %>%
  rename(mother = sample_1, infant = sample_2) %>%
  select(mother, infant, N_shared, same_dyad, unrelated_m_b) %>%
  pivot_longer(!c("mother", "infant", "N_shared")) %>%
  filter(value == T) %>%
  select(-value) %>%
  left_join(smeta, by = c("infant" = "Sequencing_ID")) %>%
  mutate(perc_shared = N_shared / vir_richness_cf * 100)

mom_share_stat <- map_dfr(c('same_dyad', 'unrelated_m_b'), function(mom){
  
  F1 <- as.formula(paste0("perc_shared ~ exact_age_months_at_collection + (1|Modified_NEXT_ID_without_preg_number)"))
  
  model <- lmer(F1, data = caring[caring$name == mom,])
  
  summary(model)$coefficients %>%
    as.data.frame() %>%
    rownames_to_column(var = "rowname") %>%
    filter(rowname != "(Intercept)") %>%
    mutate(Mother_relation = mom)
  
})


writexl::write_xlsx(mom_share_stat, '07.RESULTS/vOTU_sharing_to_mom.xlsx')

vOTU_share_mom <- caring %>%
  rename(Mother = name) %>%
  mutate(Mother = ifelse(Mother == "same_dyad", "Related", "Unrelated")) %>%
  ggplot(aes(Timepoint_new, perc_shared, fill = Mother, color = Mother)) +
  ggrastr::rasterise(geom_jitter(position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.9), size=0.01, alpha=1), dpi = 300)+
  scale_color_manual(values = c("#FFCAD4", "#B7C9F2")) +
  ggnewscale::new_scale_color() +
  geom_violin(aes(color = Mother),position = position_dodge(width = 0.9), alpha = 0.3, linewidth = 0.4, scale = "width") + 
  scale_color_manual(values = c("darkred", "darkblue")) +
  ggnewscale::new_scale_color() +
  geom_boxplot(aes(color = Mother),position = position_dodge(width = 0.9), alpha = 0.3, width = 0.3, outlier.shape = NA) +
  scale_color_manual(values = c("darkred", "darkblue")) +
  scale_fill_manual(values = c("#FFCAD4", "#B7C9F2")) +
  labs(x = "Timepoint", y="% of infant richness shared to mother") +
  theme_bw() +
  theme(axis.title = element_text(size=9),
        axis.text = element_text(size=8),
        legend.title = element_text(size=9),
        legend.text = element_text(size=8),
        legend.position = "bottom")
  
ggsave('05.PLOTS/06.DYNAMICS/vOTU_sharing_to_mom_over_time_VLP.pdf',
       vOTU_share_mom,  "pdf", width=6.5, height=8, units="cm", dpi = 300)

#############################################################
# 3.2 Analysis: vOTU drivers of convergence
#############################################################
votu_sparse <- Matrix(as.matrix(VLP > 0), sparse = TRUE)  # logical sparse

# maternal lookup
mom_votu <- lapply(unique(smeta$FAMILYupd), function(dyad) {
  mom_cols <- smeta$Sequencing_ID[smeta$Type == "M"]
  
  which(rowSums(votu_sparse[, mom_cols, drop = FALSE]) > 0) #indices, not names
})
names(mom_votu_sets) <- unique(smeta$FAMILYupd)

# infant sparse matrix
shared_sparse <- Matrix(0, nrow = nrow(votu_sparse), 
                        ncol = length(smeta$Sequencing_ID[smeta$Type == "K"]),
                        sparse = TRUE, dimnames = list(rownames(VLP), smeta$Sequencing_ID[smeta$Type == "K"]))

for (inf in smeta$Sequencing_ID[smeta$Type == "K"]) {
  
  dyad <- smeta$FAMILYupd[smeta$Sequencing_ID == inf]
  
  if (length(mom_votu_sets[[dyad]]) == 0) next # no mom
  
  mom_idx <- mom_votu_sets[[dyad]]           
  infant_col  <- votu_sparse[, inf, drop = FALSE] # sparse column

  shared_idx  <- mom_idx[as.numeric(infant_col[mom_idx, 1]) > 0]
  if (length(shared_idx) == 0) next
  shared_sparse[shared_idx, inf] <- 1L
  

}

# back to df:
shared_df <- as.data.frame(as.matrix(shared_sparse))

# exclude seq orphans:
seq_orphans <- smeta %>%
  group_by(FAMILYupd) %>%
  summarise(n_members = length(unique(NEXT_ID)),
            n_types = length(unique(Type))) %>%
  filter(n_members == 1 | n_types == 1) %>%
  pull(FAMILYupd)

shared_df <- shared_df[,colnames(shared_df) %in% (smeta %>% filter(!FAMILYupd %in% seq_orphans) %>% pull(Sequencing_ID))]
shared_df <- shared_df[rowSums(shared_df) >0,]

faster <- ETOF_vOTUr %>%
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
  ) %>%
  mutate(Host_Genus = ifelse(Host_Genus %in% c("g__Bacteroides", "g__Phocaeicola"), "Bacteroides_Phocaeicola", Host_Genus))

check_drivers <- shared_df %>%
  rownames_to_column(var = "New_CID") %>%
  left_join(faster) 

by_simple_host <- check_drivers %>%
  select(all_of(colnames(shared_df)), Host) %>%
  group_by(Host) %>%
  summarise(across(colnames(shared_df), sum)) %>%
  filter(Host != "Unknown") %>%
  pivot_longer(!Host, names_to = "Sequencing_ID") %>%
  left_join(smeta) %>%
  ggplot(aes(Timepoint_new, value, fill = Host)) + 
  geom_boxplot() # this happens primarily due to sharing bacteriophages


by_hostgen <- check_drivers %>%
  filter(Host == "bacteria") %>%
  select(all_of(colnames(shared_df)), Host_Genus) %>%
  group_by(Host_Genus) %>%
  summarise(across(colnames(shared_df), sum)) %>%
  filter(!Host_Genus %in% c("Unclassified", "Unassigned", "Unknown", "g__Unclassified")) %>%
  column_to_rownames(var = "Host_Genus")

# are bacteroides enriched in vOTU sharing?
VLP_infant <- VLP[,colnames(VLP) %in% colnames(shared_df)] # only infants that have moms
VLP_infant <- VLP_infant[rowSums(VLP_infant) > 0, ]

all_votus <- rownames(VLP_infant) 

shared_any <- rownames(shared_df)
is_shared <- all_votus %in% shared_any

is_bacteroides <- faster$Host_Genus[match(all_votus, faster$New_CID)] %in% c("Bacteroides_Phocaeicola")

fisher_bacteroides <- fisher.test(table(is_shared, is_bacteroides))
fisher_bacteroides$p.value

table(is_shared, is_bacteroides)

# are bacteroides phages the main driver?
  
dfer <- by_hostgen %>%
  rownames_to_column(var = "Host_Genus") %>%
  pivot_longer(!Host_Genus, names_to = "Sequencing_ID") %>%
  left_join(smeta) %>%
  left_join(caring %>% filter(name == "same_dyad") %>% select(infant, N_shared), by = c("Sequencing_ID" = "infant"))
  
results <- map_dfr(unique(dfer$Host_Genus), function(Host) {
  
  F1 <- as.formula("N_shared ~ value + Timepoint_new + (1|NEXT_ID)")
  
  model <- lmer(F1, data = dfer[dfer$Host_Genus == Host,])
  
  summary(model)$coefficients %>%
    as.data.frame() %>%
    rownames_to_column(var = "rowname") %>%
    filter(rowname == "value") %>%
    mutate(Host = Host)
  
})
results$FDR <- p.adjust(results$`Pr(>|t|)`, "BH")

prev <- rowSums(by_hostgen > 0) %>%
  as.data.frame() %>%
  rownames_to_column("Host") %>%
  rename("Prevalence"=2)

results <- results %>%
  left_join(prev) %>%
  select(-rowname) %>%
  relocate(Host, Prevalence)

writexl::write_xlsx(results, '07.RESULTS/Convergence_to_mom_Host_genus_drivers.xlsx')

results_filtered <- results %>%
  filter(FDR < 0.05, Prevalence >= 0.1*ncol(shared_df)) %>%
  mutate(
    ci_lo = Estimate - 1.96 * `Std. Error`,
    ci_hi = Estimate + 1.96 * `Std. Error`,
    host_clean = gsub("g__", "", Host),
    is_bact = grepl("Bacteroides", Host)
  ) %>%
  mutate(host_clean = ifelse(grepl("Bacteroides", host_clean), "Bacteroides/Phocaeicola", host_clean)) %>%
  mutate(host_clean = paste0("*", host_clean, "*")) %>%
  arrange(Estimate)

results_filtered_plot <- ggplot(results_filtered, aes(x = Estimate, y = reorder(host_clean, Estimate))) +
  geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi),
                 height = 0.2, linewidth = 0.4, alpha = 0.5,
                 color = ifelse(results_filtered$is_bact, "firebrick", "#132440")) +
  geom_point(aes(size = Prevalence),
             color = ifelse(results_filtered$is_bact, "firebrick", "#132440"),
             alpha = 0.8) +
  scale_size_continuous(range = c(1, 8), name = "Prevalence\n(N samples)") +
  labs(x = "Effect size (beta)", y = NULL) +
  theme_bw() +
  theme(
    axis.text.y = ggtext::element_markdown(size = 8),
    axis.title = element_text(size=9),
    legend.title = element_text(size=9),
    legend.text = element_text(size = 8),
    legend.position    = "bottom"
  )

ggsave('05.PLOTS/06.DYNAMICS/Convergence_drivers.pdf',
       results_filtered_plot, "pdf", width=13, height=13, units="cm", dpi = 300)

#############################################################
# 3.3 Analysis: Bacteroides bacteria expansion in the infant gut
#############################################################
prepik <- by_hostgen %>%
  rownames_to_column(var = "Host_Genus") %>%
  filter(Host_Genus == "Bacteroides_Phocaeicola") %>%
  pivot_longer(!Host_Genus, names_to = "Sequencing_ID", values_to = "N_bac_ph") %>%
  left_join(smeta)

host_bacteroides <- metaphlan %>%
  rownames_to_column(var = "Genus") %>%
  filter(grepl("g__Bacteroides", Genus)) %>%
  pivot_longer(!Genus, names_to = "corresp_MGS", values_to = "Bac_ab") %>%
  right_join(prepik)

# bacterial abundance as mediator between timepoint and phage sharing
med_fit <-lme4::lmer(Bac_ab ~ exact_age_months_at_collection + (1 | NEXT_ID), data = host_bacteroides)

out_fit <- lme4::lmer(N_bac_ph ~ exact_age_months_at_collection + Bac_ab + (1 | NEXT_ID),
                data = host_bacteroides)

med_out <- mediate(
  med_fit,
  out_fit,
  treat = "exact_age_months_at_collection",
  mediator = "Bac_ab",
  control.value = 0,
  treat.value = 12,
  sims = 1000
) # explains 12% only?

summary(lmer(Bac_ab ~ Timepoint_new + (1 | NEXT_ID), data = host_bacteroides))
