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
library(emmeans)
#library(Matrix)
#library(mediation) # messes dplyr
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
sharing_rel <- sharing %>%
  filter(same_dyad == T) %>%
  mutate(
    temp_s1 = sample_1,
    sample_1 = ifelse(Timepoint_new_1 != "Mother", sample_2, sample_1),
    sample_2 = ifelse(Timepoint_new_1 != "Mother", temp_s1, sample_2)
  ) %>%
  select(-temp_s1) %>%
  rename(mother = sample_1, infant = sample_2) %>%
  select(infant, N_shared) %>%
  left_join(smeta, by = c("infant" = "Sequencing_ID")) %>%
  mutate(perc_shared = N_shared / vir_richness_cf * 100) %>%
  mutate(name = "same_dyad")

# decided to take the mean, otherwise - inflated?

sharing_unr <- sharing %>%
  filter(unrelated_m_b == T) %>%
  mutate(
    temp_s1 = sample_1,
    sample_1 = ifelse(Timepoint_new_1 != "Mother", sample_2, sample_1),
    sample_2 = ifelse(Timepoint_new_1 != "Mother", temp_s1, sample_2)
  ) %>%
  select(-temp_s1) %>%
  rename(mother = sample_1, infant = sample_2) %>%
  select(mother, infant, N_shared, same_dyad, unrelated_m_b) %>%
  group_by(infant) %>%
  summarise(N_shared = round(mean(N_shared), 0)) %>%
  left_join(smeta, by = c("infant" = "Sequencing_ID")) %>%
  mutate(perc_shared = N_shared / vir_richness_cf * 100) %>%
  mutate(name = "unrelated_m_b")

caring <- bind_rows(sharing_rel, sharing_unr)

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
        legend.position = "bottom",
        legend.key.size = unit(0.7, 'line'))
  
ggsave('05.PLOTS/06.DYNAMICS/vOTU_sharing_to_mom_over_time_VLP.pdf',
       vOTU_share_mom,  "pdf", width=6.5, height=7, units="cm", dpi = 300)

#############################################################
# 3.2 Analysis: vOTU drivers of convergence
#############################################################
votu_sparse <- Matrix(as.matrix(VLP > 0), sparse = TRUE)  # logical sparse

# maternal lookup
mom_votu_sets <- lapply(unique(smeta$FAMILYupd), function(dyad) {
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

# saving shared_df (15,659 obs., 749 variables) for PPV-related analysis
write.table(shared_df, "06.CLEAN_DATA/Intermediate/Mom-shared_vOTUs_table.txt", sep='\t', quote=F)

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

prev <- rowSums(by_hostgen > 0) %>%
  as.data.frame() %>%
  rownames_to_column("Host") %>%
  rename("Prevalence"=2) %>%
  filter(Prevalence > 0.1*749)

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
  filter(Host_Genus %in% c(prev$Host, "Bacteroides_Phocaeicola")) %>%
  pivot_longer(!Host_Genus, names_to = "Sequencing_ID") %>%
  left_join(smeta) %>%
  left_join(caring %>% filter(name == "same_dyad") %>% select(infant, N_shared), by = c("Sequencing_ID" = "infant"))

model_full <- lmer(value ~ Host_Genus + Timepoint_new + (1 | NEXT_ID),
                   data = dfer,
                   REML = FALSE)

model_null <- lmer(value ~ Timepoint_new + (1 | NEXT_ID),
                   data = dfer,
                   REML = FALSE)

anova(model_full, model_null) # 61 339500 340027 -169689   339378 17844 55  < 2.2e-16 ***

emm_host <- emmeans(model_full, ~ Host_Genus, type = "response")

emm_host_df <- as.data.frame(emm_host) %>%
  arrange(desc(emmean)) %>%
  left_join(prev, by = c("Host_Genus" = "Host"))

###
emm_host <- emmeans(model_full, ~ Host_Genus)

host_levels <- emm_host@grid$Host_Genus

coef_bact <- ifelse(
  host_levels == "Bacteroides_Phocaeicola",
  1,
  -1 / (length(host_levels) - 1)
)

bact_vs_others <- contrast(
  emm_host,
  method = list("Bacteroides vs mean of other hosts" = coef_bact)
)

summary(bact_vs_others, infer = TRUE)
che <- as.data.frame(bact_vs_others)

writexl::write_xlsx(emm_host_df, '07.RESULTS/Convergence_to_mom_Host_genus_drivers.xlsx')

results_filtered_plot <- emm_host_df %>%
  mutate(host_clean = gsub("g__", "", Host_Genus)) %>%
  mutate(host_clean = ifelse(grepl("Bacteroides", host_clean), "Bacteroides/Phocaeicola", host_clean)) %>%
  mutate(host_clean = paste0("*", host_clean, "*")) %>%
  filter(Prevalence > 200) %>% # to keep top 10-most prevalent
ggplot(aes(x = emmean, y = reorder(host_clean, emmean))) +
  geom_errorbarh(aes(xmin = asymp.LCL, xmax = asymp.UCL),
                 height = 0.2, linewidth = 0.4, alpha = 0.5,
                 color = "#132440") +
  geom_point(aes(size = Prevalence),
             color = "#132440",
             alpha = 0.8) +
  scale_size_continuous(range = c(1, 6), name = "N infant samples sharing\nphage-host groups") +
  labs(x = "Adjusted mean value of shared phages", y = "Predicted host") +
  theme_bw() +
  theme(
    axis.text.y = ggtext::element_markdown(size = 8),
    axis.title = element_text(size=9),
    legend.title = element_text(size=9),
    legend.text = element_text(size = 8),
    legend.position = "bottom"
  )

ggsave('05.PLOTS/06.DYNAMICS/Convergence_drivers.pdf',
       results_filtered_plot, "pdf", height=13, width=12, units="cm", dpi = 300)

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

med_out <- mediation::mediate( # calling mediation this way because it messes dplyr
  med_fit,
  out_fit,
  treat = "exact_age_months_at_collection",
  mediator = "Bac_ab",
  control.value = 0,
  treat.value = 12,
  sims = 1000
) # explains 12% only?

summary(med_out)

summary(lmer(log(Bac_ab + 0.0001)~ Timepoint_new + (1 | NEXT_ID), data = host_bacteroides))$coefficients

min(host_bacteroides$Bac_ab[host_bacteroides$Bac_ab > 0])
