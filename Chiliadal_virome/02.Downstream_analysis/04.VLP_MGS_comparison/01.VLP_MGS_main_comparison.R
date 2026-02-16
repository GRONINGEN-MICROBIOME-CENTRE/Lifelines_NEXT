setwd('~/Desktop/Projects_2021/NEXT_virome/09.DATA_ANALYSIS/')

#############################################################
# Comparison of MGS and VLP metaviromes:
# - at the level of vOTU recovery
# - at the population level
#############################################################

#############################################################
# 0. Used files source
#############################################################

#############################################################
# 1. Functions
#############################################################
# rounds
smart_round <- function(x){
  
  ifelse(x < 0.01,
        as.numeric(formatC(x, format = "e", digits = 2)),
                round(x, 2))

}

# identifies the DB origin of the vOTU cluster member
# writes it to column "DB_member"
get_member_origin <- function(df, member_col, dbs, class_col="DB_member"){
  
  df[class_col] <- NA
  
  for (db in dbs) {
    
    idx <- grepl(db, df[[member_col]], fixed = T)
    
    df[[class_col]][idx] <- gsub('-', '', db, fixed = T)
    
  }
  df[[class_col]][grepl('NEXT_V', df[[member_col]], fixed = T)] <- "NEXT_VLP"
  df[[class_col]][grepl('NEXT_M', df[[member_col]], fixed = T)] <- "NEXT_MGS"
  
  df
}

#############################################################
# 1. Loading libraries
#############################################################
library(tidyverse)
library(ggplot2)
library(dplyr)
library(see)
library(UpSetR)
library(lme4)
library(lmerTest)
library(purrr)
library(broom)
library(patchwork)
#############################################################
# 2. Load Input Data
#############################################################

clean_RPKM <- read.table('06.CLEAN_DATA/02.FINAL/RPKM_table_VLP_MGS_dec99ANI_ALL_CS_ab3kbp_2220_samples.txt', sep='\t', header=T)

smeta <- read.delim('06.CLEAN_DATA/02.FINAL/Chiliadal_meta_VLP_MGS_matched_v05_suppl_w_virmetrics.txt', sep='\t', header=T)
smeta$Timepoint_new <- factor(smeta$Timepoint_new, levels=c("M1", "M3", "M6", "M12", "Mother"), ordered = T)

ETOF <- read.table('06.CLEAN_DATA/VLP_MGS_ETOF_full_rep.txt', sep='\t', header = T)

ETOF_vOTUr <- read.table('06.CLEAN_DATA/02.FINAL/Working_ETOF_120997vOTUr_ab3kbp_in_2200_VLP_MGS.txt', sep='\t', header=T)

vOTU_clustering <- read.table('06.CLEAN_DATA/NEXT_viral_clusters_MGS_VLP_long_format.txt', sep='\t', header=T)
keeper_clustering <- vOTU_clustering
vOTU_clustering <- vOTU_clustering[vOTU_clustering$Representative %in% ETOF_vOTUr$New_CID,]

vOTU_cluster_size <- read.table('06.CLEAN_DATA/NEXT_viral_clusters_MGS_VLP_size.txt', sep='\t', header=T)
vOTU_cluster_size <- vOTU_cluster_size[vOTU_cluster_size$Representative %in% ETOF_vOTUr$New_CID,]

VLP <- read.table('06.CLEAN_DATA/02.FINAL/VLP_only_RPKM_table_VLP_MGS_dec99ANI_ALL_CS_ab3kbp_1110_samples.txt', sep='\t', header=T)

MGS <- read.table('06.CLEAN_DATA/02.FINAL/MGS_only_RPKM_table_VLP_MGS_dec99ANI_ALL_CS_ab3kbp_1110_samples.txt', sep='\t', header=T)
#############################################################
# 3.1 Analysis: staton on MGS vs VLP genomics
#############################################################
genomics_compare <- c('raw_reads', 'human_reads', 'clean_reads',
                      'contigs_0_bp', 'contigs_1000_bp', 'sc_enrichment')

results_genomics <- map_dfr(genomics_compare, function(genomics) {
    
    print(paste("Data:", genomics))
  
  f1 <- paste("log10(", genomics, ")", "~ seq_type")
    
    formula <- as.formula(paste0(f1, "+ Timepoint_new + (1 | NEXT_ID)"))
    
    model <- lmer(
      formula,
      REML = FALSE,
      data = smeta
    )
   
    eff_res <- effsize::cohen.d(as.formula(f1), data = smeta) 
    
    model_summary <- summary(model)$coefficients %>%
      as.data.frame() %>%
      rownames_to_column() %>%
      filter(rowname != "(Intercept)") %>%
      mutate(Cohens_D = abs(eff_res$estimate)) %>%
      mutate(rowname = gsub('seq_type', '', rowname)) %>%
      rename(`Metavirome type` = rowname) %>%
      mutate(`Genomic metric` = genomics) %>%
      relocate(`Genomic metric`)
    
    model_summary
    
  }) 

results_genomics <- results_genomics %>%
  mutate(p_adj = p.adjust(`Pr(>|t|)`, "BH")) %>%
  mutate(across(where(is.numeric), ~ smart_round(.)) )

write.table(results_genomics, '07.RESULTS/Compare_genomic_metrics_MGS_VLP.txt', sep='\t', quote=F, row.names=F)
rm(genomics_compare)
#############################################################
# 3.2 Analysis: VLP and MGS contribution to the NEXT virome
#############################################################

# all DBs used for the dereplication:
dbs <- c('NEXT_V', 'NEXT_M', 'Guerin', 
         'Yutin', 'NCBI_CrAss', 'NL_crAss',
         'Benler', 'COPSAC', 'GVD', 
         'IMGVR', 'MGV-', 'GPD', 'VREF')

vOTU_clustering <- get_member_origin(vOTU_clustering, "Cluster_member", dbs)

# calculating db contributions:
vOTU_by_source <- vOTU_clustering %>%
  mutate(DB_member = factor(DB_member)) %>%
  count(Representative, DB_member) %>%
  pivot_wider(names_from = DB_member, values_from = n, values_fill = 0) %>%
  left_join(vOTU_cluster_size, by = "Representative") %>%
  rename(N_genomes = Cluster_size, vOTU_representative  = Representative) %>%
  left_join(ETOF_vOTUr %>% select(New_CID, POST_CHV_length, miuvig_quality, vOTU_cluster_type),
            by = c("vOTU_representative" = "New_CID"))

# cumulative external db input:
vOTU_by_source$DB_sum <- vOTU_by_source$N_genomes - vOTU_by_source$NEXT_MGS - vOTU_by_source$NEXT_VLP

listInput <- list(DB=vOTU_by_source[vOTU_by_source$DB_sum>0,]$vOTU_representative,
                  MGS=vOTU_by_source[vOTU_by_source$NEXT_MGS>0,]$vOTU_representative,
                  VLP=vOTU_by_source[vOTU_by_source$NEXT_VLP>0,]$vOTU_representative)


upset_all <- upset(fromList(listInput), order.by = "freq", sets.bar.color = c("#71C9CE", "#A6E3E9", "#CBF1F5"), 
                   number.angles = 0,
                   sets.x.label = "N vOTUs", scale.sets = "identity",
                   text.scale = c(1, 1, 1, 0.7, 1, 1))

png('05.PLOTS/05.VLP_MGS/UpSet_plot_all_vOTUs.png', width=10, height=7, units="cm", res = 300)
upset_all
dev.off()

pdf('05.PLOTS/05.VLP_MGS/UpSet_plot_all_vOTUs.pdf', width=4, height=3.5)
upset_all
dev.off()

catcont <- ETOF_vOTUr %>%
  group_by(vOTU_cluster_type) %>%
  summarise(sum = n()) %>%
  mutate(perc = sum/nrow(ETOF_vOTUr) * 100)

rm(dbs)
#############################################################
# 3.3 Analysis: VLP vs MGS genome recovery taking into
# account QoL
#############################################################

# sub-analysis: how frequent it is for vOTUs w "Genome-fragment" miuvig quality" to have high quality members
vOTU_clustering <- vOTU_clustering %>%
  left_join(ETOF %>% select(New_CID, miuvig_quality, POST_CHV_length), by = c("Cluster_member" = "New_CID"))

cluster_summary <- vOTU_clustering %>%
  group_by(Representative) %>%
  summarise(
    N_high = sum(miuvig_quality == "High-quality"),
    N_low  = sum(miuvig_quality == "Genome-fragment")
  ) %>%
  left_join(ETOF_vOTUr %>% select(New_CID, miuvig_quality, POST_CHV_length), by=c("Representative" = "New_CID"))

dim(cluster_summary[cluster_summary$N_high > 0 & cluster_summary$miuvig_quality=="Genome-fragment",])
#there are only 229 cases like that, some of them are caused by CheckV erroneously
##considering some short genomes HQ (as w vOTU==NEXT_V0156_N33_L31446_K2.8_E0_P0_F0)
#######DECISION: IGNORE IT

### hq only:
HQ_vOTU_by_source <- vOTU_by_source %>%
  filter(miuvig_quality == 'High-quality')

# UpSet plot (HQ vOTUs)
listInputHQ <- list(DB=HQ_vOTU_by_source[HQ_vOTU_by_source$DB_sum > 0,]$vOTU_representative,
                    MGS=HQ_vOTU_by_source[HQ_vOTU_by_source$NEXT_MGS > 0,]$vOTU_representative,
                    VLP=HQ_vOTU_by_source[HQ_vOTU_by_source$NEXT_VLP > 0,]$vOTU_representative)


upset_allhq <- upset(fromList(listInputHQ), order.by = "freq", sets.bar.color = c("#71C9CE", "#A6E3E9", "#CBF1F5"), 
                     number.angles = 0,
                     sets.x.label = "N (near-)complete\nvOTUs", scale.sets = "identity",
                     text.scale = c(1, 1, 1, 0.7, 1, 1))

png('05.PLOTS/05.VLP_MGS/UpSet_plot_HQ_vOTUs.png', width=8, height=7, units="cm", res = 300)
upset_allhq
dev.off()

pdf('05.PLOTS/05.VLP_MGS/UpSet_plot_HQ_vOTUs.pdf', width=4, height=3.5)
upset_allhq
dev.off()

clustering_hq_sum <- HQ_vOTU_by_source %>%
  group_by(vOTU_cluster_type) %>%
  summarise(sum = n(), .groups = "drop") %>%
  mutate(perc = sum/nrow(HQ_vOTU_by_source) * 100)

sum(clustering_hq_sum$sum[!grepl('DB', clustering_hq_sum$vOTU_cluster_type)])

### CALCUATING TRUE NOVELS taking into account the QoL:
# derive the Quality of the Longest (QoL) contig per contributing source per vOTUr,
# to truly understand which source is indispensible
sources <- c("DB_QoL", "MGS_QoL", "VLP_QoL")

hq_recovery_source <- vOTU_clustering %>%
  filter(Representative %in% HQ_vOTU_by_source$vOTU_representative) %>%
  mutate(DB_member = if_else(DB_member %in% c("NEXT_MGS", "NEXT_VLP"), DB_member, "DB")) %>%
  arrange(desc(POST_CHV_length), desc(miuvig_quality)) %>%
  group_by(Representative, DB_member) %>%
  slice_head(n = 1) %>%    
  ungroup() %>%
  select(Representative, DB_member, miuvig_quality) %>%
  pivot_wider(
    names_from = DB_member,
    values_from = miuvig_quality,
    names_prefix = ""
  ) %>%
  rename(DB_QoL = DB, MGS_QoL=NEXT_MGS, VLP_QoL=NEXT_VLP) %>%
  replace_na(list(
    DB_QoL = "Absent",
    MGS_QoL = "Absent",
    VLP_QoL = "Absent"
  )) %>% # at this point: per-vOTU (HQ) whether the longest sequence coming from 3 sources is also of High-quality
  mutate(recoverable_by = apply(across(all_of(sources)), 1, function(r) {
    matches <- sources[r == "High-quality"]
    if(length(matches) > 0) paste(matches, collapse = "+") else NA
  })) %>%
  left_join(ETOF_vOTUr %>% select(New_CID, vOTU_cluster_type, genome), by = c("Representative" = "New_CID")) 

saver_hq_recovery <- hq_recovery_source 

# UpSet plot (HQ vOTUs recoverability)
listInputHQrecover <- list(DB=hq_recovery_source[hq_recovery_source$DB_QoL=="High-quality",]$Representative,
                    MGS=hq_recovery_source[hq_recovery_source$MGS_QoL=="High-quality",]$Representative,
                    VLP=hq_recovery_source[hq_recovery_source$VLP_QoL=="High-quality",]$Representative)


upset_allhq_rec <- upset(fromList(listInputHQrecover), order.by = "freq", sets.bar.color = c("#71C9CE", "#A6E3E9", "#CBF1F5"), 
                     number.angles = 0,
                     sets.x.label = "N (near-)complete\nvOTUs", scale.sets = "identity",
                     text.scale = c(1, 1, 1, 0.7, 1, 1))

png('05.PLOTS/05.VLP_MGS/UpSet_plot_HQ_vOTUs_recovery.png', width=8, height=7, units="cm", res = 300)
upset_allhq_rec
dev.off()


summary_hq_recovery <- hq_recovery_source %>%
  filter(!grepl('DB', vOTU_cluster_type)) %>%
  group_by(recoverable_by) %>%
  summarise(sum = n()) %>%
  mutate(perc = sum/sum(sum)*100)

#### genome composition of these:

by_genome_recovery <- hq_recovery_source %>%
  filter(!grepl('DB', vOTU_cluster_type)) %>%
  group_by(recoverable_by, genome) %>%
  summarise(sum = n(), .groups = "drop_last") %>%
  mutate(within_rec_perc = sum/sum(sum)*100) %>%
  ungroup()
#############################################################
# 3.4 Analysis: VLP vs MGS vOTU richness
#############################################################
rich_abs <- map_dfr(c("VLP", "MGS"), function(method){
  
 map_dfr(c("vir_richness_cf", "temperate_richness"), function(Richness_type){
   
   rich_cf <- smeta %>%
     filter(seq_type == method) %>%
     pull(Richness_type) 
   
   rich_cf %>%
     summary(na.rm =T) %>%
     tidy() %>% 
     mutate(method = method,
            richness_type = Richness_type,
            sd = sd(rich_cf, na.rm = T))
   
 })
  
})

rich_abs <- rich_abs %>%
  mutate(across(where(is.numeric), ~ round(., 0)))


rich_abs_diff <- map_dfr(c("vir_richness_cf", "temperate_richness"), function(Richness_type){
  
  rich_diff <- smeta %>%
  select(Universal_ID, Richness_type, seq_type) %>%
    pivot_wider(names_from = seq_type,
                values_from = Richness_type) %>%
    mutate(difference = MGS - VLP) %>%
    pull(difference) 
  
  rich_diff %>%
    summary(difference, na.rm = T) %>%
    tidy() %>% 
    mutate(richness_type = Richness_type,
           sd = sd(rich_diff, na.rm = T))
  
})

diffs_assoc <- smeta %>%
  select(Universal_ID, vir_richness_cf, temperate_richness, seq_type, NEXT_ID, Timepoint_new) %>%
  pivot_wider(names_from = seq_type,
              values_from = c(vir_richness_cf, temperate_richness)) %>%
  mutate(difference_total = vir_richness_cf_MGS - vir_richness_cf_VLP,
         difference_temperate = temperate_richness_MGS - temperate_richness_VLP)


models <- c(
  "VLP ~ MGS" = "vir_richness_cf_VLP ~ vir_richness_cf_MGS + Timepoint_new + (1|NEXT_ID)",
  "difference_total ~ difference_temperate" = "difference_total ~ difference_temperate + Timepoint_new + (1|NEXT_ID)"
)

model_sums <- imap_dfr(models, function(form, name) {
  lmer(as.formula(form), data = diffs_assoc, REML = FALSE) %>%
    summary() %>%
    .$coefficients %>%
    as.data.frame() %>%
    rownames_to_column("rowname") %>%
    filter(rowname != "(Intercept)", !grepl('Timepoint', rowname)) %>%
    mutate(model = name) %>%
    relocate(model)
})

model_sums <- model_sums %>%
  mutate(across(where(is.numeric), ~ smart_round(.)))

rich_abs_compare <- map_dfr(c('vir_richness_cf', 'temperate_richness'), function(rich_type) {
  
  f1 <- paste0(rich_type, "~ seq_type")
  formula <- as.formula(paste0(f1, "+ Timepoint_new + (1|NEXT_ID)"))
  
  model <- lmer(
    formula,
    REML = FALSE,
    data = smeta
  )
  
  eff_res <- effsize::cohen.d(as.formula(f1), data = smeta) 
  
  model_summary <- summary(model)$coefficients %>%
    as.data.frame() %>%
    rownames_to_column() %>%
    filter(rowname != "(Intercept)") %>%
    mutate(Cohens_D = abs(eff_res$estimate)) %>%
    mutate(rowname = gsub('type.*', 'type', rowname)) %>%
    filter(!grepl('Timepoint', rowname)) %>%
    mutate(richness_type = rich_type)
  
})

rich_abs_compare <- rich_abs_compare %>%
  mutate(p_adj = p.adjust(`Pr(>|t|)`, "BH")) %>%
  mutate(across(where(is.numeric), ~ smart_round(.)))

writexl::write_xlsx(rich_abs_compare, '07.RESULTS/Compare_sample_richness_by_seq_type.xlsx')
writexl::write_xlsx(model_sums, '07.RESULTS/Compare_assoc_richness_temperate_richness.xlsx')

vOTU_rich <- smeta %>%
  mutate(Richness_type = "vOTU richness") %>%
  ggplot(aes(seq_type, vir_richness_cf)) +
  geom_violinhalf(aes(seq_type, vir_richness_cf, fill = seq_type, color=seq_type), scale = "width", flip=1, lwd = 0.4) +
  geom_boxplot(aes(seq_type, vir_richness_cf, fill = seq_type, color=seq_type),width = 0.2, outlier.shape = NA, lwd = 0.4) +
  geom_rect(xmin=1.002, xmax=1.5, ymin = -0.001, ymax=4000, fill="white") +
  geom_rect(xmin=1.5, xmax=1.998, ymin= -0.001, ymax=4000, fill="white") +
  geom_line(aes(group = Universal_ID),alpha=0.1, color="darkgrey", linewidth = 0.1) +
  geom_point(aes(seq_type, vir_richness_cf, color=seq_type), size=0.3, alpha=0.2) + 
  facet_wrap(~Richness_type, scales = "free") +
  ggsignif::geom_signif(comparisons = list(c("MGS", "VLP")),
                        map_signif_level = TRUE, textsize = 4) +
  theme_bw() +
  theme(legend.position = "none",
        strip.background = element_rect(NA),
        legend.text = element_text(size=7),
        legend.title = element_text(size = 8),
        axis.title = element_text(size=8),
        axis.text = element_text(size=7),
        strip.text = element_text(size = 8)) +
  labs(y = "N vOTUs", x = "Metavirome type") +
  scale_fill_manual(values = c(MetBrewer::met.brewer("Kandinsky")[1], 
                               MetBrewer::met.brewer("Kandinsky")[2])) +
  scale_color_manual(values = c("#F4991A", "#492828")) +
  ylim(c(0,4700))

temp_rich <- smeta %>%
  mutate(Richness_type = "Temperate\nvOTU richness") %>%
  ggplot(aes(seq_type, temperate_richness)) +
  geom_violinhalf(aes(seq_type, temperate_richness, fill = seq_type, color=seq_type), scale = "width", flip=1, lwd = 0.4) +
  geom_boxplot(aes(seq_type, temperate_richness, fill = seq_type, color=seq_type),width = 0.2, outlier.shape = NA, lwd = 0.4) +
  geom_rect(xmin=1.002, xmax=1.5, ymin = -0.001, ymax=800, fill="white") +
  geom_rect(xmin=1.5, xmax=1.998, ymin= -0.001, ymax=800, fill="white") +
  geom_line(aes(group = Universal_ID),alpha=0.1, color="darkgrey", linewidth = 0.1) +
  geom_point(aes(seq_type, temperate_richness, color=seq_type), size=0.3, alpha=0.2) + 
  facet_wrap(~Richness_type, scales = "free") +
  ggsignif::geom_signif(comparisons = list(c("MGS", "VLP")),
                        map_signif_level = TRUE, textsize = 4) +
  theme_bw() +
  theme(legend.position = "none",
        strip.background = element_rect(NA),
        legend.text = element_text(size=7),
        legend.title = element_text(size = 8),
        axis.title = element_text(size=8),
        axis.text = element_text(size=7),
        strip.text = element_text(size = 8)) +
  labs(y = "N vOTUs", x = "Metavirome type") +
  scale_fill_manual(values = c(MetBrewer::met.brewer("Kandinsky")[1], 
                               MetBrewer::met.brewer("Kandinsky")[2])) +
  scale_color_manual(values = c("#F4991A", "#492828")) +
  ylim(c(0,860))
  

both_rich <- (vOTU_rich | temp_rich ) + plot_layout(axis_titles = "collect")

ggsave('05.PLOTS/05.VLP_MGS/Richness_compare.png',
       both_rich,  "png", width=12, height=12, units="cm", dpi = 300)

ggsave('05.PLOTS/05.VLP_MGS/Richness_compare.pdf',
       both_rich,  "pdf", width=8, height=7, units="cm", dpi = 300)
#############################################################
# 3.5 Analysis: VLP vs MGS vOTU contribution to sample 
# composition
#############################################################
VLP_by_vOTU_type <- VLP %>%
  rownames_to_column("New_CID") %>%
  left_join(ETOF_vOTUr %>% select(New_CID, vOTU_cluster_type)) %>% 
  group_by(vOTU_cluster_type) %>%
  summarise(across(where(is.numeric), ~ sum(.x > 0)), .groups = "drop") %>% 
  mutate(across(where(is.numeric), ~ .x / sum(.x))) %>%
  mutate(vOTU_cluster_type = gsub("_sum|NEXT_", "", vOTU_cluster_type)) %>%
  pivot_longer(!"vOTU_cluster_type",
               names_to = "Sample",
               values_to = "Perc_by_type") %>%
  mutate(`Metavirome type` = "VLP")

MGS_by_vOTU_type <- MGS %>%
  rownames_to_column("New_CID") %>%
  left_join(ETOF_vOTUr %>% select(New_CID, vOTU_cluster_type)) %>% 
  group_by(vOTU_cluster_type) %>%
  summarise(across(where(is.numeric), ~ sum(.x > 0)), .groups = "drop") %>% 
  mutate(across(where(is.numeric), ~ .x / sum(.x))) %>%
  mutate(vOTU_cluster_type = gsub("_sum|NEXT_", "", vOTU_cluster_type)) %>%
  pivot_longer(!"vOTU_cluster_type",
               names_to = "Sample",
               values_to = "Perc_by_type") %>%
  mutate(`Metavirome type` = "MGS")


wo_db_stat <- map_dfr(list(VLP_by_vOTU_type, MGS_by_vOTU_type), function(dataset){
  
  wo_db <- dataset %>%
    filter(!grepl('DB', vOTU_cluster_type)) %>%
    group_by(Sample) %>%
    summarise(wo = sum(Perc_by_type), .groups = "drop")%>%
    pull(wo) 
  
  wo_db %>%
    summary(na.rm =T) %>%
    tidy() %>% 
    mutate(method = unique(dataset$`Metavirome type`),
           vOTU_cluster_type = "No_DB_at_all",
           sd = sd(wo_db, na.rm = T))
  
})

sum_by_type <- map_dfr(list(VLP_by_vOTU_type, MGS_by_vOTU_type), function(dataset){
  
  map_dfr(unique(VLP_by_vOTU_type$vOTU_cluster_type), function(vOTU_type){
    
    jahsd <- dataset %>%
      filter(vOTU_cluster_type == vOTU_type) %>%
      pull(Perc_by_type) 
    
    jahsd %>%
      summary(na.rm =T) %>%
      tidy() %>% 
      mutate(method = unique(dataset$`Metavirome type`),
             vOTU_cluster_type = vOTU_type,
             sd = sd(jahsd, na.rm = T))
    
  }
    
  )
  
})

sum_by_type <- sum_by_type %>%
  bind_rows(wo_db_stat) %>%
  mutate(across(where(is.numeric), ~ smart_round(. * 100)) )

# testing all vs all
rich_compare <- VLP_by_vOTU_type %>%
  bind_rows(MGS_by_vOTU_type) %>%
  mutate(sid = paste0(`Metavirome type`, '_',Sample)) %>%
  left_join(smeta %>% 
              mutate(sid = paste0(seq_type, '_', Universal_ID)) 
            %>% select(sid, Timepoint_new, NEXT_ID), by = c('sid' = "sid")) %>%
  select(-sid)

combos <- as.data.frame(t(combn(c(paste0('VLP_', unique(VLP_by_vOTU_type$vOTU_cluster_type)),
                                  paste0('MGS_', unique(VLP_by_vOTU_type$vOTU_cluster_type))), 2))) %>%
  mutate(method1 = gsub("_.*", "", V1),
         method2 = gsub("_.*", "", V2)
         ) %>%
  relocate(method1, V1, method2, V2) %>%
  mutate(Estimate = NA,
         `Std. Error` = NA,
         df = NA,
         `t value` = NA,
         `Pr(>|t|)` = NA,
         Cohens_D = NA)

for (combo in 1:nrow(combos)) {
  
  feature1 <-combos[combo, "V1"]
  method1 <- combos[combo, "method1"]
  feature2 <- combos[combo, "V2"]
  method2 <- combos[combo, "method2"]
  
  dat <- rich_compare %>%
    mutate(vOTU_cluster_type = paste0(`Metavirome type`, '_', vOTU_cluster_type)) %>%
    filter(vOTU_cluster_type == feature1 & `Metavirome type` == method1) %>%
    bind_rows(rich_compare %>% 
                mutate(vOTU_cluster_type = paste0(`Metavirome type`, '_', vOTU_cluster_type)) %>% 
                filter(vOTU_cluster_type == feature2 & `Metavirome type` == method2))
  
  print(paste("Testing", feature1,  method1, "versus", feature2, method2))
  
  f1 <- paste("Perc_by_type ~ vOTU_cluster_type")
  
  formula <- as.formula(paste0(f1, "+ Timepoint_new + (1 | NEXT_ID)"))

  model <- lmer(
    formula,
    REML = FALSE,
    data = dat
  )
  
  eff_res <- effsize::cohen.d(as.formula(f1), data = dat) 
  
  model_summary <- summary(model)$coefficients %>%
    as.data.frame() %>%
    rownames_to_column() %>%
    filter(rowname != "(Intercept)") %>%
    mutate(Cohens_D = abs(eff_res$estimate)) %>%
    mutate(rowname = gsub('type.*', 'type', rowname)) %>%
    filter(!grepl('Timepoint', rowname))
  
  combos[combo, 5:10] <- model_summary[-1]
}

combos_round <- combos %>%
  mutate(V1 = gsub('.*_', '', V1),
         V2 = gsub('.*_', '', V2)) %>%
  mutate(p_adj = p.adjust(`Pr(>|t|)`, "BH")) %>%
  mutate(across(where(is.numeric), ~ smart_round(.)) )
  

writexl::write_xlsx(combos_round, '07.RESULTS/Compare_contrib_sample_richness_by_vOTU_source.xlsx')

dat_text <- data.frame(
  "Metavirome type" = c("MGS", "MGS", "VLP", "VLP"),
  start = c("VLP+MGS", "MGS", "VLP+MGS", "VLP"),
  end = c("MGS", "VLP",  "VLP", "MGS"),
  y = c(0.8, 0.9, 0.8, 0.9),
  label = c("***", "***", "n.s.", "***")
) %>%
  rename( "Metavirome type" = Metavirome.type)

prop_richness <- VLP_by_vOTU_type %>%
  bind_rows(MGS_by_vOTU_type) %>%
  mutate(vOTU_cluster_type = factor(vOTU_cluster_type, levels = c("VLP+MGS", "MGS", "VLP", "VLP+MGS+DB", "MGS+DB", "DB", "VLP+DB"), ordered = T)) %>%
  ggplot(aes(vOTU_cluster_type, Perc_by_type)) +
  ggrastr::rasterise(geom_jitter(aes(color = vOTU_cluster_type), width = 0.4, alpha = 0.4, size = 1), dpi = 300) +
  geom_boxplot(aes(fill = vOTU_cluster_type), outlier.shape = NA, alpha = 0.5) + 
  facet_wrap(~`Metavirome type`) + 
  labs(y = "Proportion of sample richness", x = "Cluster Type", fill = "Cluster Type", color = "Cluster Type") +
  scale_color_manual(values = MetBrewer::met.brewer("Cross")) + 
  scale_fill_manual(values = MetBrewer::met.brewer("Cross")) +
  theme_bw() +
  theme(legend.position = "none", 
        strip.background = element_rect(NA),
        axis.text.x = element_text(angle=30, hjust = 0.8, vjust = 0.9),
        axis.text = element_text(size=7),
        axis.title = element_text(size=8),
        strip.text = element_text(size=8)) +
  ggsignif::geom_signif(data = dat_text,
                        aes(xmin = start, xmax = end, annotations = label, y_position = y),
                        textsize = 4,
                        manual = T)

ggsave('05.PLOTS/05.VLP_MGS/Proportion_richness_by_source.png',
       prop_richness,  "png", width=16, height=12, units="cm", dpi = 300)

ggsave('05.PLOTS/05.VLP_MGS/Proportion_richness_by_source.pdf',
       prop_richness,  "pdf", width=12, height=10, units="cm", dpi = 300)

#############################################################
# 3.6 Analysis: HQ richness in MGS vs VLP
#############################################################

# For HQs:
hq_votus <- ETOF_vOTUr[ETOF_vOTUr$miuvig_quality=="High-quality",]$New_CID
smeta$vir_richness_hq <- colSums(clean_RPKM[row.names(clean_RPKM) %in% hq_votus,] > 0)[match(smeta$Sequencing_ID, colnames(clean_RPKM))]

amq_votus <- ETOF_vOTUr$New_CID[ETOF_vOTUr$checkv_quality %in% c("Complete", "High-quality", "Medium-quality")]
smeta$vir_richness_amq <- colSums(clean_RPKM[row.names(clean_RPKM) %in% amq_votus,] > 0)[match(smeta$Sequencing_ID, colnames(clean_RPKM))]


ggplot(smeta[smeta$seq_type == "MGS",], aes(Timepoint_new, vir_richness_amq)) + 
  geom_boxplot()

p_VLP_MGS_richnessHQ <- smeta %>%
  mutate(Richness_type = "(Near-)complete vOTUs richness") %>%
  ggplot(aes(seq_type, vir_richness_hq)) +
  geom_violinhalf(aes(seq_type, vir_richness_hq, fill = seq_type, color=seq_type), flip=1) +
  geom_boxplot(aes(seq_type, vir_richness_hq, fill = seq_type, color=seq_type),width = 0.1, outlier.shape = NA) +
  geom_rect(xmin=1, xmax=1.1, ymin=10, ymax=50, fill="white") +
  geom_rect(xmin=1.9, xmax=2, ymin=10, ymax=100, fill="white") +
  geom_line(aes(group = Universal_ID),alpha=0.1, color="darkgrey") +
  geom_point(aes(seq_type, vir_richness_hq, color=seq_type), size=1.5, alpha=0.3) + 
  labs(y="vir_richness_hq of HQ vOTUs") +
  scale_x_discrete(labels = c("MGS", "VLP")) +
  theme_bw() +
  theme(legend.position = "none") +
  facet_wrap(~Richness_type, scales = "free") +
  ggsignif::geom_signif(comparisons = list(c("MGS", "VLP")),
                        map_signif_level = TRUE, textsize = 5) +
  theme_bw() +
  theme(legend.position = "none",
        strip.background = element_rect(NA)) +
  labs(y = "N vOTUs", x = "Metavirome type") +
  scale_fill_manual(values = c(MetBrewer::met.brewer("Kandinsky")[1], 
                               MetBrewer::met.brewer("Kandinsky")[2])) +
  scale_color_manual(values = c("#F4991A", "#492828")) +
  ylim(c(0,150))

ggsave('05.PLOTS/05.VLP_MGS/HQ_richness_comparison.png',
       p_VLP_MGS_richnessHQ,  "png", width=10, height=10, units="cm", dpi = 300)

modelhq <- lmer(vir_richness_hq ~ seq_type + Timepoint_new + (1|NEXT_ID),
  REML = FALSE,
  data = smeta
)

summary(modelhq)$coefficients # beta = 8.1, p-value = 5.6e-47

# explanation why HQs are higher in VLP vs MGS:
hq_votus_ssDNA <- ETOF_vOTUr[ETOF_vOTUr$miuvig_quality=="High-quality" &
                               ETOF_vOTUr$genome=="ssDNA",]$New_CID
smeta$clean_richness_HQss <- colSums(clean_RPKM[row.names(clean_RPKM) %in% hq_votus_ssDNA,] > 0)[match(smeta$Sequencing_ID, colnames(clean_RPKM))]

modelhqss <- lmer(clean_richness_HQss ~ seq_type + Timepoint_new + (1|NEXT_ID),
                REML = FALSE,
                data = smeta
)

summary(modelhqss)$coefficients # beta = 9.1, p-value = 5.8e-153

table(ETOF_vOTUr[ETOF_vOTUr$vOTU_cluster_type=="NEXT_VLP" &
                   ETOF_vOTUr$miuvig_quality=="High-quality" &
                   ETOF_vOTUr$genome %in% c('ssDNA', 'RNA'), "genome"]) #122 RNA, 2773 ssDNA

  
#############################################################
# 4. OUTPUT
#############################################################
# table with per-contig HQ vOTU recovery info by source
write.table(hq_recovery_source, './06.CLEAN_DATA/03.RESULTS_to_plot/HQ_vOTU_recoverable_by_source.txt', sep='\t', quote=F, row.names=F)



