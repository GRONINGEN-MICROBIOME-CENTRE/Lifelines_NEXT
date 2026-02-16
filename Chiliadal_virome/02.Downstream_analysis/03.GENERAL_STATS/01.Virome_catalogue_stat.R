setwd('~/Desktop/Projects_2021/NEXT_virome/09.DATA_ANALYSIS/')

#############################################################
# Descriptive statistics on NEXT virome catalogue
# Staton on novels
#############################################################

#############################################################
# 0. Used files source
#############################################################

#############################################################
# 1. Functions
#############################################################
get_iqr <- function(vec) {
  
  SMR <- summary(vec)
  
  print(paste0(round(SMR[3], 2), ' (', round(SMR[2], 2), ' - ', round(SMR[5], 2), ')')) 
  
}

# saturation_stat_fast is a function for creating "plot_data" 
# necessary for the viral entities saturation plot
# plot_data is a data frame with X number of rows, where X is 
# the number of samples and 3 columns: Samples (1 - X), 
# Mean_Detected_ent, standard deviation
# input: count_table (rows - entities, columns - samples, row.names - entities names)
# output: list with plot data, summarized permutation data, stats
saturation_stat_fast <- function(count_table, n_permutations) {
  
  bin_mat <- count_table > 0
  
  # 2) Helper to get a single richness curve
  richness_curve <- function(bm) {
    cum_pres <- t(apply(bm, 1, cumsum))
    colSums(cum_pres > 0)
  }
  
  # 3) Compute original
  original <- richness_curve(bin_mat)
  
  # 4) Compute all perms at once
  permuted_matrix <- replicate(
    n_permutations,
    richness_curve(bin_mat[, sample(ncol(bin_mat))])
  )
  
  # 5) Summarize
  mean_votus <- rowMeans(permuted_matrix)
  sd_votus   <- apply(permuted_matrix, 1, sd)
  
  plot_data <- data.frame(
    Samples             = seq_along(original),
    Mean_Detected_vOTUs = mean_votus,
    SD                  = sd_votus
  )
  
  # 6) Marginal‐gain stats on the flattened diffs
  marginal_gains <- apply(permuted_matrix, 2, diff)      # matrix: (n_samples-1) × n_permutations
  all_gains      <- as.vector(marginal_gains)           # flatten to one vector
  
  stat_tab <- data.frame(
    Metric = c(names(summary(all_gains)), "SD"),
    Value  = c(as.numeric(summary(all_gains)), sd(all_gains))
  )
  
  list(
    original = data.frame(Samples = seq_along(original),
                          Detected_vOTUs = original),
    permuted = plot_data,
    stat     = stat_tab
  )
}

# gets the summaries of different curve fittings
podgonian <- function(DF, x, y){
  
  logger <- list()
  
  F1 <- as.formula(paste0("log(", y, ") ~ log(", x, ")"))
  
  fit_power <- lm(F1, data = DF)
  logger[["lm_fit_power"]] <- summary(fit_power)
  
  b <- unname(fit_power$coefficients[2])
  median_index <- round(median(DF[,x]))
  a <- sqrt(DF[median_index,y])
  
  # power function
  F2 <- as.formula(paste(y, " ~ a * ", x, "^b"))
  
  fit_root <- nls(F2, 
                  data = DF, 
                  start = list(a=a,b=b))
  logger[["fit_root"]] <- summary(fit_root)
  
  root_params <- coef(fit_root)
  MDR <- root_params["a"] * root_params["b"] * (nrow(DF)^(root_params["b"]-1))
  logger[["MDR"]] <- MDR
  
  # plateau:
  F3 <- as.formula(paste(y, " ~ ", x))
  fit_plat <- drc::drm(F3, # drc is called here in the isolated way, otherwise it masks "select" from dplyr
                       data = DF, 
                       fct = drc::L.4())
  
  logger[["fit_plat"]] <- summary(fit_plat)
  
  logger[["plat_vs_root"]] <- AIC(fit_plat, fit_root)
  
  # asymptotic model
  F4 <- as.formula( paste(y, " ~ SSasymp(", x, ", Asym, R0, lrc)"))
  
  fit_asym <- nls(F4, 
                  data = DF)
  
  logger[["fit_asym"]] <- summary(fit_asym)
  
  logger[["asym_vs_root"]] <- AIC(fit_asym, fit_root)
  logger[["asym_vs_plat"]] <- AIC(fit_asym, fit_plat)
  
  
  params <- coef(fit_asym)
  Asym <- params["Asym"]
  R0   <- params["R0"]
  lrc  <- params["lrc"]
  
  rate <- exp(lrc)
  # first derivative
  slope <- (Asym - R0) * rate * exp(-rate * nrow(DF))
  
  logger[["slope"]] <- slope
  
  return(logger)
  
  
}
#############################################################
# 1. Loading libraries
#############################################################
library(dplyr)
library(purrr)
library(tidyverse)
#############################################################
# 2. Load Input Data
#############################################################
ETOF <- read.table('06.CLEAN_DATA/VLP_MGS_ETOF_full_rep.txt', sep='\t', header = T)

ETOF_only <- ETOF %>%
  filter(grepl('NEXT_', ETOF$New_CID)) %>%
  mutate(sample = gsub('_.*', '', Original_CID)) %>%
  mutate( method = ifelse(grepl('NEXT_V', New_CID), 'VLP', 'MGS') )

ETOF_vOTUr <- read.table('06.CLEAN_DATA/02.FINAL/Working_ETOF_120997vOTUr_ab3kbp_in_2200_VLP_MGS.txt', sep='\t', header=T)

smeta <- read.delim('06.CLEAN_DATA/02.FINAL/Chiliadal_meta_VLP_MGS_matched_v05_suppl_w_virmetrics.txt', sep='\t', header=T)
smeta$Timepoint_new <- factor(smeta$Timepoint_new, levels=c("M1", "M3", "M6", "M12", "Mother"), ordered = T)

VLP <- read.table('06.CLEAN_DATA/02.FINAL/VLP_only_RPKM_table_VLP_MGS_dec99ANI_ALL_CS_ab3kbp_1110_samples.txt', sep='\t', header=T)

MGS <- read.table('06.CLEAN_DATA/02.FINAL/MGS_only_RPKM_table_VLP_MGS_dec99ANI_ALL_CS_ab3kbp_1110_samples.txt', sep='\t', header=T)

vOTU_clustering <- read.table('06.CLEAN_DATA/NEXT_viral_clusters_MGS_VLP_long_format.txt', sep='\t', header=T)

#############################################################
# 3.1 Analysis: redundant & hq genomes per method
#############################################################

genomic_stat <- ETOF_only %>%
  group_by(method) %>%
  summarise(N_total = n()) %>%
  left_join(ETOF_only %>%
              filter(miuvig_quality == "High-quality") %>%
              group_by(method) %>%
              summarise(N_HQ = n()) ) %>%
  left_join(ETOF_only %>%
              filter(miuvig_quality == "High-quality" & POST_CHV_length >= 10000) %>%
              group_by(method) %>%
              summarise(N_HQ10 = n())) %>%
  pivot_longer(cols = !method) %>%
  pivot_wider(names_from = method,
              values_from = value) %>%
  rename(feature = name) %>%
  mutate(across(where(is.numeric), as.character ))

sample_stat <- ETOF_only %>%
  group_by(sample) %>%
  summarise(N_discovered = n()) %>%
  left_join(ETOF_only %>%
              filter(miuvig_quality == "High-quality") %>%
              group_by(sample) %>%
              summarise(N_discovered_HQ = n())) %>%
  left_join(ETOF_only %>%
              filter(miuvig_quality == "High-quality" & POST_CHV_length >= 10000) %>%
              group_by(sample) %>%
              summarise(N_discovered_HQ10 = n()))

smeta <- smeta %>%
  left_join(sample_stat, by = c("Sequencing_ID" = "sample"))

disc_staton <- map_dfr(c('MGS', 'VLP'), function(metavirome){
  
  lengths <- ETOF_only %>%
    filter(miuvig_quality == "High-quality" & method == metavirome) %>%
    pull(POST_CHV_length) %>%
    get_iqr()
  
  lengths_df <- data.frame(median = lengths,
                           method = metavirome,
                           feature = "Genome_length")
  
  metrics_df <- map_dfr(c("N_discovered", "N_discovered_HQ", "N_discovered_HQ10"), function(metric) {
    
    med_and_iqr <- smeta %>%
    filter(seq_type == metavirome) %>%
      pull(metric) %>%
      get_iqr()
    
    data.frame(median = med_and_iqr,
               method = metavirome,
               feature = metric)
    
  })
  
  rbind(lengths_df, metrics_df)
})

disc_staton <- genomic_stat %>%
  bind_rows(disc_staton %>%
              pivot_wider(names_from = method, 
                          values_from = median) %>%
              mutate(feature = if_else(feature != "Genome_length", paste0(feature, "_per_sample"), feature)))

writexl::write_xlsx(disc_staton, '07.RESULTS/Virus_discovery_summary_stat.xlsx')
# tidying up
rm(list=c("ETOF_only", "genomic_stat", "sample_stat"))

#############################################################
# 3.2 Analysis: N detected vOTUs saturation curves
#############################################################
# since there is a lot of interdependency (paired samples, 
# longitudinal samples) no sense to make saturation over 2,220
# -> merging to holovirome table

##### creating holovirome binary table
colnames(VLP) <- smeta$Universal_ID[match(colnames(VLP), smeta$Sequencing_ID)]
colnames(MGS) <- smeta$Universal_ID[match(colnames(MGS), smeta$Sequencing_ID)]

## ordering
MGS <- MGS[, colnames(VLP)]  # if needed

all_votus <- ETOF_vOTUr$New_CID

# initiate matrix
m1 <- matrix(0,
             nrow = length(all_votus),
             ncol = ncol(VLP),
             dimnames = list(all_votus, colnames(VLP)))
m2 <- m1

# populating matrices
m1[rownames(VLP), ] <- (as.matrix(VLP) > 0) * 1
m2[rownames(MGS), ] <- (as.matrix(MGS) > 0) * 1

## unified presence/absence
unified <- ((m1 + m2) > 0) * 1

holo <- as.data.frame(unified)

holo_saver <- holo

write.table(holo_saver, "06.CLEAN_DATA/Intermediate/Holovirome_RPKM_1110samples_120997vOTUs.txt", sep='\t', quote=F)

rm(list = c("m1", "m2", "all_votus", "unified", "holo_saver"))

gen_holo <- holo %>%
  rownames_to_column("New_CID") %>%
  left_join(ETOF_vOTUr %>% select(New_CID, Genus_OTUr)) %>% # library(drc) masks "select", wtf!  
  group_by(Genus_OTUr) %>%
  summarise(across(where(is.numeric), sum)) %>% # 27,017 genera detected
  column_to_rownames("Genus_OTUr")

fam_holo <- holo %>%
  rownames_to_column("New_CID") %>%
  left_join(ETOF_vOTUr %>% select(New_CID, Family_OTUr)) %>% # library(drc) masks "select", wtf!  
  group_by(Family_OTUr) %>%
  summarise(across(where(is.numeric), sum)) %>% # 7,596 families detected
  column_to_rownames("Family_OTUr")


start.time <- Sys.time()
all_vOTU_cumulative <- saturation_stat_fast(holo, 1000)
all_gOTU_cumulative <- saturation_stat_fast(gen_holo, 1000)
all_fOTU_cumulative <- saturation_stat_fast(fam_holo, 1000)
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken # Time difference of 5.643255 hours

all_cumulative <- all_vOTU_cumulative[["permuted"]] %>%
  mutate(Rank = "Species") %>%
  bind_rows(all_gOTU_cumulative[["permuted"]] %>% mutate(Rank = "Genus")) %>%
  bind_rows(all_fOTU_cumulative[["permuted"]] %>% mutate(Rank = "Family"))

rm(list = c("gen_holo", "fam_holo", "all_vOTU_cumulative", "all_gOTU_cumulative", "all_fOTU_cumulative"))

all_satu <- ggplot(all_cumulative, aes(x = Samples, y = Mean_Detected_vOTUs, color=Rank)) +
  geom_line() +
  geom_point() +
  geom_errorbar(aes(ymin = Mean_Detected_vOTUs - SD, ymax = Mean_Detected_vOTUs + SD), width = 0.1) +
  labs(x = "Number of fecal samples", y = "Number of OTUs") +
  scale_color_manual(values = c("#7C3E66", "#A5BECC", "#243A73")) +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave('05.PLOTS/05.VLP_MGS/all_OTUs_detected_holovirome.png',
       all_satu,  "png", width=10, height=8, units="cm", dpi = 300)

ggsave('05.PLOTS/05.VLP_MGS/all_OTUs_detected_holovirome.pdf',
       all_satu,  "pdf", width=10, height=8, units="cm", dpi = 300)

# looking for approximation:
# Species curve looks like something * sqrt()

all_cumulative %>%
  filter(Rank == "Species") %>%
  ggplot(aes(log(Samples), log(Mean_Detected_vOTUs))) +
  geom_point() +
  theme_minimal() +
  labs(y="log(mean detected species)") # ehm kind of straight line?

sp_stat <- podgonian(all_cumulative[all_cumulative$Rank == "Species",],
                     "Samples",
                     "Mean_Detected_vOTUs") # best is root square

gn_stat <- podgonian(all_cumulative[all_cumulative$Rank == "Genus",],
                     "Samples",
                     "Mean_Detected_vOTUs") # best is root square, but gain is only 8 or 2 genera ps -> plateau

fm_stat <- podgonian(all_cumulative[all_cumulative$Rank == "Family",],
                     "Samples",
                     "Mean_Detected_vOTUs") # best is root square, but gain is only 2 or (0.5) ps -> def plateau

# stat for NEXT vs external DB-clustered seqs
holo_next <- holo %>%
  filter(row.names(.) %in% ETOF_vOTUr$New_CID[ETOF_vOTUr$vOTU_novelty == "novel"])

holo_db <- holo %>%
  filter(row.names(.) %in% ETOF_vOTUr$New_CID[ETOF_vOTUr$vOTU_novelty == "described"])

start.time <- Sys.time()
NEXT_vOTU_cumulative <- saturation_stat_fast(holo_next, 1000)

DB_vOTU_cumulative <- saturation_stat_fast(holo_db, 1000)
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken


cumulative_by_source <- NEXT_vOTU_cumulative[["permuted"]] %>%
  mutate(Source = "Study-derived") %>%
  bind_rows(DB_vOTU_cumulative[["permuted"]] %>% mutate(Source = "DB-matched"))

satu_by_source <- ggplot(cumulative_by_source, aes(x = Samples, y = Mean_Detected_vOTUs, color=Source)) +
  geom_line() +
  geom_point() +
  geom_errorbar(aes(ymin = Mean_Detected_vOTUs - SD, ymax = Mean_Detected_vOTUs + SD), width = 0.1) +
  labs(x = "Number of fecal samples", y = "Number of vOTUs") +
  scale_color_manual(values = c(MetBrewer::met.brewer("Derain")[3],
                                MetBrewer::met.brewer("Derain")[7])) +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave('05.PLOTS/05.VLP_MGS/all_vOTUs_detected_holovirome_by_source.png',
       satu_by_source,  "png", width=10, height=8, units="cm", dpi = 300)

ggsave('05.PLOTS/05.VLP_MGS/all_vOTUs_detected_holovirome_by_source.pdf',
       satu_by_source,  "pdf", width=10, height=8, units="cm", dpi = 300)

# stat:
next_stat <- podgonian(cumulative_by_source[cumulative_by_source$Source == "Study-derived",],
                       "Samples",
                       "Mean_Detected_vOTUs") # best is root square

db_stat <- podgonian(cumulative_by_source[cumulative_by_source$Source == "DB-matched",],
                     "Samples",
                     "Mean_Detected_vOTUs") # best is sigmoid, gain is 5 vOTUs per new sample, plateau
#############################################################
# 3.3 Analysis: descriptive stat NEXT virome catalog
#############################################################

table(ETOF_vOTUr$miuvig_quality) # 15,969 complete or nearly complete (>90%)
min(ETOF_vOTUr[ETOF_vOTUr$miuvig_quality=="High-quality",]$POST_CHV_length) #3,000
max(ETOF_vOTUr[ETOF_vOTUr$miuvig_quality=="High-quality",]$POST_CHV_length) #396,781 
# % genomes, % temperate etc was derived from summary stat of ETOF

virnov <- ETOF_vOTUr  %>%
  group_by(vOTU_novelty) %>%
  summarise(sum=n(), .groups = "drop") %>%
  mutate(perc = round(sum/nrow(ETOF_vOTUr) * 100, 1))

virclass <- ETOF_vOTUr  %>%
  separate(
    col    = tax_ictv_aai,
    into = c('Life', 'Realm', 'Kingdom', 'Phylum',  'Class', 'Order', 'Family',  'Genus',  'Species', 'Strain'),
    sep    = ";",
    fill   = "right",               
    remove = FALSE                  
  ) %>%
  group_by(Class) %>%
  summarise(sum=n()) %>%
  mutate(perc = round(sum / nrow(ETOF_vOTUr) * 100, 1))

virhost_simple <- ETOF_vOTUr  %>%
  group_by(Host_simple) %>%
  summarise(sum = n(), .groups = "drop") %>%
  mutate(perc = round(sum/nrow(ETOF_vOTUr) * 100, 1))

virhost <- ETOF_vOTUr  %>%
  filter(Host_simple!="Eukaryote") %>%
  select(Host_taxonomy) %>%
  separate(
    col    = Host_taxonomy,
    into = c('Domain', 'Phylum', 'Class', 'Order', 'Family',  'Genus'),
    sep    = ";",
    fill   = "right",               
    remove = FALSE                  
  ) %>%
  group_by(Genus) %>%
  summarise(sum=n(), .groups = "drop") %>%
  mutate(perc = round(sum / nrow(ETOF_vOTUr) * 100, 1))

#############################################################
# 3.4 Analysis: novels' descriptive
#############################################################
# define novel genera and families:
ETOF_vOTUr <- ETOF_vOTUr %>%
  group_by(Genus_OTUr) %>%
  mutate(Genus_OTU_novelty = if_else(all(vOTU_novelty == "novel"), "novel", "described")) %>%
  ungroup() %>%
  group_by(Family_OTUr) %>%
  mutate(Family_OTU_novelty = if_else(all(vOTU_novelty == "novel"), "novel", "described")) %>%
  ungroup()

# count novel genera and families:
novel_votus <- ETOF_vOTUr %>% 
  filter(miuvig_quality == "High-quality", vOTU_novelty == "novel")

# stats on the number of novel genera and novel families
stats_nov_higher <- list(
  n_vOTUs_in_novel_genera = novel_votus %>% filter(Genus_OTU_novelty == "novel") %>% nrow(),
  n_novel_genera         = novel_votus %>% filter(Genus_OTU_novelty == "novel") %>% pull(Genus_OTUr) %>% n_distinct(),
  n_novel_families       = novel_votus %>% filter(Family_OTU_novelty == "novel") %>% pull(Family_OTUr) %>% n_distinct(),
  n_genera_in_novel_fams = novel_votus %>% filter(Family_OTU_novelty == "novel") %>% pull(Genus_OTUr) %>% n_distinct()
)

novel_by_genome <- novel_votus %>%
  group_by(genome) %>%
  summarise(sum = n(), .groups = "drop") %>%
  mutate(perc = round(sum/nrow(novel_votus) * 100, 2))

novel_by_Class <- novel_votus %>%
  separate(
    col    = tax_ictv_aai,
    into = c('Life', 'Realm', 'Kingdom', 'Phylum',  'Class', 'Order', 'Family',  'Genus',  'Species', 'Strain'),
    sep    = ";",
    fill   = "right",               
    remove = FALSE                  
  ) %>%
  group_by(Class, genome) %>%
  summarise(sum=n(), .groups = "drop") %>%
  mutate(perc = round(sum / nrow(novel_votus) * 100, 2))

novel_by_Order <- novel_votus %>%
  separate(
    col    = tax_ictv_aai,
    into = c('Life', 'Realm', 'Kingdom', 'Phylum',  'Class', 'Order', 'Family',  'Genus',  'Species', 'Strain'),
    sep    = ";",
    fill   = "right",               
    remove = FALSE                  
  ) %>%
  group_by(Order, genome) %>%
  summarise(sum=n(), .groups = "drop") %>%
  mutate(perc = round(sum / nrow(novel_votus) * 100, 2))

novel_virhost_simple <- novel_votus  %>%
  group_by(Host_simple) %>%
  summarise(sum=n(), .groups = "drop") %>%
  mutate(perc = round(sum / nrow(novel_votus) * 100, 1))

novel_virhost <- novel_votus  %>%
  filter(Host_simple!="Eukaryote") %>%
  select(Host_taxonomy) %>%
  separate(
    col    = Host_taxonomy,
    into = c('Domain', 'Phylum', 'Class', 'Order', 'Family',  'Genus'),
    sep    = ";",
    fill   = "right",               
    remove = FALSE                  
  ) %>%
  group_by(Family) %>%
  summarise(sum=n(), .groups = "drop") %>%
  mutate(perc = round(sum / nrow(novel_votus) * 100, 1))

prev_in_holo <- holo %>%
  filter(row.names(holo) %in% novel_votus$New_CID) %>%
  mutate(total_presence = rowSums(select(., where(is.numeric)) > 0))

# stats on prevalence of novel viruses in our dataset (holovirome)
stats_nov_prevalence <- list(
  novel_singletons = sum(prev_in_holo$total_presence == 1),
  perc_novel_singletons = sum(prev_in_holo$total_presence == 1)/nrow(novel_votus),
  novel_perv_5_perc       = sum(prev_in_holo$total_presence >= 0.05*ncol(holo)),
  perc_novel_prev_5_perc = sum(prev_in_holo$total_presence >= 0.05*ncol(holo))/nrow(novel_votus)*100
)

# external DB ssDNA & RNA content:
ETOF_DB <- ETOF %>%
  filter(!grepl('NEXT_V|NEXT_M', New_CID) & miuvig_quality == "High-quality" & New_CID %in% vOTU_clustering$Representative) %>% # 60k vOTUs
  separate(
    col    = taxonomy,
    into = c('Life', 'Realm', 'Kingdom', 'Phylum',  'Class', 'Order', 'Family'), # nothing below Family is available
    sep    = ";",
    fill   = "right",               
    remove = FALSE                  
  ) %>% # only 138 is total Unclassified from 60k genomes -> looks OK for assessing the representation of higher tax ranks
  group_by(Realm) %>% # using Realm here as a proxy for genome; RERUN W ICTV AS RUSHECHKA
  summarise(sum=n(), .groups = "drop") %>%
  mutate(perc = sum/sum(sum) * 100) %>%
  mutate(genome = case_when(Realm %in% c('Monodnaviria', 'Anelloviridae') ~ "ssDNA",
                            Realm == "Riboviria" ~ "RNA",
                            .default = "dsDNA")) %>%
  group_by(genome) %>%
  summarise(sum = sum(sum), .groups = "drop")

contrib_to_db <- novel_by_genome %>%
  filter(genome!="Unclassified") %>%
  select(genome, sum) %>%
  mutate(source = "NEXT") %>%
  bind_rows(ETOF_DB %>% select(genome, sum) %>% mutate(source = "DB") ) %>%
  mutate(genome = factor(genome, levels = c('dsDNA', 'ssDNA', 'RNA'), ordered = T),
         source = factor(source, levels = c('NEXT', 'DB'), ordered = T)) %>%
  ggplot(aes(x = genome, y = sum, fill = source)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +
  geom_text(aes(label = sum), position = position_stack(vjust = 0.5), size = 3, color = "black") +
  scale_fill_manual(values = c("NEXT" = "#E64B35FF", "DB" = "gray70"), labels = c("Novel", "DB-deposited")) + 
  theme_bw() +
  labs(y = "N of (near-)complete vOTUs", x = "Genome Type", fill = "Data Source") +
  facet_wrap(~genome, scales = "free") +
  theme(strip.background = element_rect(NA))


ggsave('05.PLOTS/04.VIRAL_DB/Novel_contribution_to_extDB.png',
       contrib_to_db,  "png", width=12, height=12, units="cm", dpi = 300)

ggsave('05.PLOTS/04.VIRAL_DB/Novel_contribution_to_extDB.pdf',
       contrib_to_db,  "pdf", width=12, height=12, units="cm", dpi = 300)

#############################################################
# 4. OUTPUT
#############################################################
# graphs section-wise
