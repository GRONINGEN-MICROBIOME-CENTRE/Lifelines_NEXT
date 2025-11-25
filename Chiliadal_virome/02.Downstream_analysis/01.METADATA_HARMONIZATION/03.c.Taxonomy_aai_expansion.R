setwd('~/Desktop/Projects_2021/NEXT_virome/09.DATA_ANALYSIS/')

#############################################################
# Expanding consensus taxonomy according to AAI clustering
#############################################################

#############################################################
# 0. Used files source
#############################################################

#############################################################
# 1. Functions
#############################################################
ranks <- c("Strain","Species","Genus","Family","Order","Class","Phylum","Kingdom","Realm", "Life")

# adds singleton families for vOTUs lost during AAI due to 
# the blast & clustering logic
add_lost_taxa <- function(rank_clustering, # table (data.frame) in the long format that can be aggregated according to the rank column
                          cluster_size, # table in the wide format, indicates the cluster sizes for every unique rank
                          all_subtaxa, # all vOTUs that should have been there, character vector
                          rank # character, (Genus or Family in this case)
                          ) {
  
  n_taxa <- nrow(cluster_size)
  
  lost_taxa_vec <- setdiff(all_subtaxa, rank_clustering$Cluster_member)
  
  N_lost_taxa <- length(lost_taxa_vec)
  
  if (N_lost_taxa > 0) {
    
    lost_taxa <- data.frame(Cluster_member = lost_taxa_vec, 
                            stringsAsFactors = FALSE)
    
    lost_taxa$Representative <- paste0(
      rank, '_',
      (n_taxa + 1) : (n_taxa + N_lost_taxa)
    )
    
    rank_clustering_UPD <- rbind(rank_clustering, lost_taxa)
  } else {
    
    rank_clustering_UPD <- rank_clustering
    
  }
  
  cluster_size_UPD <- as.data.frame(table(rank_clustering_UPD$Representative))
  colnames(cluster_size_UPD) <- colnames(cluster_size)
  
  list(
    rank_clustering_UPD = rank_clustering_UPD,
    cluster_size_UPD    = cluster_size_UPD
  )
  
}

# helper, identifies the first rank with classified values (not to use outside the main function!)
rank_seeker <- function(df, lineage) {
  
  for (rnk in lineage) {
    
    first_classified <- any(df[,rnk] != "Unclassified")
    
    if (first_classified) { break }
    
  }
  
  return(rnk)
  
}

# helper, chooses the best propagator based on the majority of assignments, in case of ties chooses
# the assignment with the longest representative
# further ties are resolved by priority
propagator_chooser <- function(df, prop_rank, candidates, lineage) {
  
  df_filtered <- df %>%
    filter(!!sym(prop_rank) %in% candidates[[prop_rank]]) %>%
    group_by(across(all_of(lineage))) %>%
    summarise(
      n = n(),                       # number of rows per lineage combo
      POST_CHV_length = max(POST_CHV_length),  # longest length in that group
      .groups = "drop"
    ) %>%
    arrange(desc(n), desc(POST_CHV_length)) %>%
    slice_head(n=1)
  
  return(df_filtered[,all_of(lineage)])
}

# helper, runs rank-wise comparison of the given taxonomy and the taxonomy of the propagator
# only compares them at given lineage
conflict_seeker <- function(s1, propagator, lineage) {
  
  conflict <- FALSE
  
  for (rnk in rev(lineage)) {
    
    if ( s1[rnk] == propagator[rnk] | s1[rnk] == "Unclassified" ) {
      
      next
      
    } else {
      
      conflict <- TRUE
      #print(paste0("conflict at ", rnk))
      
    }
    
    if (conflict) {
      break
    }
    
  }
  
  return(conflict)
}

# propagates taxonomy
# it might overdrive "Unassigned", but it does not change the taxonomy composition / distribution
propagate_taxonomy_patched <- function(df, # df
                                       cluster_col, # name of the AAI cluster column (tested and written for Genus and Family only)
                                       taxon_col # name of the taxon where to start propagator seeker
                                       ) { 
  
  n <- grep(taxon_col, ranks)
  rank_cols <- if (n < length(ranks)) ranks[(n + 1):length(ranks)] else character(0)
  lineage <- c(rev(rank_cols), taxon_col)
  
  for ( i in unique(df[[cluster_col]]) ) {
    
    OTU <- df %>%
      filter(!!sym(cluster_col) == i)
    
    most_freq_rows <- OTU[,lineage] %>%
      count(across(everything())) %>%   # count each unique combination of columns
      arrange(desc(n))
    
    if ( nrow(most_freq_rows) > 1 ) { # if there is something to propagate or somewhere to propagate to
      
      prop_rank <- rank_seeker(OTU, rev(lineage)) # rank from which propagation should start
      
      candidates <- most_freq_rows %>% # at this rank, is there a single candidate or many?
        filter(!!sym(prop_rank) != "Unclassified")
      
      if ( nrow(candidates) == 1) { # if only one -> it is the propagator
        
        propagator <- candidates[,lineage]
        
      } else { # if more than 1, choose propagator based on the majority of assignments, in case of ties: first length, then priority
        
        propagator <- propagator_chooser(OTU, prop_rank, candidates,lineage)
        
      }
      
      # preparing propagator (not to overwrite things, just in case):
      k <- grep(prop_rank, ranks)
      rank_cols2 <- if (k < length(ranks)) rev(ranks[(k + 1):length(ranks)]) else character(0)
      prop_lineage <- c(rank_cols2, prop_rank) # select the propagator lineage based on the prop rank, not the entire lineage (i.e., with genus)
      downstream <- if (k > 1) rev(ranks[(k - 1):2]) else character(0)
      
      for (votur in OTU$New_CID) {
        
        current_taxonomy <- OTU[OTU$New_CID == votur, lineage]
        
        if (conflict_seeker(current_taxonomy, propagator, lineage) ) { # if there is a conflict -> skip this line
          
          next
          
        } else { # else - propagate
          
          if ( all(df[df$New_CID == votur, downstream] == "Unclassified") ) { # if no downstream conflicts
            df[df$New_CID == votur, prop_lineage] <- propagator[prop_lineage]
          }
          
        }
        
      }
      
      
    }
    else next # keep as it is
  }
  return(df)
}

# reused & modified from 3b
# finds cases when the rank of interest has multiple assigned upstream ranks
check_fd <- function(df, 
                     lower, 
                     higher, 
                     drop_na_higher = TRUE,
                     ignore_lower_values = c("Unclassified"),
                     ignore_case = TRUE) {
  stopifnot(lower %in% names(df), higher %in% names(df))
  
  # build comparator for ignored lower values
  norm_case <- function(x) if (ignore_case) str_to_lower(x) else x
  ignore_set <- norm_case(ignore_lower_values)
  
  tmp <- df %>%
    select(all_of(c(lower, higher))) %>%
    mutate(across(everything(), ~ .x %>% str_squish() %>% na_if("")))
  
  if (drop_na_higher) {
    tmp <- tmp %>% filter(!is.na(.data[[higher]]))
  }
  
  # Count distinct higher per lower
  summary <- tmp %>%
    filter(.data[[higher]] != ignore_lower_values) %>%
    group_by(.data[[lower]]) %>%
    summarize(
      n = n(),
      n_distinct_higher = n_distinct(.data[[higher]], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    rename(!!lower := 1L)
  
  # Filter offenders (exclude NA and ignored lower values)
  offenders_keys <- summary %>%
    mutate(.lower_norm = norm_case(!!sym(lower))) %>%
    filter(
      n_distinct_higher > 1,
      !is.na(!!sym(lower)),
      !.lower_norm %in% ignore_set
    ) %>%
    select(all_of(lower))
  
  offenders <- tmp %>%
    semi_join(offenders_keys, by = setNames(lower, lower)) %>%
    count(.data[[lower]], .data[[higher]], name = "assignments") %>%
    arrange(.data[[lower]], desc(assignments)) %>%
    rename(
      !!lower := 1L,
      !!higher := 2L
    )
  
  list(
    lower = lower,
    higher = higher,
    summary = summary,      # includes Unclassified rows for visibility
    offenders = offenders,  # excludes Unclassified (and other ignored) conflicts
    ok = nrow(offenders) == 0
  )
}
#############################################################
# 1. Loading libraries
#############################################################
library(dplyr)
library(tidyverse)
#############################################################
# 2. Load Input Data
#############################################################
taxa <- read.table('06.CLEAN_DATA/Consensus_taxonomy_121062_ICTV_curated.txt', sep='\t', header=T)

ETOF_vOTUr <- read.table('06.CLEAN_DATA/02.FINAL/Basic_ETOF_121062vOTUr_ab3kbp_in_2200_VLP_MGS.txt', sep='\t', header=T)

genus_clusters <- read.table('06.CLEAN_DATA/genus_clusters_121062vOTUr_labeled_long_format.txt', sep = '\t', header=T)
genus_size <- read.table('06.CLEAN_DATA/genus_clusters_121062vOTUr_labeled_size.txt', sep='\t', header=T)

family_clusters <- read.table('06.CLEAN_DATA/family_clusters_121062vOTUr_labeled_long_format.txt', sep='\t', header = T)
family_size <- read.table('06.CLEAN_DATA/family_clusters_121062vOTUr_labeled_size.txt', sep='\t', header=T)

#############################################################
# 3.1 Analysis: adding genera and families lost by AAI
#############################################################
# genera
res <- add_lost_taxa(genus_clusters, genus_size, ETOF_vOTUr$New_CID, 'Genus')

genus_clusters <- res$rank_clustering_UPD
genus_size <- res$cluster_size_UPD

# genera tables are ready to be saved

# families
res <- add_lost_taxa(family_clusters, family_size, ETOF_vOTUr$New_CID, 'Family')

family_clusters <- res$rank_clustering_UPD
family_size <- res$cluster_size_UPD

# family tables need further curation

rm(res)
#############################################################
# 3.2 Analysis: solving conflict between AAI genera and families
#############################################################

df <- taxa %>%
  left_join(genus_clusters, by=c("New_CID" = "Cluster_member")) %>%
  rename(Genus_OTUr = Representative) %>%
  left_join(family_clusters, by=c("New_CID" = "Cluster_member")) %>%
  rename(Family_OTUr = Representative)

# check how many families there are per single genus
genus_family_check <- df %>%
  distinct(Genus_OTUr, Family_OTUr) %>%
  count(Genus_OTUr) %>%
  filter(n > 1) # -> in the current release, 16.9% of AAI genera have at least 2 Family assignments

# choosing a major family (random/priority tiebreaker):
major_family <- df %>%
  group_by(Genus_OTUr, Family_OTUr) %>%
  tally(name = "count") %>%
  slice_max(order_by = count, n = 1, with_ties = FALSE) %>%
  ungroup()

# fixed (no more than 1 family per genus):
df_fixed <- df %>%
  select(-Family_OTUr) %>%            
  left_join(major_family, by = "Genus_OTUr") %>%
  select(-count) 

family_clusters_UPD <- df_fixed[,c("Family_OTUr", "New_CID")]
colnames(family_clusters_UPD) <- colnames(family_clusters)

family_size_UPD <- as.data.frame(table(family_clusters_UPD$Representative))
colnames(family_size_UPD) <- colnames(family_size)
# families are ready to be saved

rm(list=c("family_size", "family_clusters", "genus_family_check",
          "family_clusters_UPD", "family_size_UPD",
          'genus_clusters', 'genus_size', 'major_family', "df"))
#############################################################
# 3.3 Analysis: checking AAI ranks vs (higher) taxa ranks conflicts
#############################################################
pairs_to_check <- tidyr::expand_grid(
  lower  = c('Genus', 'Family'),
  higher = ranks
) %>%
  filter(match(lower, ranks) <= match(higher, ranks))

pairs_to_check$lower <- gsub('Genus', 'Genus_OTUr', pairs_to_check$lower)
pairs_to_check$lower <- gsub('Family', 'Family_OTUr', pairs_to_check$lower)

fd_results <- pairs_to_check %>%
  pmap(function(lower, higher) check_fd(df_fixed, lower, higher))

# reporting offender checker results:
report <- map_df(fd_results, function(res) {
  tibble(
    rule = paste0(res$lower, " → ", res$higher),
    offenders = n_distinct(res$offenders[[res$lower]] %||% character()),
    ok = res$ok
  )
})

print(report, n=100)

# To inspect specific violations, e.g. genus_OTUr -> genus
get_res <- function(results, lower, higher) {
  keep(results, ~ .x$lower == lower && .x$higher == higher)[[1]]
}

res_genus_family <- get_res(fd_results, "Genus_OTUr", "Genus")
res_genus_family$offenders 
# after some checks, I can leave with these conflicts (higher ranks as well)
# taxonomy will be propagated upwards for Unclassified ranks only & if there
# are no upper rank conflicts. Offenders will not be overwritten.
rm(list= c("pairs_to_check", "fd_results", "res_genus_family"))
#############################################################
# 3.4 Analysis: check if Unclassified can be assigned by majority
#############################################################
top12_wide <- df_fixed %>%

  filter(Genus != "Unclassified") %>%
  
  group_by(Genus_OTUr, Genus) %>%
  summarise(n = n(), .groups = "drop_last") %>%

  mutate(
    total_in_rep = sum(n),
    prop         = n / total_in_rep
  ) %>%
  ungroup() %>%

  group_by(Genus_OTUr) %>%
  filter(n_distinct(Genus) > 1) %>%

  arrange(Genus_OTUr, desc(n)) %>%
  mutate(rank = row_number()) %>%
  
  # keep only top 2
  filter(rank <= 2) %>%
  ungroup() %>%
  
  # go wide
  select(Genus_OTUr, rank, Genus, n, prop) %>%
  pivot_wider(
    names_from  = rank,
    values_from = c(Genus, n, prop),
    names_glue  = "top{rank}_{.value}"
  )

View(top12_wide[top12_wide$top1_n == top12_wide$top2_n,]) # there are cases with ties, tie-breaker -> length of the genome(s)
# if genome length tie-breaker does not resolve it, it will be a priority resolution (or I'll think about it later)
rm(top12_wide)
#############################################################
# 3.5 Analysis: expanding taxonomy assignment based on genus & family
#############################################################
# @ AAI genus level:
df_updated <- propagate_taxonomy_patched(
  df_fixed,
  cluster_col = "Genus_OTUr",
  taxon_col   = "Genus"
)

# @ AAI family level:
df_final <- propagate_taxonomy_patched(
  df_updated,
  cluster_col = "Family_OTUr",
  taxon_col   = "Family"
)

# crassvirales (I did it this way only because I first checked they are all Caudoviricetes with unclassified lower ranks OR
# already Crassvirales):
df_final$Order[grepl('Guerin|Yutin|NCBI_CrAss|NL_crAss', df_final$New_CID)] <- 'Crassvirales'

rm(list = c("df_fixed", "ETOF_vOTUr"))
#############################################################
# 3.6 Analysis: check if there is a difference in N assigned
#############################################################
# compared to the previous calculations, this one is done at the level of a single vOTU, not unique taxa, so numbers are different

taxa_change_number <- as.data.frame(matrix(NA, nrow=10, ncol=4))
colnames(taxa_change_number) <- c('Domain', "perc_classified_before", "perc_classified_after_gen", "perc_classified_after_fam")
taxa_change_number$Domain <- ranks

for (rank in taxa_change_number$Domain) {
  
  taxa_change_number[taxa_change_number$Domain==rank,"perc_classified_before"] <- 
    round(sum(!taxa[,rank] %in% c("Unclassified")) / 121062 * 100, 2)
  
  taxa_change_number[taxa_change_number$Domain==rank,"perc_classified_after_gen"] <- 
    round(sum(!df_updated[,rank] %in% c("Unclassified")) / 121062 * 100, 2)
  
  taxa_change_number[taxa_change_number$Domain==rank,"perc_classified_after_fam"] <- 
    round(sum(!df_final[,rank] %in% c("Unclassified")) / 121062 * 100, 2)
} # -> Family expansion significantly increases the number of assignments, overshot?

#############################################################
# 4. OUTPUT
#############################################################
# updated AAI genus clustering:
write.table(genus_clusters, '06.CLEAN_DATA/02.FINAL/genus_clusters_121062vOTUr_curated_long_format.txt', sep='\t', row.names = F, quote = F)

# updated genus cluster size:
write.table(genus_size, '06.CLEAN_DATA/02.FINAL/genus_clusters_121062vOTUr_curated_size.txt', sep='\t', row.names = F, quote = F)

# updated AAI family clustering:
write.table(family_clusters_UPD, '06.CLEAN_DATA/02.FINAL/family_clusters_121062vOTUr_curated_long_format.txt', sep='\t', row.names = F, quote = F)

# updated family cluster size:
write.table(family_size_UPD, '06.CLEAN_DATA/02.FINAL/family_clusters_121062vOTUr_curated_size.txt', sep='\t', row.names = F, quote = F)

# expanded taxonomy:
write.table(df_final, "06.CLEAN_DATA/02.FINAL/Consensus_taxonomy_121062_ICTV_curated_AAI_expanded.txt", sep='\t', row.names = F, quote = F)

# stat on increase in taxonomy assignment:
write.table(taxa_change_number, "06.CLEAN_DATA/Table_change_in_Unassigned_after_AAI_expansion.txt", sep='\t', row.names = F, quote = F)

# report on AAI ranks vs taxon conflicts:
write.table(report, "06.CLEAN_DATA/Table_N_offenders_per_rank_in_AAI_vs_taxa_conflict.txt", sep='\t', row.names = F, quote = F)
