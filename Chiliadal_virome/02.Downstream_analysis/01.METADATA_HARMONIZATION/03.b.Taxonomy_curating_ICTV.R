setwd('~/Desktop/Projects_2021/NEXT_virome/09.DATA_ANALYSIS/')

#############################################################
# Refining consensus taxonomy according to the latest ICTV
# release
#############################################################

#############################################################
# 0. Used files source
#############################################################

#############################################################
# 1. Functions
#############################################################
ranks <- c("Species","Genus","Family","Order","Class","Phylum","Kingdom","Realm")

rename_and_update_patched <- function(df, 
                                      ictv, 
                                      renames, 
                                      rank_col, 
                                      upstream_cols) {
  rank_col_sym <- sym(rank_col)
  
  # candidates to rename at this rank
  candidates <- df %>%
    filter(
      .data[[rank_col]] != "Unclassified",
      !.data[[rank_col]] %in% ictv[[rank_col]],
      .data[[rank_col]] %in% renames$old_tax_name
    ) %>%
    distinct(.data[[rank_col]]) %>%
    rename(old_tax_name = !!rank_col_sym)
  
  if (nrow(candidates) == 0) return(df)
  
  # old -> new map
  map_tbl <- candidates %>%
    left_join(renames %>% distinct(old_tax_name, new_tax_name), by = "old_tax_name")
  
  # 1) rename the focal rank
  df2 <- df %>%
    left_join(map_tbl, by = setNames("old_tax_name", rank_col)) %>%
    mutate(
      !!rank_col_sym := if_else(!is.na(new_tax_name), new_tax_name, !!rank_col_sym)
    ) %>%
    select(-new_tax_name)
  
  # 2) attach ICTV lookup (keyed by the *new* rank)
  ictv_lookup <- ictv %>%
    #filter(!is.na(.data[[rank_col]])) %>%
    select(all_of(c(rank_col, upstream_cols))) %>%
    distinct()
  
  
  unassignders <- df2 %>%
    filter(!!rank_col_sym == "Unassigned")
  
  df2 <- df2 %>%
    filter(!!rank_col_sym != "Unassigned") %>%
    left_join(ictv_lookup, by = rank_col, suffix = c("", ".ictv"))
  
  # rows that actually got renamed
  renamed_keys <- unique(na.omit(map_tbl$new_tax_name))
  renamed_row  <- df2[[rank_col]] %in% renamed_keys
  
  # 3) overwrite upstream ranks only for renamed rows and non-NA ICTV values
  for (u in upstream_cols) {
    ictv_u <- paste0(u, ".ictv")
    if (ictv_u %in% names(df2)) {
      take_ictv <- renamed_row
      df2[[u]][take_ictv] <- df2[[ictv_u]][take_ictv]
      df2[[ictv_u]] <- NULL
    }
  }
  
  df2 <- bind_rows(df2, unassignders)
  df2
}

# Driver: iterate ranks and call existing rename_and_update() each time
rename_all_ranks <- function(df, 
                             ictv, 
                             ICTV_renaming, 
                             ranks_vec = ranks) {
  out <- df
  for (i in seq_along(ranks_vec)) {
    focal    <- ranks_vec[i]
    upstream <- if (i < length(ranks_vec)) ranks_vec[(i + 1):length(ranks_vec)] else character(0)
    out <- rename_and_update_patched(
      df            = out,
      ictv          = ictv,
      renames       = ICTV_renaming,
      rank_col      = focal,
      upstream_cols = upstream
    )
  }
  out
}

# identifies placeholders introduced by VITAP
# renames them according to the new ICTV release
# identifies ranks named according to the previous ICTV releases,
# changes their taxonomy to ICTV latest release
rename_vitap_placeholders <- function(df,
                                      ictv,
                                      ICTV_renaming,
                                      ranks) {
  
  # --- 1) Old → new taxon map as a named vector ---
  rename_vec <- ICTV_renaming %>%
    distinct(old_tax_name, new_tax_name) %>%
    { setNames(.$new_tax_name, .$old_tax_name) }
  
  # --- 2) Per-rank lookup for ICTV rows (taxon → row index) ---
  ictv_by_rank <- lapply(stats::setNames(ranks, ranks), function(rk) {
    x   <- ictv[[rk]]
    idx <- which(!is.na(x))
    split(idx, x[idx])   # named list: taxon -> row indices in ictv
  })
  
  # --- 3) Main loop over ranks ---
  for (i in seq_along(ranks)) {
    focal <- ranks[i]
    downstream <- if (i > 1) ranks[seq_len(i - 1)] else character(0)
    
    # Find candidate rows just once
    candidates_raw <- df %>%
      filter(
        if (length(downstream)) if_all(all_of(downstream), ~ .x == "Unclassified") else TRUE,
        str_detect(.data[[focal]], "\\[")
      )
    
    if (nrow(candidates_raw) == 0L) next
    
    # Unique placeholders at this focal rank
    candidates <- candidates_raw %>%
      distinct(placeholder = .data[[focal]]) %>%
      mutate(
        taxon = placeholder %>%
          str_remove_all("\\[|\\]") %>%
          str_replace(".*_", ""),        # keep part after last '_'
        taxon = dplyr::if_else(
          taxon %in% names(rename_vec),
          rename_vec[taxon],
          taxon
        )
      )
    
    # If there are no preceding ranks, you cannot match anything upstream
    if (!length(downstream)) next
    
    # Only consider taxa that exist in ICTV at any downstream rank
    taxa_in_ictv <- purrr::map_lgl(
      candidates$taxon,
      ~ any(map_lgl(downstream, function(rk) .x %in% names(ictv_by_rank[[rk]])))
    )
    if (!any(taxa_in_ictv)) next
    
    candidates_ok <- candidates[taxa_in_ictv, , drop = FALSE]
    
    # For each preceding rank, update *all* matching placeholders in one go
    for (matchrank in downstream) {
      
      # taxa that appear at this specific rank in ictv
      taxa_here <- intersect(candidates_ok$taxon, names(ictv_by_rank[[matchrank]]))
      if (!length(taxa_here)) next
      
      # Build a small lookup table: taxon -> ICTV row (one row per taxon)
      map_df <- map_dfr(taxa_here, function(tx) {
        idx <- ictv_by_rank[[matchrank]][[tx]][1]  # if multiple rows, take the first
        out <- ictv[idx, ranks, drop = FALSE]
        out$taxon <- tx
        out
      })
      
      # Join with placeholders (focal value) for this rank
      placeholder_map <- candidates_ok %>%
        inner_join(map_df, by = "taxon") # has columns: placeholder, taxon, ranks...
      
      # Which ranks to update from this focal point upward?
      ranks_to_update <- rev(setdiff(ranks, downstream))
      
      # Match df rows to placeholders
      m <- match(df[[focal]], placeholder_map$placeholder)
      rows <- which(!is.na(m))
      if (!length(rows)) next
      
      df[rows, ranks_to_update] <- placeholder_map[m[rows], ranks_to_update]
    }
  }
  
  #df[is.na(df)] <- "Unclassified"
  
  df
}

# identifies ranks not in ictv, renames them to "Unclassified"
# identifies ranks in ICTV and renames their upstream ranks to
# ICTV-approved ranks
fill_from_ictv_patched <- function(df,
                                   ictv,
                                   ranks) {
  # df:   taxonomy table with entries split by ranks
  # ictv: ictv table split by ranks
  
  for (i in seq_along(ranks)) {
    rnk <- ranks[i]
    
    # ranks "above" this one (or below, depending on order of 'ranks')
    upstream <- if (i < length(ranks)) ranks[(i + 1):length(ranks)] else character(0)
    
    # 1) mark taxa not present in ictv as Unclassified at this rank
    is_known <- df[[rnk]] %in% ictv[[rnk]]
    df[[rnk]][!is_known] <- "Unclassified" # WHAT WILL HAPPEN HERE WHEN FAMILY IS UNASSIGNED??
    
    # if no upstream ranks, nothing more to do at this level
    if (!length(upstream)) next
    
    # 2) build ICTV lookup once per rank
    ictv_lookup <- ictv %>%
      filter(!is.na(.data[[rnk]])) %>%
      select(all_of(c(rnk, upstream))) %>%
      distinct()
    
    # 3) join ICTV info onto df for all rows at this rank
    unassigneders <- df %>%
      filter(!!sym(rnk) == "Unassigned")
    
    df <- df %>%
      filter(!!sym(rnk) != "Unassigned") %>%
      left_join(ictv_lookup, by = rnk, suffix = c("", ".ictv"))
    
    # 4) overwrite upstream ranks only where ICTV has a non-NA value
    for (u in upstream) {
      ictv_u <- paste0(u, ".ictv")
      if (!ictv_u %in% names(df)) next
      
      replace_idx <- !is.na(df[[ictv_u]])
      df[[u]][replace_idx] <- df[[ictv_u]][replace_idx]
      
      df[[ictv_u]] <- NULL
    }
    df <- bind_rows(df, unassigneders)
  }
  
  df
}



# finds cases when the rank of interest has multiple upstream ranks
check_fd <- function(df, lower, higher, drop_na_higher = TRUE,
                     ignore_lower_values = c("Unclassified", "Unassigned"),
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
library(stringr)
library(purrr)
library(tidyverse)
#############################################################
# 2. Load Input Data
#############################################################
taxa <- read.table('06.CLEAN_DATA/Consensus_taxonomy_121062_not_curated.txt', sep='\t', header=T)
unique_taxa <- data.frame("consensus_lineage" = unique(taxa$consensus_lineage))

ictv <- readxl::read_xlsx('06.CLEAN_DATA/ICTV_Master_Species_List_2024_MSL40.v2.xlsx', sheet=2, col_types = "text")
ictv <- ictv[,colnames(ictv) %in% c('Realm', 'Kingdom', 'Phylum',  'Class', 'Order', 'Family',  'Genus',  'Species', 'Genome')]
ictv[is.na(ictv)] <- "Unassigned"
# previous ICTV releases:
msl38 <- readxl::read_xlsx('06.CLEAN_DATA/Previous_ICTV_MSL/ICTV_Master_Species_List_2022_MSL38.v3.xlsx', sheet=4, col_types = "text")
msl39 <- readxl::read_xlsx('06.CLEAN_DATA/Previous_ICTV_MSL/ICTV_Master_Species_List_2023_MSL39.v4.xlsx', sheet=4, col_types = "text")
msl40 <- readxl::read_xlsx('06.CLEAN_DATA/ICTV_Master_Species_List_2024_MSL40.v2.xlsx', sheet=4, col_types = "text")

ICTV_renaming <- rbind(msl38, msl39, msl40[,c('Rank', 'Old Name', 'New Name', 'proposal')])
ICTV_renaming$old_tax_name <- gsub('.*;', '', ICTV_renaming$`Old Name`)
ICTV_renaming$new_tax_name <- gsub('.*;', '', ICTV_renaming$`New Name`)
ICTV_renaming <- ICTV_renaming[!is.na(ICTV_renaming$Rank),]
table(ICTV_renaming$Rank)

rm(list = c("msl38", "msl39", "msl40"))
#############################################################
# 3.1 Analysis: renaming ranks to reflect recent ICTV changes
#############################################################
# splitting to ranks:
unique_taxa <- unique_taxa %>%
  separate(
    col    = consensus_lineage,
    into = c('Life', 'Realm', 'Kingdom', 'Phylum',  'Class', 'Order', 'Family',  'Genus',  'Species', 'Strain'),
    sep    = ";",
    fill   = "right",               
    remove = FALSE                  
  ) %>%
  mutate(across(
    c('Life', 'Realm', 'Kingdom', 'Phylum',  'Class', 'Order', 'Family',  'Genus',  'Species', 'Strain'),
    ~ coalesce(na_if(str_trim(.), ""), "Unclassified")
  ))

# complex cases:
# family>order
unique_taxa$Order[unique_taxa$Family == "Autographiviridae"] <- "Autographivirales"
unique_taxa$Order[unique_taxa$Family == "Lavidaviridae"] <- "Lavidavirales"
# maybe I'll regret about it later:
unique_taxa$Family[unique_taxa$Family == "Autographiviridae"] <- "Unclassified"
unique_taxa$Family[unique_taxa$Family == "Lavidaviridae"] <- "Unclassified"

# merging cases (not documented in MSL?)
unique_taxa$Species[unique_taxa$Species == "Torque teno virus 8"] <- "Torque teno virus 7"
unique_taxa$Species[unique_taxa$Species == "Torque teno virus 16"] <- "Torque teno virus 15"
unique_taxa$Species[unique_taxa$Species == "Megavirus chiliensis"] <- "Megavirus chilense"

# Somewhat auto-renaming:
modified_taxa <- unique_taxa %>% ungroup()
modified_taxa <- rename_all_ranks(modified_taxa, ictv, ICTV_renaming)

#modified_taxa[is.na(modified_taxa)] <- "Unassigned"

#############################################################
# 3.2 Analysis: renaming placeholders introduced by VITAP
#############################################################
modified_taxa <- rename_vitap_placeholders(
  modified_taxa,
  ictv,
  ICTV_renaming,
  ranks
)

#############################################################
# 3.3 Analysis: keeping only ICTV-approved ranks
#############################################################
modified_taxa_clean <- fill_from_ictv_patched(modified_taxa,
                                      ictv,
                                      ranks)

#############################################################
# 3.4 Analysis: check that there is no conflicts
#############################################################

# all combos "given rank - its upstream ranks"
pairs_to_check <- tidyr::expand_grid(
  lower  = ranks,
  higher = ranks
) %>%
  filter(match(lower, ranks) < match(higher, ranks))

fd_results <- pairs_to_check %>%
  pmap(function(lower, higher) check_fd(modified_taxa_clean, lower, higher))

# reporting offender checker results:
report <- map_df(fd_results, function(res) {
  tibble(
    rule = paste0(res$lower, " → ", res$higher),
    offenders = n_distinct(res$offenders[[res$lower]] %||% character()),
    ok = res$ok
  )
})

print(report, n=100) 

# To inspect specific violations, e.g. species -> genus (not the case anymore)
get_res <- function(results, lower, higher) {
  keep(results, ~ .x$lower == lower && .x$higher == higher)[[1]]
}

res_family_order <- get_res(fd_results, "Family", "Order") # caused by Unassigned -> changed the ignored values

rm(list = c("pairs_to_check", "res_family_order"))
#############################################################
# 3.5 Analysis: check if there is a difference in N assigned
#############################################################
taxa_change_number <- as.data.frame(matrix(NA, nrow=10, ncol=3))
colnames(taxa_change_number) <- c('Domain', "perc_classified_before", "perc_classified_after")
taxa_change_number$Domain <- c('Life', rev(ranks), 'Strain')

for (rank in taxa_change_number$Domain) {
  
  taxa_change_number[taxa_change_number$Domain==rank,"perc_classified_before"] <- 
  round(sum(unique_taxa[,rank]!="Unclassified") / 1328 * 100, 2)
  
  taxa_change_number[taxa_change_number$Domain==rank,"perc_classified_after"] <- 
    round(sum(modified_taxa_clean[,rank]!="Unclassified") / 1328 * 100, 2)
  
} # -> there definitely is, but it was expected

#############################################################
# 3.6 Analysis: write new taxonomy for vOTUs
#############################################################
order_taxa <- colnames(taxa)

taxa <- merge(taxa, modified_taxa_clean, by="consensus_lineage")

taxa <- taxa[,c(order_taxa, colnames(modified_taxa_clean)[-1])]
#############################################################
# 4. OUTPUT
#############################################################
write.table(taxa, "06.CLEAN_DATA/Consensus_taxonomy_121062_ICTV_curated.txt", 
            sep='\t', quote=F, row.names=F)

write.table(taxa_change_number, "06.CLEAN_DATA/Table_change_in_Unassigned_after_ICTV_curation.txt",
            sep='\t', quote=F, row.names=F)


