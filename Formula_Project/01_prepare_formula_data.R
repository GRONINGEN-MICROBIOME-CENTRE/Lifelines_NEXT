# Load libraries
library(dplyr)
library(tidyverse)
library(stringr)

# Define answer lists
Negative_answers <- c('Geen', 'geen', 'N.v.t.','Geen flesvoeding, enkel BV','Geen flesvoeding','geen flesvoeding','Nvt','niet van toepassing','Borstvoeding ','Geen flesvoeding','Bv','bv','niet ','borstvoeding','Nog niks gekocht ','niks krijgt borstvoedign','Nvt alleen de eerste dagen in ziekenhuis bijgevoed ','Geen, maar moet dit blijkbaar invullen','Geen, want borstvoeding','niet van toepassing, alleen borstvoeding','geen. alleen borstvoeding','Niet ','Alleen borstvoeding','Niet van toepassing','Nvt, uitsluitend borstvoeding','Niks ')

answer_lists <- list(
  AH_standaard_1 = c("Ah zuigelingenmelk","Voeding van AH", "Eigen merk Albert heijn", "Albert Heijn", "AH eigen merk", "Ah zuigelingen flesvoeding","Ah zuigelingenmelk", "Albert Heijn", "Albert Heijn eigen merk standaard 1 met melkvet", "Eigen merk AH", "AH huismerk", "Albert Heijn huismerk", "AH standaard 1", "Albert Heijn Huismerk", "Albert Heijn standaard 1", "Albert Heijn zuigelingenmelk 1", "AH zuigelingenmelk 1", "AH huismerk", "Albertheijn standaard 1"),
  Kruidvat_1 = c("Kruidvat","Standaard 1 Volledige Zuigelingenvoeding","Kruidvat standaard 1","Kruidvat standaard","Kruidvat huismerk 1", "Kruidvat 1", "kruidvat", "Kruidvat standaard 1 complete zuigelingenvoedinf", "kruitvat huismerk", "Kruitvat standaard 1","Kruidvat eigen merk 1","Kruidvat eigen merk","Kruitvat nummer 1","Kruidvat huismerk nr 1","Kruidvat 0-6 maand","Kruidvat zuigelingenvoeding 1","volledige zuigelingenvoeding van kruidvat eigen merk", "Kruidvat standaard 1 volledige zuigelingenvoeding","kruidvat 1","kruidvat standaard 1","Kruidvat volledige zuigelingenvoeding 1", "kruidvat eigen merk","Kruidvat eigen merk", "Kruidvat Standaard 1","Kruidvat Standaard 1 Volledige Zuigelingenvoeding"),
  Nutrilon_AR_1 = "Nutrilon_AR",
  Jumbo_standaard_1 = c("Jumbo standaard 1", "Jumbo eigen merk", "Jumbo", "Jumbo zuigelingenmelk", "huismerk jumbo", "Jumbo zuigeling", "Jumbo huismerk", "jumbo eigen merk", "Jumbo zuigelingsmelk 0-6 maand"),
  HIPP_PRE_HA = c("Hipp pre HA", "Hipp Pre", "Hipp Pre HA", "Hipp bio PRE", "Hipp Pre HA"),
  Nutramigen_1_LG = c("Nutramigen lgg", "Nutramigen 1 lgg", "Nutramigen LGG Hypoallergeen voeding", "Nutramigen 1Lgg"),
  AH_bio_1 = c("AH Biologisch", "Ah huismerk bio. 1 week x als bijvoeding gegeven", "Albert Heijn biologisch", "AH bio", "Ah biologisch 1"),
  Hero_HA = c("Hero Hypo Allergeen","Hero hypo allergeen")
)

# Function to process each dataset
process_formula_data <- function(df, timepoint, formula_col, answer_lists, negative_answers) {
  df <- df %>%
    rename(other_formula = !!sym(formula_col)) %>%
    mutate(other_formula = ifelse(other_formula %in% negative_answers, 0, other_formula))
  
  for (formula in names(answer_lists)) {
    col_name <- paste0(formula, "_", timepoint)
    df[[col_name]] <- 0
    df[[col_name]] <- ifelse(df$other_formula %in% answer_lists[[formula]], 1, df[[col_name]])
    df$other_formula <- ifelse(df[[col_name]] == 1, "", df$other_formula)
  }
  
  df %>%
    mutate(timepoint = timepoint,
           other_formula = ifelse(other_formula != "" & other_formula != "0", 1, other_formula),
           other_formula = na_if(other_formula, ""),
           other_formula = as.numeric(other_formula)) %>%
    pivot_longer(cols = -c("NG_ID", "timepoint"), names_to = "formula", values_drop_na = TRUE) %>%
    mutate(formula = str_replace(formula, "_[^_]*$", ""))
}

# Load data
W2_raw <- read.table("data/20250211_binary_table_W2.txt", header = TRUE, sep = "\t")
M1_raw <- read.table("data/20250211_binary_table_M1.txt", header = TRUE, sep = "\t")
M2_raw <- read.table("data/20250211_binary_table_M2.txt", header = TRUE, sep = "\t")
M3_raw <- read.table("data/20250211_binary_table_M3.txt", header = TRUE, sep = "\t")

# Process each timepoint
W2 <- process_formula_data(W2_raw, "W2", "CHFFQ5TXT_W2", answer_lists, Negative_answers)
M1 <- process_formula_data(M1_raw, "M1", "CHFFQ5TXT_M1", answer_lists, Negative_answers)
M2 <- process_formula_data(M2_raw, "M2", "CHFFQ5TXT_M2", answer_lists, Negative_answers)
M3 <- process_formula_data(M3_raw, "M3", "CHFFQ5TXT_M3", answer_lists, Negative_answers)

# Combine and save
data_combined <- bind_rows(W2, M1, M2, M3)
saveRDS(data_combined, file = "output/df_data_combined.rds")
