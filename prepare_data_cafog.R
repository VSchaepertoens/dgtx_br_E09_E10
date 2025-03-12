#prepare data for cafog analysis
library(tidyverse)
library(dplyr)
library(stringr)
library(fs)

# load abundance data -----------------------------------------------------

# load("analysis/abundance_data_cpb_notcorrected.RData")
load("analysis/abundance_data_cpb_notcorrected_ref.RData")


glycosylation <-  abundance_data_averaged  %>%
  separate_wider_delim(experiment_enzyme,
                       delim = "_",
                       names = c("experiment","enzyme"),
                       cols_remove = FALSE
                       ) %>%
  filter(enzyme == "CPB")

# load glycation data -----------------------------------------------------


load("analysis/abundance_data_cpb_pngase.RData")

glycation <-  abundance_data_averaged %>%
  separate_wider_delim(experiment_enzyme,
                       delim = "_",
                       names = c("experiment","enzyme1", "enzyme2"),
                       cols_remove = FALSE
  ) 


# Define composition mapping ------------------------------------------------------------------

composition_mapping <- list(
  'A2G2F +1 Hex' = "6 Hex, 4 HexNAc, 1 Fuc",
  'A2G2F +2 Hex' = "7 Hex, 4 HexNAc, 1 Fuc",
  'A2G2F +3 Hex' = "8 Hex, 4 HexNAc, 1 Fuc",
  none = "0 Hex"
)


# for each coef make files for cafog analysis --------------------------------------
coefs <-  unique(glycosylation$experiment)

fs::dir_create(paste0("analysis/cafog/",coefs))

#for loop to create cafog files for each CHO_cell_variant_bio_repliacte
for (coef in coefs) {
  print(coef)
  glycosylation %>%
    filter(experiment == coef) %>%
    select(modcom_name, frac_abundance, error) %>%
    rename(`#glycoform` = modcom_name,
           abundance = frac_abundance) %>%
    write_csv(paste0("analysis/cafog/",coef,"/glycosylation.csv"),
              col_names = TRUE)

  glycation %>%
    filter(experiment == coef) %>%
    select(modcom_name, frac_abundance, error) %>%
    rename(`#count` = modcom_name,
           abundance = frac_abundance) %>%
    mutate(`#count` = as.character(`#count`)) %>%
    mutate(`#count` = str_replace_all(`#count`, c("3xHex" = "3","2xHex" = "2","1xHex" = "1","none" = "0"))) %>%
    write_csv(paste0("analysis/cafog/",coef,"/glycation.csv"),
              col_names = TRUE)
  glycosylation %>%
    filter(experiment == coef) %>%
    select(modcom_name) %>%
    separate(modcom_name, into = c("glycoform_1", "glycoform_2"), sep = "/") %>%
    pivot_longer(cols = c("glycoform_1", "glycoform_2"), names_to = "names", values_to = "glycoforms") %>%
    select(glycoforms) %>%
    unique()  %>%
    mutate(
      composition = case_when(
        # glycoforms == "G1F" ~ composition_mapping[["G1F"]],
        glycoforms == "none" ~ composition_mapping[["none"]],
        glycoforms == "A2G2F +1 Hex" ~ composition_mapping[["A2G2F +1 Hex"]],
        glycoforms == "A2G2F +2 Hex" ~ composition_mapping[["A2G2F +2 Hex"]],
        glycoforms == "A2G2F +3 Hex" ~ composition_mapping[["A2G2F +3 Hex"]],
        TRUE ~ NA_character_  # Handle unmatched cases
      )) %>%
    write_csv(paste0("analysis/cafog/",coef,"/glycan_library.csv"),
                        col_names = TRUE)
} 

## Remove NAs from glycan library manually
## Run CAFOG analysis using subprocess_cafog.ipynb from Anaconda --> vs studio
## Continue with plotting the corrected results --> plot_cafog_corrected.R
