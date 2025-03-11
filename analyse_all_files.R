## ---------------------------
##
## Script name: Quantification of fractional abundances of glycans found in 
## BOKU Herceptin project 2023
##
## Purpose of script: Using fragquaxi package to quantify glycans in samples
##
## Author: Dr. Veronika Schäpertöns
##
## Date Created: 01.02.2024 
##
## Copyright (c) Veronika Schäpertöns, 2023
## Email: veronika.schaepertoens@plus.ac.at
##
## ---------------------------
##
## Notes:
##   
##
## ---------------------------

library(fs)
library(tidyverse, warn.conflicts = FALSE)
library(janitor, warn.conflicts = FALSE)
library(fragquaxi)

# define to either run analysis of sample mab or pngase digested
product <- "sample" #OR product <- "pngase"  #OR product <- "sample" 

directory <- "data/20240209_Nistmab_firstMeasurements/"

# define mab_sequence  --------------------------------

mab_sequence <- "mab_sequence/Nistmab_RM8671_noLys.fasta"

proteins <- define_proteins(
  cNISTmab = mab_sequence,
  .disulfides = 16
)

# depending on product, define regexp pattern and modcoms -----------------

if (product == "sample") {
  #specify pattern to match regexp for filename
  #pattern <-  "ambr.*\\.mzML"
  pattern <-  "E.*\\.mzML"
  
  #obtained from Kathi's script e_Fragquaxi_intact.R
  modcoms <- tribble(
    ~modcom_name    , ~Hex, ~HexNAc, ~Neu5Gc, ~Fuc,  ~PYRRO,
    "none/A2G0F", 3, 4, 0, 1, 2, # "3 Hex, 4 HexNAc, 0 Neu5Gc, 1 Fuc, 2 Pyrro"
    "none/A2G1F", 4, 4, 0, 1, 2, # "4 Hex, 4 HexNAc, 0 Neu5Gc, 1 Fuc, 2 Pyrro"
    "none/A2G2F", 5, 4, 0, 1, 2, # "4 Hex, 4 HexNAc, 0 Neu5Gc, 1 Fuc, 2 Pyrro"
    #"A1G1F/A1G0", 7, 6, 0, 1, 2, # "7 Hex, 6 HexNAc, 0 Neu5Gc, 1 Fuc, 2 Pyrro"
    # "FA1G0/A1G0", 6, 6, 0, 1, 2, # "6 Hex, 6 HexNAc, 0 Neu5Gc, 1 Fuc, 2 Pyrro"
    #"A2G0F/A1G0F", 6, 7, 0, 2, 2,  # "6 Hex, 7 HexNAc, 0 Neu5Gc, 2 Fuc, 2 Pyrro"
    "A2G0/A2G0F", 6, 8, 0, 1, 2, 
    "A2G0/A2G0", 6, 8, 0, 0, 2, 
    "A2G0F/A2G0F", 6, 8, 0, 2, 2, 
    "A2G1F/A2G0F", 7, 8, 0, 2, 2, 
    "A2G1F/A2G1F", 8, 8, 0, 2, 2, 
    "A2G2F/A2G1F", 9, 8, 0, 2, 2, 
    "A2G2F/A2G2F", 10, 8, 0, 2, 2, 
    "A2G2F/A2G2F +1 Hex", 11, 8, 0, 2, 2, 
    "A2G2F/A2G2F +2 Hex", 12, 8, 0, 2, 2, 
    "A2G2F/A2G2F +3 Hex", 13, 8, 0, 2, 2, # "13 Hex, 8 HexNAc, 0 Neu5Gc, 2 Fuc, 2 Pyrro"
    # "FA2G2/FA2G2 +4 Hex - PYRRO", 14, 8, 0, 2, 0, # "14 Hex, 8 HexNAc, 0 Neu5Gc, 2 Fuc, 0 Pyrro"
  ) %>% 
    define_ptm_compositions()
   
 } else if (product == "pngase")  {
   #specify pattern to match regexp for filename 
   pattern <- "PNGaseF.*\\.mzML"
   
   #specify modifications composition
   modcoms <- tribble(
     ~modcom_name, ~Hex, ~HexNAc, ~Fuc, ~Neu5Ac,
     "none", 0, 0, 0, 0,
     "1xHex", 1, 0, 0, 0,
     "2xHex", 2, 0, 0, 0,
     "3xHex", 3, 0, 0, 0, 
   ) %>%
     define_ptm_compositions()
   
 } else if (product == "reference")  {
   #specify pattern to match regexp for filename 
   pattern <- "nist.*\\.mzML"
   
   #specify modifications composition
   #obtained from Kathi's script e_Fragquaxi_intact.R
   modcoms <- tribble(
     ~modcom_name    , ~Hex, ~HexNAc, ~Neu5Gc, ~Fuc,  ~PYRRO,
     "none/A2G0F", 3, 4, 0, 1, 2, # "3 Hex, 4 HexNAc, 0 Neu5Gc, 1 Fuc, 2 Pyrro"
     "none/A2G1F", 4, 4, 0, 1, 2, # "4 Hex, 4 HexNAc, 0 Neu5Gc, 1 Fuc, 2 Pyrro"
     "A1G1F/A1G0", 7, 6, 0, 1, 2, # "7 Hex, 6 HexNAc, 0 Neu5Gc, 1 Fuc, 2 Pyrro"
     # "FA1G0/A1G0", 6, 6, 0, 1, 2, # "6 Hex, 6 HexNAc, 0 Neu5Gc, 1 Fuc, 2 Pyrro"
     "A2G0F/A1G0F", 6, 7, 0, 2, 2,  # "6 Hex, 7 HexNAc, 0 Neu5Gc, 2 Fuc, 2 Pyrro"
     "A2G0F/A2G0F", 6, 8, 0, 2, 2, 
     "A2G1F/A2G0F", 7, 8, 0, 2, 2, 
     "A2G1F/A2G1F", 8, 8, 0, 2, 2, 
     "A2G2F/A2G1F", 9, 8, 0, 2, 2, 
     "A2G2F/A2G2F", 10, 8, 0, 2, 2, 
     "A2G2F/A2G2F +1 Hex", 11, 8, 0, 2, 2, 
     "A2G2F/A2G2F +2 Hex", 12, 8, 0, 2, 2, 
     "A2G2F/A2G2F +3 Hex", 13, 8, 0, 2, 2, # "13 Hex, 8 HexNAc, 0 Neu5Gc, 2 Fuc, 2 Pyrro"
     # "FA2G2/FA2G2 +4 Hex - PYRRO", 14, 8, 0, 2, 0, # "14 Hex, 8 HexNAc, 0 Neu5Gc, 2 Fuc, 0 Pyrro"
   ) %>%
     define_ptm_compositions()
}

# specify paths ------------------------------------------------------------

df <- tibble(mzml_full_path = dir_ls(path = directory,regexp = pattern),) %>%
  separate(mzml_full_path,
           into = c("data", "exp_month","filename"),
           sep = "/",
           remove = FALSE) %>%
  mutate(analysis_path = fs::path("analysis",
                                  exp_month,
                                  gsub("\\..*$", "", filename))) 
  
fs::dir_create(df$analysis_path)

# wrangle out info ---------------------------------------------------

df <- separate_wider_delim(df,
                           filename,
                           delim = "_",
                           names = c("day","month", "year", "experiment", "technical_replicate", "enzyme","aquisition_number"),
                           too_few = "debug") %>%
  unite(experiment_enzyme,
        c("experiment","enzyme"),
        remove = FALSE)
  # unite(CHO_cell_variant_bio_replicate, 
  #       c("CHO_cell_variant","bio_replicate"), 
  #       remove = FALSE) %>%
  # mutate(tech_replicate = rep(c(1,2,3), 14)) #perhaps is redundant

# write_csv(x = df,
#           file = "data/20240209_Nistmab_firstMeasurements/overview_sample.csv")

# import data with information on charge states and rt limits  --------

cs_rt_data <- read_csv('data/20240209_Nistmab_firstMeasurements/rt_seconds_SA.csv') %>%
  separate(sample_name,
           into = c("day", "month", "year", "experiment","technical_replicate","enzyme", "aquisition_number"),
           sep = "_",
           remove = FALSE) %>%
  # filter(subunit == "intact") %>%
  unite(experiment_enzyme,
        c("experiment","enzyme"),
        remove = FALSE)


# merge data sample with data cs and rt -----------------------------------
data_merged <- df %>% 
  left_join(cs_rt_data, by = "experiment_enzyme") 
  #cs_rt_data %>%  select(CHO_cell_variant_bio_replicate, )
  # filter(CHO_cell_variant.x != "A25") # nearly no signal

write_csv(x = data_merged,
          file = "data/20240209_Nistmab_firstMeasurements/overview_sample_merged.csv")

# fragquaxi analysis ------------------------------------------------------


## custom function ---------------------------------------------------------
calculate_abundance <- function(mzml_full_path,
                                rt_start_sec,
                                rt_end_sec, 
                                analysis_path,
                                scan_start,
                                scan_end,
                                ...){
  
  ms_data <- mzR::openMSfile(mzml_full_path)
  print(ms_data)
  print(c(rt_start_sec,rt_end_sec))
  
  pfm_ions <-
    assemble_proteoforms(proteins, modcoms) %>%
    ionize(charge_states = c(42:53), ppm = 300)
  print(dim(pfm_ions)) #check that for every file the correct # of charge states was used
  
  extracted_filename <- str_extract(ms_data@fileName, "(?<=data\\/).*(?=\\.mzML)")
  # single charge states
  plot_ions(
    ms_data,
    ions = pfm_ions,
    scans = scan_start:scan_end,
    xlim = c(3320, 3400)
  )
  
  ggsave(filename = paste0("figures/",extracted_filename,"one_cs.png"),    
         height = 200,
         width = 300,
         units = "mm",
         dpi = 600)
  
  # three charge states
  plot_ions(
    ms_data,
    ions = pfm_ions,
    scans = scan_start:scan_end,
    xlim = c(3110, 3320)
  )
  
  ggsave(filename = paste0("figures/",extracted_filename,"three_cs.png"),    
         height = 200,
         width = 300,
         units = "mm",
         dpi = 600)
  
  abundances <- quantify_ions(ms_data,
                              ions = pfm_ions,
                              rt_limits = c(rt_start_sec,rt_end_sec)
                              ) %>%
    as_tibble() %>% 
    mutate(modcom_name = factor(modcom_name) %>% 
    fct_inorder()
    ) %>% 
    unnest(abundance_data) %>% 
    group_by(modcom_name) %>%
    summarise(abundance = sum(abundance)) %>% 
    mutate(frac_ab = abundance / sum(abundance) * 100,
           file_name = mzml_full_path)

    write_csv(x = abundances,
              file = paste(analysis_path,"frac_ab_tb_cs50.csv",sep = "/")
              )
  
  print('Analysis finished')
}


## apply custom function to dfr --------------------------------------------

pwalk(data_merged[,], calculate_abundance, .progress = TRUE)






          
