library(tidyverse, warn.conflicts = FALSE)
library(RColorBrewer)
library("scales")    
library(ComplexHeatmap)
library(circlize)

# load an overview table of data & analysis paths -------------------------

samples_table <- read_csv("data/20240209_Nistmab_firstMeasurements/overview_sample_merged.csv") 
reference_table <- read_csv("data/20240209_Nistmab_firstMeasurements/overview_reference_merged.csv") %>%
  mutate(tech_replicate = c(1,2,1,2,3,4,5))


# load abundances using a for loop  ---------------------------------------

abundance_data <- NULL
for (i in 1:nrow(samples_table)) {
  file_path <- paste0(samples_table[i, "analysis_path"], "/frac_ab_tb_cs50.csv")

abundance_data <- rbind(abundance_data,
                        read_csv(file_path) %>%
                          mutate(experiment_enzyme = samples_table$experiment_enzyme[i],
                                 tech_replicate = samples_table$technical_replicate.x[i]
                                 )
                        )
}

abundance_data_ref <- NULL
for (i in 1:nrow(reference_table)) {
  file_path <- paste0(reference_table[i, "analysis_path"], "/frac_ab_tb_cs50.csv")
  
  abundance_data_ref <- rbind(abundance_data_ref,
                          read_csv(file_path) %>%
                            mutate(experiment_enzyme = "NISTmab_CPB",
                                   tech_replicate = reference_table$tech_replicate[i]
                                   )
  )
}
#to remove two datapoints from Tom's measurements
abundance_data_ref <- abundance_data_ref[-(1:24),]

abundance_data_sample_ref <- rbind(abundance_data, abundance_data_ref)

# plot heatmap ------------------------------------------------------------

## wrangle dataframe into matrix -------------------------------------------
data.matrix <- abundance_data %>%
  mutate(sample_name = paste(CHO_cell_variant_bio_replicate,tech_replicate, sep = "_")) %>%
  select("sample_name", "modcom_name", "frac_ab") %>%
  spread(key = sample_name, value = frac_ab) %>%
  column_to_rownames("modcom_name") %>%
  as.matrix()

data.matrix <- abundance_data_sample_ref %>%
  filter(!grepl("_I", experiment_enzyme)) %>%
  mutate(sample_name = paste(experiment_enzyme,tech_replicate, sep = "_")) %>%
  select(modcom_name, frac_ab, sample_name) %>%
  pivot_wider(names_from = sample_name, values_from = frac_ab) %>%
  column_to_rownames("modcom_name") %>%
  replace(is.na(.), 0) %>%
  arrange(c(15, 14,13, 8, 9, 11, 7, 6, 5, 4, 3, 2, 1, 12,10)) %>%
  as.matrix()

## calculate z-score & plot heatmap -------------------------------------

scaled.data.matrix = t(scale(t(data.matrix))) # for scaling by row  

#check for sanity
mean(data.matrix[1,])
sd(data.matrix[1,])
(data.matrix[1] - mean(data.matrix[1,]))/sd(data.matrix[1,])
(data.matrix[1,2] - mean(data.matrix[1,]))/sd(data.matrix[1,])


BASE_TEXT_SIZE_PT <- 5
ht_opt(
  simple_anno_size = unit(1.5, "mm"),
  COLUMN_ANNO_PADDING = unit(1, "pt"),
  DENDROGRAM_PADDING = unit(1, "pt"),
  HEATMAP_LEGEND_PADDING = unit(1, "mm"),
  ROW_ANNO_PADDING = unit(1, "pt"),
  TITLE_PADDING = unit(2, "mm"),
  heatmap_row_title_gp = gpar(fontsize = BASE_TEXT_SIZE_PT),
  heatmap_row_names_gp = gpar(fontsize = BASE_TEXT_SIZE_PT),
  heatmap_column_title_gp = gpar(fontsize = BASE_TEXT_SIZE_PT),
  heatmap_column_names_gp = gpar(fontsize = BASE_TEXT_SIZE_PT),
  legend_labels_gp = gpar(fontsize = BASE_TEXT_SIZE_PT),
  legend_title_gp = gpar(fontsize = BASE_TEXT_SIZE_PT),
  legend_border = FALSE
)

#set the correct color scheme
min(scaled.data.matrix)
max(scaled.data.matrix)
f1 = colorRamp2(seq(-max(abs(scaled.data.matrix)),
                    max(abs(scaled.data.matrix)),
                    length = 9),
                c("seagreen4",
                  "seagreen3",
                  "seagreen2",
                  "seagreen1",
                  "gold",
                  "darkorchid1",
                  "darkorchid2",
                  "darkorchid3",
                  "darkorchid4"),
                space = "RGB")
#set the correct color scheme
png(filename = "figures/heatmap_cpb_notcorrected_ref.png",    
    height = 6.5,
    width = 8,
    units = "cm",
    res = 600)


Heatmap(scaled.data.matrix,
        col = f1,
        cluster_rows = FALSE,
        rect_gp = gpar(col = "white", lwd = 2),
        name = "z-score",
        row_gap = unit(2, "pt"),
        column_gap = unit(2, "pt"),
        width = unit(2, "mm") * ncol(scaled.data.matrix) + 5 * unit(2, "pt"), # to make each cell a square
        height = unit(2, "mm") * nrow(scaled.data.matrix) + 5 * unit(2, "pt"), # to make each cell a square
        show_row_names = TRUE
        )

dev.off()


# calculate mean and sd and plot ------------------------------------------
abundance_data_averaged <- abundance_data_sample_ref %>% 
  group_by(modcom_name, experiment_enzyme) %>%
  summarise(frac_abundance = mean(frac_ab),
            error = sd(frac_ab)) 

save(abundance_data_ref, 
     abundance_data_averaged, 
     file = "analysis/abundance_data_cpb_notcorrected_ref.RData")

# load("analysis/Jan_2024/abundance_data.RData")
# load("analysis/Jan_2024/abundance_data_selected_glycans_intact.RData")

# plot bar chart ----------------------------------------------------------

abundance_data_averaged <- abundance_data_averaged %>%
  filter(!grepl("_I", experiment_enzyme)) %>%
  separate(experiment_enzyme, 
           into = c("experiment", "enzyme"), 
           sep = "_", 
           remove = FALSE) %>%
  as.data.frame() %>%
  mutate(modcom_name = factor(modcom_name, levels = c("A2G2F/A2G2F +3 Hex",
                                                      "A2G2F/A2G2F +2 Hex",
                                                      "A2G2F/A2G2F +1 Hex",  
                                                      "A2G2F/A2G2F",
                                                      "A2G2F/A2G1F",
                                                      "A2G1F/A2G1F",
                                                      "A2G1F/A2G0F",
                                                      "A2G0F/A2G0F",
                                                      "A2G0/A2G0F",
                                                      "A2G0F/A1G0F",
                                                      "A2G0/A2G0",
                                                      "A1G1F/A1G0",
                                                      "none/A2G2F", 
                                                      "none/A2G1F", 
                                                      "none/A2G0F")))


abundance_data_averaged %>%
  ggplot(aes(modcom_name, frac_abundance)) +
  geom_col(
    aes(y = frac_abundance, fill = experiment),
    # position = position_dodge(.9),
    position = position_dodge2(preserve = "single"),
    # width = 0.5
  ) +
  geom_errorbar(
    aes(
      ymin = frac_abundance - error,
      ymax = frac_abundance + error,
      group = experiment
    ),
    # position = position_dodge(.9),
    position = position_dodge2(preserve = "single"),
    # width = .5,
    linewidth = .25
  ) +
  scale_fill_brewer(palette = "Dark2") +
  scale_y_continuous(breaks = c(0,10,20,30,35),
                     labels = number_format(accuracy = 1)) +
  xlab("") +
  ylim(0, 75) +
  ylab("fractional abundance (%)") +
  labs(title = "") +
  geom_hline(yintercept = 0, linewidth = .35) +
  coord_flip() +
  theme_bw() +
  theme(text = element_text(size = 16, 
                            face = "bold", 
                            family = "sans"),
        axis.text = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black", hjust = .5),
        axis.ticks.y = element_blank(),
        legend.title = element_blank(),
        panel.border = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
)


ggsave(filename = "figures/fractional_abundance_cpb_notcorrected_ref.png",    
       height = 160,
       width = 160,
       units = "mm",
       dpi = 600)

