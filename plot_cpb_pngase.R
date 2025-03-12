library(tidyverse, warn.conflicts = FALSE)
library(RColorBrewer)
library("scales")    
library(ComplexHeatmap)
library(circlize)

# load an overview table of data & analysis paths -------------------------

samples_table <- read_csv("data/20240209_Nistmab_firstMeasurements/overview_cpb_pngase_merged.csv") 

# load abundances using a for loop  ---------------------------------------

abundance_data <- NULL
for (i in 1:nrow(samples_table)) {
  file_path <- paste0(samples_table[i, "analysis_path"], "/frac_ab_tb_cs50.csv")

abundance_data <- rbind(abundance_data,
                        read_delim(file_path) %>%
                          mutate(experiment_enzyme = samples_table$experiment_enzyme[i],
                                 tech_replicate = samples_table$technical_replicate.x[i]
                                 )
                        )
}


# calculate mean and sd and plot ------------------------------------------
abundance_data_averaged <- abundance_data %>% 
  group_by(modcom_name, experiment_enzyme) %>%
  summarise(frac_abundance = mean(frac_ab),
            error = sd(frac_ab)) %>%
  mutate(experiment_enzyme = case_when(
    experiment_enzyme == "mAbPAC-RP-2d1mm_deglyco_intact" ~ "NISTmab_CPB_PNGF",
    # experiment_enzyme == "old_value2" ~ "new_value2",
    TRUE ~ experiment_enzyme  # Keep other values unchanged
  ))


save(abundance_data, 
     abundance_data_averaged, 
     file = "analysis/abundance_data_cpb_pngase.RData")

# plot bar chart ----------------------------------------------------------

abundance_data_averaged <- abundance_data_averaged %>% 
  separate_wider_delim(experiment_enzyme,
                       delim = "_",
                       names = c("experiment","enzyme", "enzyme2")) %>%
  mutate(modcom_name = factor(modcom_name, levels = c("3xHex", "2xHex","1xHex","none")))

abundance_data_averaged %>%
  ggplot(aes(modcom_name, frac_abundance)) +
  geom_col(
    aes(y = frac_abundance, fill = experiment),
    position = position_dodge(.9),
  ) +
  geom_errorbar(
    aes(
      ymin = frac_abundance - error,
      ymax = frac_abundance + error,
      group = experiment
    ),
    position = position_dodge(.9),
    width = .5,
    linewidth = .25
  ) +
  scale_fill_brewer(palette = "Dark2") +
  scale_y_continuous(breaks = c(0,10,20,30,35),
                     labels = number_format(accuracy = 1)) +
  xlab("") +
  ylim(0, 100) +
  ylab("fractional abundance (%)") +
  labs(title = "") +
  geom_hline(yintercept = 0, linewidth = .35) +
  coord_flip() +
  theme_bw() +
  theme(text = element_text(size = 16, 
                            face = "bold", 
                            family = "sans"),
        axis.text = element_text(colour = "black"),
        axis.ticks.y = element_blank(),
        legend.title = element_blank(),
        panel.border = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
)


ggsave(filename = "figures/fractional_abundance_cpb_pngase.png",    
       height = 160,
       width = 160,
       units = "mm",
       dpi = 600)

