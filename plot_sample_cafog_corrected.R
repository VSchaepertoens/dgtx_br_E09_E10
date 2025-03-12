library(tidyverse)

# load cafog corrected data -----------------------------------------------
# load abundances using a for loop

abundance_data <- NULL
coefs <- c("E09", "E10", "NISTmab")

for (coef in coefs) {
  file_path <- paste0("analysis/cafog/",coef,"/results.csv")
  
  abundance_data <- rbind(abundance_data,
                          read_csv(file_path,
                                   n_max = 13) %>%
                            mutate(experiment = coef)
  )
}

glycoforms_ordered <- c("A2G0F/A2G0F",
                "A2G1F/A2G0F",
                "A2G0/A2G0F",
                "A2G1F/A2G1F",
                "none/A2G0F",
                "A2G0/A2G0",
                "A2G2F/A2G1F", 
                "A2G2F/A2G2F",
                "none/A2G1F",
                "A2G2F/A2G2F +1 Hex",
                "A2G2F/A2G2F +2 Hex",
                "A2G2F/A2G2F +3 Hex",
                "none/A2G2F",
                "A2G0F/A2G0F", #starting second rep
                "A2G1F/A2G0F",
                "A2G0/A2G0F",
                "A2G1F/A2G1F",
                "none/A2G0F",
                "A2G0/A2G0",
                "A2G2F/A2G1F", 
                "A2G2F/A2G2F",
                "none/A2G1F",
                "A2G2F/A2G2F +1 Hex",
                "A2G2F/A2G2F +2 Hex",
                "A2G2F/A2G2F +3 Hex",
                "none/A2G2F",
                "A2G1F/A2G0F", # third rep
                "A2G1F/A2G1F",
                "A2G0F/A2G0F",
                "A2G2F/A2G1F",
                "A2G0F/A1G0F",
                "A2G2F/A2G2F",
                "A2G2F/A2G2F +1 Hex", 
                "A2G2F/A2G2F +2 Hex",
                "A1G1F/A1G0",
                "none/A2G1F",
                "none/A2G0F",
                "A2G2F/A2G2F +3 Hex",
                "none/A2G2F"
)

data_to_plot <- abundance_data %>%
# mutate(modcom_name = rep(glycoforms_ordered, length(coefs))) %>%
  mutate(modcom_name = glycoforms_ordered) %>%
  select(modcom_name, corr_abundance, corr_abundance_error, experiment) %>%
  filter(!str_detect(modcom_name, "Hex")) %>%
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


data_to_plot %>%
  ggplot(aes(modcom_name, corr_abundance)) +
  geom_col(
    aes(y = corr_abundance, fill = experiment),
    position = position_dodge(.9),
  ) +
  geom_errorbar(
    aes(
      ymin = corr_abundance - corr_abundance_error,
      ymax = corr_abundance + corr_abundance_error,
      group = experiment
    ),
    position = position_dodge(.9),
    width = .5,
    linewidth = .25
  ) +
  scale_fill_brewer(palette = "Dark2") +
  xlab("") +
  ylim(0, 75) +
  ylab("fractional abundance (%)") +
  labs(title = "Hexose bias corrected glycoforms - intact") +
  geom_hline(yintercept = 0, linewidth = .35) +
  coord_flip() +
  theme_bw() +
  theme(text = element_text(size = 16, 
                            face = "bold", 
                            family = "sans"),
        axis.text.y = element_text(colour = "black", hjust = 0.5),
        axis.text = element_text(colour = "black"),
        axis.ticks.y = element_blank(),
        plot.title = element_text(hjust = 1),
        legend.title = element_blank(),
        legend.position = "top",
        panel.border = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
  )


ggsave(filename = "figures/fractional_abundance_cpb_corrected_ref.png",    
       height = 160,
       width = 160,
       units = "mm",
       dpi = 600)


# # add data for reference mab and plot  ------------------------------------
# this step not necessary, as in the end I already integrate after the correction, but just in case it is still here
# 
# load("analysis/abundance_data_cpb_notcorrected_ref.RData")
# 
# abundance_data_averaged <- abundance_data_averaged %>%
#   filter(experiment_enzyme %in% "NISTmab_CPB") %>%
#   rename(corr_abundance = frac_abundance, corr_abundance_error = error) %>%
#   mutate(experiment = "NISTmab") %>%
#   select(modcom_name, experiment, corr_abundance, corr_abundance_error)
# 
# data_to_plot <- rbind(data_to_plot, abundance_data_averaged) %>%
#   mutate(modcom_name = factor(modcom_name, levels = c("none/A2G0F",
#                                                       "none/A2G2F",
#                                                       "none/A2G1F",
#                                                       "A1G1F/A1G0",
#                                                       "A2G0/A2G0",
#                                                       "A2G0F/A1G0F",
#                                                       "A2G0/A2G0F",
#                                                       "A2G0F/A2G0F",
#                                                       "A2G1F/A2G0F",
#                                                       "A2G1F/A2G1F",
#                                                       "A2G2F/A2G1F",
#                                                       "A2G2F/A2G2F",
#                                                       "A2G2F/A2G2F +1 Hex",
#                                                       "A2G2F/A2G2F +2 Hex",
#                                                       "A2G2F/A2G2F +3 Hex" 
#                                                       )))
# 
# data_to_plot %>%
#   ggplot(aes(modcom_name, corr_abundance)) +
#   geom_col(
#     aes(y = corr_abundance, fill = experiment),
#     position = position_dodge2(preserve = "single"),
#   ) +
#   geom_errorbar(
#     aes(
#       ymin = corr_abundance - corr_abundance_error,
#       ymax = corr_abundance + corr_abundance_error,
#       group = experiment
#     ),
#     position = position_dodge2(preserve = "single"),
#     # width = .5,
#     linewidth = .25
#   ) +
#   scale_fill_brewer(palette = "Dark2") +
#   xlab("") +
#   ylim(0, 75) +
#   ylab("fractional abundance (%)") +
#   labs(title = "Hexose bias corrected glycoforms - intact") +
#   geom_hline(yintercept = 0, linewidth = .35) +
#   # coord_flip() +
#   theme_bw() +
#   theme(text = element_text(size = 16, 
#                             face = "bold", 
#                             family = "sans"),
#         axis.text.y = element_text(colour = "black", hjust = 0.5),
#         axis.text = element_text(colour = "black"),
#         axis.ticks.x = element_blank(),
#         axis.text.x = element_text(angle = 90, vjust = .5, hjust = 1),
#         plot.title = element_text(hjust = 0.5),
#         legend.title = element_blank(),
#         panel.border = element_blank(),
#         panel.grid.major.x = element_blank(),
#         # panel.grid.minor = element_blank(),
#   )
# 
# 
# ggsave(filename = "figures/frac_ab_barplot_cafog_corrected_samples_ref_reordered_horizontal.png",    
#        height = 150,
#        width = 200,
#        units = "mm",
#        dpi = 600)






  