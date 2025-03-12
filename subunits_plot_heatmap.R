library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(viridis)


# load_data ---------------------------------------------------------------

data <- read_csv('data/RelQuantIntact01.csv') %>%
  column_to_rownames(var = "...1")


# heatmap -----------------------------------------------------------------
## plot heatmap of raw data ------------------------------------------------

data.matrix <- as.matrix(data)

Heatmap(data.matrix)
Heatmap(data.matrix, col = plasma(100))
Heatmap(data.matrix, col = rev(rainbow(10)))

## calculate z-score & plot heatmap -------------------------------------------------------

scaled.data.matrix = t(scale(t(data.matrix))) # for scaling by row  

#check for sanity
mean(data.matrix[1,])
sd(data.matrix[1,])
(data.matrix[1] - mean(data.matrix[1,]))/sd(data.matrix[1,])
(data.matrix[1,2] - mean(data.matrix[1,]))/sd(data.matrix[1,])

min(scaled.data.matrix)
max(scaled.data.matrix)
col_fun = colorRamp2(c(-1.3, 0, 1.6), c("green", "white", "red"))
Heatmap(scaled.data.matrix, col = col_fun)
Heatmap(scaled.data.matrix, col = plasma(100))
Heatmap(scaled.data.matrix, col = rev(rainbow(10)))


# build correct table -----------------------------------------------------

data <- data %>%
  rownames_to_column("experiment") %>%
  separate_wider_delim("experiment",
                       delim = "_",
                       names = c("day","month", "year", "experiment")
                       ) %>%
  pivot_longer(cols = colnames(data),
               names_to = c("subunit","replicate"),
               names_sep = "_",
               values_to = "peak_area")

data$replicate <- str_replace(data$replicate, "PeakArea", "")

#calculate group means
data_averaged <- data %>%
  group_by(experiment,subunit) %>%
  summarise(mean_peak_area = mean(peak_area)) %>%
  group_by(experiment) %>%
  mutate(percent = (mean_peak_area/sum(mean_peak_area))) %>%
  filter(experiment %in% c("E09", "E10"))

#convert 'subunit' to factor and specify level order ## important for plotting order
data_averaged$subunit <- factor(data_averaged$subunit, levels = c('Intact', 'LC/LC', 'LC'))

# plot stacked bar chart -------------------------------------------------

ggplot(data_averaged, aes(y = experiment, x = mean_peak_area, fill = subunit)) + 
  geom_bar(stat = "identity", position = "fill") +
  xlab("peak_area (%)") +
  scale_fill_brewer(palette = "Accent") +
  scale_y_discrete(limits = rev) +
  theme_bw() +
  theme(text = element_text(size = 12, 
                            # face = "bold", 
                            family = "sans"),
        axis.text = element_text(colour = "black"),
        legend.position = "top",
        axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank()
        ) 

ggsave("figures/subunits_peak_area_percent.png", 
       dpi = 600,
       width = 100,
       height = 50,
       units = "mm",
       bg = "white")
























