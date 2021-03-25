library(tidyverse)
library(cowplot)
#getting the average % homology for all proteins in the alignment.

homol_data <- read.csv(file = "./percent_homol.csv", header=TRUE, sep=",")

plot_a <- homol_data %>%
  ggplot(aes(x = percent_homol)) +
  geom_density(fill = "thistle4", alpha = 0.6) +
  scale_x_continuous(
    name = "% Homology per Protein Alignment",
    limits = c(0.0, 1.0),
    breaks = seq(from = 0.0, to = 1.0, by = 0.1),
    expand = c(0, 0)) +
  scale_y_continuous(
    name = "Density",
    limits = c(0, 5.5),
    breaks = seq(from = 0, to = 5, by = 1),
    expand = c(0, 0)) + 
  theme_cowplot()

ggsave(filename = "../../analysis/figures/figure_hom.png", plot = plot_a, width = 8, height = 5)
