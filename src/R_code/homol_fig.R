library(tidyverse)
library(cowplot)
library(ggforce)
#getting the average % homology for all proteins in the alignment.

homol_data <- read.csv(file = "./percent_homol.csv", header=TRUE, sep=",")

plot_a <- homol_data %>%
  ggplot(aes(x = percent_homol)) +
  geom_density(fill = "seashell3", color = "seashell4") +
  scale_x_continuous(
    name = "% Homology per Protein Alignment",
    limits = c(0.0, 1.01),
    breaks = seq(from = 0.0, to = 1.0, by = 0.1),
    expand = c(0, 0)) +
  scale_y_continuous(
    name = "Density",
    limits = c(0, 5.5),
    breaks = seq(from = 0, to = 5, by = 1),
    expand = c(0, 0)) + 
  theme_cowplot()
plot_a

ggsave(filename = "../../analysis/figures/figure_hom.png", plot = plot_a, width = 8, height = 5)


# now getting the average class percents per column in every alignment. 

class_fractions <- read.csv(file = "./stats_align_all.csv", header=TRUE, sep=",")

class_fractions <- class_fractions %>%
  select(position, gene, c(q_aliphatic:q_proline)) %>%
  pivot_longer(cols = c(q_aliphatic:q_proline), names_to = "class", values_to = "fraction") %>%
  group_by(gene, class) %>%
  summarise(mean_frac = mean(fraction))

# order the level by the mean of all proteins.
class_fractions <- class_fractions %>%
  group_by(class) %>%
  summarise(mean_all_prot = mean()) %>%
  fct_reorder(class, mean_all_prot)


plot_b <- class_fractions %>%
  ggplot(aes(y = mean_frac, x = fct_rev(fct_reorder(class, mean_frac)))) +
  geom_violin(color = "seashell4", fill = "seashell2") +
  geom_sina(size = 0.7, seed = 123) + 
  scale_x_discrete(
    name = "Amino Acid Class",
    labels = c("Aliphatic", "Polar", "Positive", "Negative", "Aromatic", "Proline")) +
  scale_y_continuous(
    name = "Fraction per Alignmnent Position \n (Averaged by Protein)",
    limits = c(0, 0.55),
    breaks = seq(from = 0.0, to = 0.5, by = 0.1),
    expand = c(0, 0)) +
  theme_bw(12) +
  theme(
    axis.text = element_text(color = "black", size = 12),
    panel.grid.minor = element_blank())

plot_b

ggsave(filename = "../../analysis/figures/class_frac.png", plot = plot_b, width = 8, height = 5)







