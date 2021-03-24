# continuation of figure 2. 
# comparing neff predicted vs. neff natural across different alignments similarities. 
library(tidyverse)
library(cowplot)

#set working directory to:
#C:\Users\avch\Desktop\Natural_var_project\output\output_PSICOV\

# reading csv files
natural_var <- read.csv(file="./stats_align_all.csv", header=TRUE, sep=",")
natural_var_20 <- read.csv(file = "./stats_align_files/stats_align_20.csv", header=TRUE, sep=",")
natural_var_40 <- read.csv(file = "./stats_align_files/stats_align_40.csv", header=TRUE, sep=",")
natural_var_60 <- read.csv(file = "./stats_align_files/stats_align_60.csv", header=TRUE, sep=",")
natural_var_80 <- read.csv(file = "./stats_align_files/stats_align_80.csv", header=TRUE, sep=",")
natural_var_100 <- read.csv(file = "./stats_align_files/stats_align_100.csv", header=TRUE, sep=",")

cnn_var <- read.csv(file="./stats_cnn.csv", header=TRUE, sep=",")
cnn_var2 <- cnn_var %>%
  select(position, gene, n_eff, n_eff_class) %>%
  mutate(group = "predicted")

#0-20
natural_var_20 <- natural_var_20 %>%
  select(position, gene, n_eff, n_eff_class) %>%
  mutate(group = "natural")

  
joined_20 <- rbind(natural_var_20, cnn_var2) %>%
  mutate(perc_sim = "(0-20%]") 

#20-40
natural_var_40 <- natural_var_40 %>%
  select(position, gene, n_eff, n_eff_class) %>%
  mutate(group = "natural")

joined_40 <- rbind(natural_var_40, cnn_var2) %>%
  mutate(perc_sim = "(20-40%]") 

#40-60
natural_var_60 <- natural_var_60 %>%
  select(position, gene, n_eff, n_eff_class) %>%
  mutate(group = "natural") 

joined_60 <- rbind(natural_var_60, cnn_var2) %>%
  mutate(perc_sim = "(40-60%]")  

#60-80
natural_var_80 <- natural_var_80 %>%
  select(position, gene, n_eff, n_eff_class) %>%
  mutate(group = "natural") 

joined_80 <- rbind(natural_var_80, cnn_var2) %>%
  mutate(perc_sim = "(60-80%]") 

#80-100
natural_var_100 <- natural_var_100 %>%
  select(position, gene, n_eff, n_eff_class) %>%
  mutate(group = "natural") 

joined_100 <- rbind(natural_var_100, cnn_var2) %>%
  mutate(perc_sim = "(80-100%]") 

  
all_joined <- rbind(joined_20, joined_40, joined_60, joined_80, joined_100)

all_joined_wide <- all_joined %>%
  select(-n_eff_class) %>%
  pivot_wider(names_from = group, values_from = n_eff)

cor <- all_joined_wide %>%
  na.omit() %>%
  group_by(gene, perc_sim) %>%
  summarise(cor = cor(natural, predicted)) 

data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}

# getting rid of the genes that are not present in all 5 groups
cor_wider <- cor %>%
  pivot_wider(names_from = perc_sim, values_from = cor)

cor_reduced <- na.omit(cor_wider)

cor_reduced <- cor_reduced %>%
  pivot_longer(cols =  c("(0-20%]", "(20-40%]", "(40-60%]", "(60-80%]", "(80-100%]"), names_to = "perc_sim", values_to = "cor")

plot_a <- cor_reduced %>%
  group_by(gene) %>%
  mutate(
    # pick y value corresponding to y3
    color_y = sum(cor * (perc_sim == "(80-100%]"))
  ) %>%
  ggplot(aes(x = perc_sim, y = cor, group = gene, color = color_y, fill = color_y)) +
  geom_path(size = 0.25, position = position_jitter(width = 0.05, height = 0, seed = 123)) +
  geom_point(
    shape = 21, color = "black",
    size = 2, position = position_jitter(width = 0.05, height = 0, seed = 123)
  ) +
  scale_x_discrete(
    name = "Percent Sequence Similarity of Alignment") +
  scale_y_continuous(
    name = "Correlation Coefficients",
    limits = c(-0.4, 0.6),
    breaks = seq(from = -0.4, to = 0.6, by = 0.1),
    expand = c(0, 0)) +
  scale_color_viridis_c(aesthetics = c("color", "fill"), option = "E", begin = 0.3, end = 1) +
  theme_bw() +
  theme(legend.position="none")

plot_a

#Now making a plot for classes

all_joined_wide2 <- all_joined %>%
  select(-n_eff) %>%
  pivot_wider(names_from = group, values_from = n_eff_class)

cor_2 <- all_joined_wide2 %>%
  na.omit() %>%
  group_by(gene, perc_sim) %>%
  summarise(cor = cor(natural, predicted)) 

b <- cor_2 %>%
  ggplot(aes(x = perc_sim, y = cor)) +
  geom_violin(fill = "#ecb613", alpha = 0.5) + 
  #geom_sina() +
  stat_summary(fun.data=data_summary) +
  #labs(title = "Comparing Predicted Neff to Natural Neff", 
       #subtitle = "Within Class Predictions") +
  theme_cowplot() + 
  theme(plot.title = element_text(hjust = 0.5), 
        plot.subtitle = element_text(hjust = 0.5),
        panel.grid.major.y = element_line(color = "grey92", size=0.5),
        legend.position = "none"
  ) +
  scale_x_discrete(
    name = "Percent Sequence Similarity of Alignment"
  ) +
  scale_y_continuous(
    name = "Correlation Coefficients",
    limits = c(-0.4, 0.6),
    breaks = seq(from = -0.4, to = 0.6, by = 0.1),
    expand = c(0, 0)) 

# change the plot below to show neff pred vs neff natural for classes:
b <- cor_2 %>%
  group_by(gene) %>%
  mutate(
    # pick y value corresponding to y3
    color_y = sum(cor * (perc_sim == "(80-100%]"))
  ) %>%
  ggplot(aes(x = perc_sim, y = cor, group = gene, color = color_y, fill = color_y)) +
  geom_path(size = 0.25, position = position_jitter(width = 0.05, height = 0, seed = 123)) +
  geom_point(
    shape = 21, color = "black",
    size = 2, position = position_jitter(width = 0.05, height = 0, seed = 123)
  ) +
  scale_x_discrete(
    name = "Percent Sequence Similarity of Alignment") +
  scale_y_continuous(
    name = "Correlation Coefficients",
    limits = c(-0.4, 0.6),
    breaks = seq(from = -0.4, to = 0.6, by = 0.1),
    expand = c(0, 0)) +
  scale_color_viridis_c(aesthetics = c("color", "fill"), option = "E", begin = 0.3, end = 1) +
  theme_bw() +
  theme(legend.position="none")


figure_1 <- plot_grid(a, b, nrow = 2, align="v", labels = c('A', 'B'))

ggsave(filename = "../../analysis/figures/figure_6a.png", plot = a, width = 8, height = 4)
ggsave(filename = "../../analysis/figures/figure_6b.png", plot = b, width = 8, height = 5)

#================================================================================================================
#making a plot for each class of amino acids in the wt

cnn_data <- read.csv(file = "./cnn_wt_max_freq.csv", header=TRUE, sep=",")
  
wt_classes <- cnn_data %>%
  filter(group == "wt") %>%
  select(gene, position, aa_class)

wt_labels <- inner_join(wt_classes, all_joined)

wt_labels_wide <- wt_labels %>%
  select(-n_eff_class) %>%
  pivot_wider(names_from = group, values_from = n_eff)

cor_3 <- wt_labels_wide %>%
  na.omit() %>%
  group_by(gene, perc_sim, aa_class) %>%
  summarise(cor = cor(natural, predicted)) 

cor3_wider <- cor_3 %>%
  pivot_wider(names_from = perc_sim, values_from = cor)

cor3_reduced <- na.omit(cor3_wider)

cor3_reduced <- cor3_reduced %>%
  pivot_longer(cols =  c("(0-20%]", "(20-40%]", "(40-60%]", "(60-80%]", "(80-100%]"), names_to = "perc_sim", values_to = "cor")

plot_c <- cor3_reduced %>%
  group_by(gene, aa_class) %>%
  mutate(
    # pick y value corresponding to y3
    color_y = sum(cor * (perc_sim == "(80-100%]"))
  ) %>%
  ggplot(aes(x = perc_sim, y = cor, group = gene, color = color_y, fill = color_y)) +
  geom_path(size = 0.25, position = position_jitter(width = 0.05, height = 0, seed = 123)) +
  geom_point(
    shape = 21, color = "black",
    size = 2, position = position_jitter(width = 0.05, height = 0, seed = 123)) +
  facet_wrap(vars(aa_class)) +
  scale_x_discrete(
    name = "Percent Sequence Similarity of Alignment") +
  scale_y_continuous(
    name = "Correlation Coefficients",
    limits = c(-0.4, 0.6),
    breaks = seq(from = -0.4, to = 0.6, by = 0.1),
    expand = c(0, 0)) +
  scale_color_viridis_c(aesthetics = c("color", "fill"), option = "E") +
  theme_bw() +
  theme(legend.position="none") 

plot_c

ggsave(filename = "../../analysis/figures/figure_6c.png", plot = plot_c, width = 10, height = 6)

#==============================================================================================
# SUPPLEMENTARY PLOT: boxplot of number of seqs per protein for each seq similarity group:
#================================================================================================








#================================================================================================================
#tests and extraneous code:

a <- cor_reduced %>%
  group_by(gene) %>%
  mutate(
    # pick y value corresponding to y3
    color_y = sum(cor * (perc_sim == "(80-100%]"))
  ) %>%
  ggplot(aes(x = perc_sim, y = cor, group = gene, color = color_y, fill = color_y)) +
  #geom_violin(fill = "#9875bd", alpha = 0.5) + 
  #geom_sina() +
  geom_path(size = 0.25, position = position_jitter(width = 0.05, height = 0, seed = 123)) +
  #geom_line(aes(group = gene), alpha = 0.5, size = 0.7) +
  #stat_summary(fun.data=data_summary) +
  #labs(title = "Comparing Predicted Neff to Natural Neff", 
  #     subtitle = "Amino Acid Predictions") +
  #theme(plot.title = element_text(hjust = 0.5), 
  #  plot.subtitle = element_text(hjust = 0.5),
  #  panel.grid.major.y = element_line(color = "grey92", size=0.5),
  #  legend.position = "none") +
  #theme_cowplot() +
  geom_point(
    shape = 21, color = "black",
    size = 2, position = position_jitter(width = 0.05, height = 0, seed = 123)
  ) +
  scale_x_discrete(
    name = "Percent Sequence Similarity of Alignment") +
  scale_y_continuous(
    name = "Correlation Coefficients",
    limits = c(-0.4, 0.6),
    breaks = seq(from = -0.4, to = 0.6, by = 0.1),
    expand = c(0, 0)) +
  #scale_color_discrete_qualitative(palette = "Dynamic")
  scale_color_viridis_c(aesthetics = c("color", "fill"), option = "E") +
  theme_bw() +
  theme(legend.position="none")



