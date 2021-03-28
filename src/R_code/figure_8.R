#figure 8- wt prob predicted by CNN ~ n_eff of alignment 
# done across the different % similarity groups

library(tidyverse)
library(cowplot)
library(broom)

# reading csv files
natural_var <- read.csv(file="./stats_align_all.csv", header=TRUE, sep=",")
natural_var_20 <- read.csv(file = "./stats_align_files/stats_align_20.csv", header=TRUE, sep=",")
natural_var_40 <- read.csv(file = "./stats_align_files/stats_align_40.csv", header=TRUE, sep=",")
natural_var_60 <- read.csv(file = "./stats_align_files/stats_align_60.csv", header=TRUE, sep=",")
natural_var_80 <- read.csv(file = "./stats_align_files/stats_align_80.csv", header=TRUE, sep=",")
natural_var_100 <- read.csv(file = "./stats_align_files/stats_align_100.csv", header=TRUE, sep=",")

cnn_data <- read.csv(file = "./cnn_wt_max_freq.csv", header=TRUE, sep=",")

cnn_data2 <- cnn_data %>%
  filter(group == "wt") %>%
  select(position, gene, freq) %>%
  mutate(freq_wt_cnn = freq) %>%
  select(position, gene, freq_wt_cnn)
  
#0-20
natural_var_20 <- natural_var_20 %>%
  select(position, gene, n_eff, n_eff_class) %>%
  mutate(group = "natural")

joined_20 <- inner_join(x = natural_var_20, y = cnn_data2, by = c('position', 'gene')) %>%
  select(position, gene, n_eff, n_eff_class, freq_wt_cnn) %>%
  mutate(perc_sim = "(0-20%]")

#20-40
natural_var_40 <- natural_var_40 %>%
  select(position, gene, n_eff, n_eff_class) %>%
  mutate(group = "natural")

joined_40 <- inner_join(x = natural_var_40, y = cnn_data2, by = c('position', 'gene')) %>%
  select(position, gene, n_eff, n_eff_class, freq_wt_cnn) %>%
  mutate(perc_sim = "(20-40%]") 

#40-60
natural_var_60 <- natural_var_60 %>%
  select(position, gene, n_eff, n_eff_class) %>%
  mutate(group = "natural") 

joined_60 <- inner_join(x = natural_var_60, y = cnn_data2, by = c('position', 'gene')) %>%
  select(position, gene, n_eff, n_eff_class, freq_wt_cnn) %>%
  mutate(perc_sim = "(40-60%]")  

#60-80
natural_var_80 <- natural_var_80 %>%
  select(position, gene, n_eff, n_eff_class) %>%
  mutate(group = "natural") 

joined_80 <- inner_join(x = natural_var_80, y = cnn_data2, by = c('position', 'gene')) %>%
  select(position, gene, n_eff, n_eff_class, freq_wt_cnn) %>%
  mutate(perc_sim = "(60-80%]") 

#80-100
natural_var_100 <- natural_var_100 %>%
  select(position, gene, n_eff, n_eff_class) %>%
  mutate(group = "natural") 

joined_100 <- inner_join(x = natural_var_100, y = cnn_data2, by = c('position', 'gene')) %>%
  select(position, gene, n_eff, n_eff_class, freq_wt_cnn) %>%
  mutate(perc_sim = "(80-100%]") 


all_joined <- rbind(joined_20, joined_40, joined_60, joined_80, joined_100)

#===============================================================================
# Finding correlation coefficients and p-values
#===============================================================================

# fitting a linear model to data (getting R^2 and p-values)
lm_summary <- all_joined %>%
  na.omit() %>%
  nest(data = -c(gene, perc_sim)) %>%
  mutate(
    fit = map(data, ~lm(freq_wt_cnn ~ n_eff, data = .x)),
    glance_out = map(fit, glance)
  ) %>%
  select(gene, perc_sim, glance_out) %>%
  unnest(cols = glance_out)
lm_summary

# getting rid of the genes that are not present in all 5 groups:
lm_summary <- lm_summary %>%
  select(gene, perc_sim, r.squared, p.value) 

# making dataframe for p-values:
p_values <- lm_summary %>%
  select(gene, perc_sim, p.value)

# removing genes that are not found across all 5 similarity groups:
p_values_wide <- p_values %>%
  pivot_wider(names_from = perc_sim, values_from = p.value)

p_values_reduced <- na.omit(p_values_wide)

# making both dataframes longer again for plotting:
p_values <- p_values_reduced %>%
  pivot_longer(
    cols = -gene, 
    names_to = "perc_sim", 
    values_to = c("p_value"))

# labeling the significant p-values:
p_values <- p_values %>%
  mutate(signif = ifelse(p_value <= 0.05, TRUE, FALSE))


#FINDING CORRELATION COEFFS:
cor <- all_joined %>%
  na.omit() %>%
  group_by(gene, perc_sim) %>%
  summarise(cor = cor(freq_wt_cnn, n_eff)) 

# getting rid of the genes that are not present in all 5 groups:
cor_wider <- cor %>%
  pivot_wider(names_from = perc_sim, values_from = cor)

cor_reduced <- na.omit(cor_wider)

cor_reduced <- cor_reduced %>%
  pivot_longer(cols =  c("(0-20%]", "(20-40%]", "(40-60%]", "(60-80%]", "(80-100%]"), names_to = "perc_sim", values_to = "cor")

# filtering for the correlations that are significant:
joined_cors <- inner_join(p_values, cor_reduced)

sig_cor <- joined_cors %>%
  filter(signif == TRUE) %>%
  select(cor, perc_sim, gene)

# filtering for the correlations that are **NOT** significant:
not_signif <- joined_cors %>%
  filter(signif == FALSE) %>%
  select(cor, perc_sim, gene)


#======================================================================================
# Plotting the data
#======================================================================================

sig_cor <- sig_cor %>%
  group_by(gene) %>%
  mutate(
    # pick y value corresponding to y3
    color_y = sum(cor * (perc_sim == "(80-100%]"))
  ) 

all_data <- cor_reduced %>%
  group_by(gene) %>%
  mutate(
    # pick y value corresponding to y3
    color_y = sum(cor * (perc_sim == "(80-100%]"))
  ) 


plot_8a <- ggplot() +
  geom_path(
    data = all_data,
    aes(x = perc_sim, y = cor, group = gene, color = color_y),
    size = 0.25, 
    position = position_jitter(width = 0.07, height = 0, seed = 123)) +
  geom_point(
    data = not_signif,
    aes(x = perc_sim, y = cor, group = gene),
    shape = 21, 
    color = "black",
    fill = "white",
    size = 2, 
    position = position_jitter(width = 0.07, height = 0, seed = 123)) +
  geom_point(
    data = sig_cor,
    aes(x = perc_sim, y = cor, group = gene, fill = color_y),
    shape = 21, 
    color = "black",
    size = 2, 
    position = position_jitter(width = 0.07, height = 0, seed = 123)) +
  scale_x_discrete(
    name = "% Sequence Similarity of Alignment") +
  scale_y_continuous(
    name = "Correlation Coefficients",
    #limits = c(-0.6, 0.4),
    #breaks = seq(from = -0.6, to = 0.4, by = 0.2),
    expand = c(0.01, 0.01)) +
  scale_color_gradient(
    aesthetics = c("color", "fill"), 
    high = "#ffd966", 
    low = "#080845") +
  theme_bw(12) +
  theme(
    legend.position = "none",
    axis.text = element_text(color = "black", size = 12),
    panel.grid.minor = element_blank())

plot_8a
ggsave(filename = "../../analysis/figures/figure_8a.png", plot = plot_8a, width = 8, height = 4)


#===============================================================================
# a jitter plot per wt class
#===============================================================================

cnn_data <- read.csv(file = "./cnn_wt_max_freq.csv", header=TRUE, sep=",")

wt_classes <- cnn_data %>%
  filter(group == "wt") %>%
  select(gene, position, aa_class)

wt_labels <- inner_join(wt_classes, all_joined)

wt_labels_wide <- wt_labels %>%
  select(-n_eff_class)

# getting the linear model:
# fitting a linear model to data (getting R^2 and p-values)
lm_summary_3 <- wt_labels_wide %>%
  na.omit() %>%
  nest(data = -c(gene, perc_sim, aa_class)) %>%
  mutate(
    fit = map(data, ~lm(freq_wt_cnn ~ n_eff, data = .x)),
    glance_out = map(fit, glance)
  ) %>%
  select(gene, perc_sim, aa_class, glance_out) %>%
  unnest(cols = glance_out)
lm_summary_3

# getting rid of the genes that are not present in all 5 groups:
p_values_3 <- lm_summary_3 %>%
  select(gene, perc_sim, aa_class, p.value) 


# removing genes that are not found across all 5 similarity groups:
p_values_wide <- p_values_3 %>%
  pivot_wider(names_from = perc_sim, values_from = p.value)

p_values_reduced <- na.omit(p_values_wide)

# making both dataframes longer again for plotting:
p_values_3 <- p_values_reduced %>%
  pivot_longer(
    cols = -c(gene, aa_class),
    names_to = "perc_sim", 
    values_to = c("p_value"))

# labeling the significant p-values:
p_values_3 <- p_values_3 %>%
  mutate(signif = ifelse(p_value <= 0.05, TRUE, FALSE))

#getting the correlation coefficients
cor_3 <- wt_labels_wide %>%
  na.omit() %>%
  group_by(gene, perc_sim, aa_class) %>%
  summarise(cor = cor(freq_wt_cnn, n_eff)) 

cor3_wider <- cor_3 %>%
  pivot_wider(names_from = perc_sim, values_from = cor)

cor3_reduced <- na.omit(cor3_wider)

cor3_reduced <- cor3_reduced %>%
  pivot_longer(cols =  c("(0-20%]", "(20-40%]", "(40-60%]", "(60-80%]", "(80-100%]"), names_to = "perc_sim", values_to = "cor")


# filtering for the correlations that are significant:
joined_cors <- inner_join(p_values_3, cor3_reduced)

sig_cor <- joined_cors %>%
  filter(signif == TRUE) %>%
  select(cor, perc_sim, gene, aa_class)

# filtering for the correlations that are **NOT** significant:
not_signif <- joined_cors %>%
  filter(signif == FALSE) %>%
  select(cor, perc_sim, gene, aa_class)

#getting values for coloring the significant dots and all lines:
sig_cor <- sig_cor %>%
  group_by(gene, aa_class) %>%
  mutate(
    # pick y value corresponding to y3
    color_y = sum(cor * (perc_sim == "(80-100%]"))
  ) 

all_data <- cor3_reduced %>%
  group_by(gene, aa_class) %>%
  mutate(
    # pick y value corresponding to y3
    color_y = sum(cor * (perc_sim == "(80-100%]"))
  ) 


plot_8c <- ggplot() +
  geom_path(
    data = all_data,
    aes(x = perc_sim, y = cor, group = gene, color = color_y),
    size = 0.25, 
    position = position_jitter(width = 0.07, height = 0, seed = 123)) +
  geom_point(
    data = not_signif,
    aes(x = perc_sim, y = cor, group = gene),
    shape = 21, 
    color = "black",
    fill = "white",
    size = 2, 
    position = position_jitter(width = 0.07, height = 0, seed = 123)) +
  geom_point(
    data = sig_cor,
    aes(x = perc_sim, y = cor, group = gene, fill = color_y),
    shape = 21, 
    color = "black",
    size = 2, 
    position = position_jitter(width = 0.07, height = 0, seed = 123)) +
  facet_wrap(vars(aa_class)) +
  scale_x_discrete(
    name = "% Sequence Similarity of Alignment") +
  scale_y_continuous(
    name = "Correlation Coefficients",
    limits = c(-0.4, 0.6),
    breaks = seq(from = -0.4, to = 0.6, by = 0.1),
    expand = c(0, 0)) +
  scale_color_gradient(
    aesthetics = c("color", "fill"), 
    high = "#ffd966", 
    low = "#080845") +
  theme_bw(14) +
  theme(
    legend.position = "none",
    axis.text = element_text(color = "black", size = 14),
    panel.grid.minor = element_blank(),
    strip.text.x = element_text(size = 16))

plot_8c

ggsave(filename = "../../analysis/figures/figure_8c.png", plot = plot_8c, width = 13.5, height = 6)




