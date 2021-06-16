# continuation of figure 2. 
# comparing neff predicted vs. neff natural across different alignments similarities. 
library(tidyverse)
library(cowplot)
library(broom)
library(stats)

box_size = "20"

# useful function for getting mean and standard of deviation (for violin plots):
data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}

#set working directory to:
#C:\Users\avch\Desktop\Natural_var_project

# reading csv files
natural_var <- read.csv(file="./output/output_PSICOV/stats_align_all.csv", header=TRUE, sep=",")
natural_var_20 <- read.csv(file = "./output/output_PSICOV/stats_align_files/stats_align_20.csv", header=TRUE, sep=",")
natural_var_40 <- read.csv(file = "./output/output_PSICOV/stats_align_files/stats_align_40.csv", header=TRUE, sep=",")
natural_var_60 <- read.csv(file = "./output/output_PSICOV/stats_align_files/stats_align_60.csv", header=TRUE, sep=",")
natural_var_80 <- read.csv(file = "./output/output_PSICOV/stats_align_files/stats_align_80.csv", header=TRUE, sep=",")
natural_var_100 <- read.csv(file = "./output/output_PSICOV/stats_align_files/stats_align_100.csv", header=TRUE, sep=",")

cnn_var <- read.csv(file= paste0("./data/PSICOV_box_",box_size,"/output/stats_cnn.csv"), header=TRUE, sep=",")
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

#===============================================================================
# Finding correlation coefficients and p-values
#===============================================================================

# fitting a linear model to data (getting R^2 and p-values)
lm_summary <- all_joined_wide %>%
  na.omit() %>%
  nest(data = -c(gene, perc_sim)) %>%
  mutate(
    fit = map(data, ~lm(natural ~ predicted, data = .x)),
    glance_out = map(fit, glance)
    ) %>%
  select(gene, perc_sim, glance_out) %>%
  unnest(cols = glance_out)
lm_summary

# getting rid of the genes that are not present in all 5 groups:
lm_summary <- lm_summary %>%
  select(gene, perc_sim, r.squared, p.value) 

# making two dataframes for correlation coefficients and p-values:
cor_coeffs <- lm_summary %>%
  select(gene, perc_sim, r.squared)

p_values <- lm_summary %>%
  select(gene, perc_sim, p.value)
p_values

# false discovery correction is applied to the p-values:
new_p_values <- p_values %>%
  nest(data = -c(perc_sim)) %>%
  mutate(
    new_p = map(data, ~p.adjust(.x$p.value, method = "fdr", n = length(.x$p.value)))
  ) %>%
  unnest(cols = c(data, new_p))



# removing genes that are not found across all 5 similarity groups:
cor_coeffs_wide <- cor_coeffs %>%
  pivot_wider(names_from = perc_sim, values_from = r.squared)

p_values_wide <- new_p_values %>%
  select(perc_sim, gene, new_p) %>%
  pivot_wider(names_from = perc_sim, values_from = new_p)

cor_coeffs_reduced <- na.omit(cor_coeffs_wide)
p_values_reduced <- na.omit(p_values_wide)

# making both dataframes longer again for plotting:
cor_coeffs <- cor_coeffs_reduced %>%
  pivot_longer(
    cols = -gene, 
    names_to = "perc_sim", 
    values_to = c("r_squared"))

p_values <- p_values_reduced %>%
  pivot_longer(
    cols = -gene, 
    names_to = "perc_sim", 
    values_to = c("p_value"))

# labeling the significant p-values:
p_values <- p_values %>%
  mutate(signif = ifelse(p_value <= 0.05, TRUE, FALSE))



#FINDING CORRELATION COEFFS:
# an alternative method for finding the correlations:
cor <- all_joined_wide %>%
  na.omit() %>%
  group_by(gene, perc_sim) %>%
  summarise(cor = cor(natural, predicted)) 

# getting rid of the genes that are not present in all 5 groups:
cor_wider <- cor %>%
  pivot_wider(names_from = perc_sim, values_from = cor)

cor_reduced <- na.omit(cor_wider)

cor_reduced <- cor_reduced %>%
  pivot_longer(cols =  c("(0-20%]", "(20-40%]", "(40-60%]", "(60-80%]", "(80-100%]"), names_to = "perc_sim", values_to = "cor")

# filtering for the correlations that are significant:
# you need this for the figure_6_box_size.R script
joined_cors_1 <- inner_join(p_values, cor_reduced)

get_x_value <- function(perc_sim){
  if (perc_sim == "(0-20%]") {
    x_value <- as.numeric(1)
  } else if (perc_sim == "(20-40%]") {
    x_value <- as.numeric()
  } else if (perc_sim == "(40-60%]") {
    x_value <- as.numeric(1)
  } else if (perc_sim == "(60-80%]") {
    x_value <- as.numeric(1)
  } else {
    x_value <- as.numeric(1)
  }
  return(x_value)
}

joined_cors_1 <- joined_cors_1 %>%
  group_by(gene) %>%
  mutate(
    # pick y value corresponding to y3
    color_y = sum(cor * (perc_sim == "(40-60%]")),
    dx = rnorm(n(), mean = 0, sd = .05),
    dy = rnorm(n(), mean = 0, sd = .05),
    x_value = as.numeric(factor(perc_sim)))  

sig_cor <- joined_cors_1 %>%
  filter(signif == TRUE) %>%
  select(cor, perc_sim, gene, dx, dy, color_y)



get_dx_dy <- function(perc_sim){
  return(rnorm(n(), mean = map(perc_sim, get_x_value), sd = 0.1))
}

# filtering for the correlations that are **NOT** significant:
not_signif <- joined_cors_1 %>%
  filter(signif == FALSE) %>%
  select(cor, perc_sim, gene, dx, dy, color_y)
  
not_signif
  


#===========================================================================================
# plot_a (correlations bw neff predicted and neff natural across seq similarity groups)
#===========================================================================================

# just in case (not used below)
# all_data <- cor_reduced %>%
#   group_by(gene) %>%
#   mutate(
#     # pick y value corresponding to y3
#     color_y = sum(cor * (perc_sim == "(80-100%]")),
#     dx = rnorm(n(), mean = 0, sd = .1),
#     dy = rnorm(n(), mean = 0, sd = .1),
#     x_value = as.numeric(factor(perc_sim))) 
#  
# all_data

plot_a <- ggplot() +
  geom_path(
    data = not_signif,
    aes(x = as.numeric(factor(perc_sim))+dx, y = cor+dy, group = gene, color = color_y),
    size = 0.25) +
  geom_path(
    data = sig_cor,
    aes(x = as.numeric(factor(perc_sim))+dx, y = cor+dy, group = gene, color = color_y),
    size = 0.25) +
  geom_point(
    data = not_signif,
    aes(x = as.numeric(factor(perc_sim))+dx, y = cor+dy, group = gene),
    shape = 21, 
    color = "black",
    fill = "white",
    size = 2) +
  geom_point(
    data = sig_cor,
    aes(x = as.numeric(factor(perc_sim))+dx, y = cor+dy, group = gene, fill = color_y),
    shape = 21, 
    color = "black",
    size = 2) +
  scale_x_continuous(
    name = "Percent sequence similarity to wild type",
    limits = c(0.5,5.5),
    labels = c("(0-20%]", "(20-40%]", "(40-60%]", "(60-80%]", "(80-100%]"),
    breaks = (seq(from = 1, to = 5, by = 1))) +
  scale_y_continuous(
    name = "Correlation Coefficients",
    limits = c(-0.25, 0.7),
    breaks = seq(from = -0.2, to = 0.6, by = 0.1),
    expand = c(0, 0)) +
  scale_color_gradient(
    aesthetics = c("color", "fill"), 
    high = "#ffd966", 
    low = "#080845") +
  theme_bw(16) +
  theme(
    legend.position = "none",
    axis.text = element_text(color = "black", size = 16),
    panel.grid.minor = element_blank())

plot_a


ggsave(filename = paste0("./analysis/figures/figure_6a_box_",box_size,".png"), plot = plot_a, width = 8, height = 4)

#===============================================================================================
#plot_b (Now making a plot for neff natural CLASSES vs neff predicted CLASSES)
#===============================================================================================

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
plot_b <- cor_2 %>%
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
    name = "Correlation Coefficients \n",
    limits = c(-0.4, 0.6),
    breaks = seq(from = -0.4, to = 0.6, by = 0.1),
    expand = c(0, 0)) +
  scale_color_viridis_c(aesthetics = c("color", "fill"), option = "E", begin = 0.3, end = 1) +
  theme_bw() +
  theme(
    legend.position="none")
plot_b

figure_1 <- plot_grid(a, b, nrow = 2, align="v", labels = c('A', 'B'))


ggsave(filename = paste0("./analysis/figures/figure_6b_box_",box_size,".png"), plot = plot_b, width = 8, height = 5)

#==================================================================================================
#================================================================================================================
# plot_c (making a plot for each class of amino acids in the wt)
#================================================================================================================
#=====================================================================================================

cnn_data <- read.csv(file = "./data/PSICOV_box_20/output/cnn_wt_max_freq.csv", header=TRUE, sep=",")

wt_classes <- cnn_data %>%
  filter(group == "wt") %>%
  select(gene, position, aa_class)

wt_labels <- inner_join(wt_classes, all_joined)

wt_labels <- wt_labels %>%
  mutate(aa_class = fct_recode(aa_class, `small polar` = "small_polar"))

wt_labels_wide <- wt_labels %>%
  select(-n_eff_class) %>%
  pivot_wider(names_from = group, values_from = n_eff)

# getting the linear model:
# fitting a linear model to data (getting R^2 and p-values)
lm_summary_3 <- wt_labels_wide %>%
  na.omit() %>%
  nest(data = -c(gene, perc_sim, aa_class)) %>%
  mutate(
    fit = map(data, ~lm(natural ~ predicted, data = .x)),
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
  summarise(cor = cor(natural, predicted)) 

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

scaleFUN <- function(x) sprintf("%.1f", x)

plot_c <- ggplot() +
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
  facet_wrap(vars(aa_class), ncol = 2) +
  scale_x_discrete(
    name = "Percent sequence similarity to wild type") +
  scale_y_continuous(
    name = "Correlation Coefficients",
    limits = c(-0.6, 1.0),
    breaks = seq(from = -0.6, to = 1.0, by = 0.2),
    expand = c(0.02, 0.02),
    labels = scaleFUN) +
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

plot_c

ggsave(filename = paste0("./analysis/figures/figure_6c_box_20.png"), plot = plot_c, width = 10, height = 9)

#==============================================================================================
# SUPPLEMENTARY PLOT: boxplot of number of seqs per protein for each seq similarity group:
#================================================================================================

seq_counts <- read.csv(file = "./output/output_PSICOV/seq_counts.csv", header=TRUE, sep=",")

plot_d <- seq_counts %>%
ggplot(aes(x = factor(group), y = seq_count)) +
  geom_boxplot(fill = "grey86") +
  scale_x_discrete(
    name = "Percent sequence similarity to wild type",
    breaks = c(20, 40, 60, 80, 100),
    labels = c("(0-20%]", "(20-40%]", "(40-60%]", "(60-80%]", "(80-100%]")) + 
  scale_y_log10(
    name = "Sequence Count per Protein",
    breaks = c(0, 10, 100, 1000, 10000, 60000),
    labels = c("0", "10", "100", "1000", "10000", "60000")) +
  theme_cowplot(14)+
  theme(
    axis.text = element_text(color = "black", size = 14),
    strip.text.x = element_text(size = 16)
  )

plot_d
  
ggsave(filename = "./analysis/figures/aln_seq_count.png", plot = plot_d, width = 7.5, height = 6)




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


#=======================================================================================
# t-tests for plor a
#=======================================================================================

for_ttests <- sig_cor %>%
  select(c(gene, perc_sim, cor))

averages <- for_ttests %>%
  group_by(perc_sim) %>%
  summarise(mean = mean(cor))
averages

for_paired_test <- for_ttests %>%
  pivot_wider(names_from = perc_sim, values_from = cor) %>%
  rename("to20" = `(0-20%]`,
         "to40" = `(20-40%]`,
          "to60" = `(40-60%]`,
          "to80" = `(60-80%]`,
         "to100" = `(80-100%]`)


results_1 <- t.test(for_paired_test$to20, for_paired_test$to40, paired = TRUE, alternative = "two.sided")
results_1
#p-value = 5.825e-07 (there is a significant difference)

results_2 <- t.test(for_paired_test$to40, for_paired_test$to60, paired = TRUE, alternative = "two.sided")
results_2
#p-value = 0.1257 (no difference)

results_3 <- t.test(for_paired_test$to60, for_paired_test$to80, paired = TRUE, alternative = "two.sided")
results_3
#p-value = 0.7802 (no difference)

results_4 <- t.test(for_paired_test$to80, for_paired_test$to100, paired = TRUE, alternative = "two.sided")
results_4
#p-value = 0.0025 (there is a significant difference)

results_5 <- t.test(for_paired_test$to40, for_paired_test$to80, paired = TRUE, alternative = "two.sided")
results_5
#p-value = 0.3758 (no difference)

results_6 <- t.test(for_paired_test$to20, for_paired_test$to60, paired = TRUE, alternative = "two.sided")
results_6
#p-value = 0.001395 (there is a significant difference)

results_7 <- t.test(for_paired_test$to60, for_paired_test$to100, paired = TRUE, alternative = "two.sided")
results_7
#p-value = 0.001549 (there is a significant difference)





