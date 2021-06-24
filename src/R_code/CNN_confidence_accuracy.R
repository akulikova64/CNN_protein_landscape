library(tidyverse)
library(cowplot)
library(ggforce)

#===============================================================================
# Analysis for protein engineering applications.
#===============================================================================
# Looking at wt/consensus accuracy as a function of CNN confidence (pred_prob or n-eff)
#===============================================================================
# ***BOX SIZE 20*****
#===============================================================================

# set working directory to: "Desktop/Natural_var_project/"
# loading data
cnn_data <- read.csv(file = "./data/PSICOV_box_20/output/cnn_wt_max_freq.csv", header=TRUE, sep=",")
natural_data <- read.csv(file = "./data/PSICOV_box_20/output/natural_max_freq_files/natural_max_freq_all.csv", header=TRUE, sep=",")

joined_data <- rbind(x = cnn_data, y = natural_data)

joined_data_trimmed <- joined_data %>%
  filter(!gene %in% c('1dbx', '1eaz', '1fvg', '1k7j', '1kq6', '1kw4', '1lpy', '1ne2', '1ny1', '1pko', '1rw1', '1vhu', '1w0h', '1wkc', '2tps'))

wide <- joined_data_trimmed %>%
  select(-c(aa_class, class_freq)) %>%
  pivot_wider(names_from = group, values_from = c(freq, aa)) 


get_pred_bin <- function(x) {
  
  if (x > 0 & x <= 0.2) {
    return("(0-0.2]")
  }
  else if (x > 0.2 & x <= 0.4) {
    return("(0.2-0.4]")
  }
  else if (x > 0.4 & x <= 0.6) {
    return("(0.4-0.6]")
  }
  else if (x > 0.6 & x <= 0.8) {
   return("(0.6-0.8]")
  }
  else if (x > 0.8 & x <= 1.0) {
    return("(0.8-1.0]")
  }
}


wide_new <- wide %>%
  na.omit() %>%
  mutate(pred_bin = map_chr(freq_predicted, get_pred_bin))
 

matches <- wide_new %>%
  mutate(match_predict_cons = aa_predicted == aa_natural_max) %>%
  mutate(match_predict_wt = aa_predicted == aa_wt)


for_plot_a <- matches %>%
  select(c(gene, position, pred_bin, match_predict_cons )) %>%
  group_by(gene, pred_bin) %>%
  summarise(freq_predict_cons = sum(match_predict_cons, na.rm = TRUE)/sum(!is.na(match_predict_cons)))


stat_data_1 <- for_plot_a %>%
  select(-gene) %>%
  group_by(pred_bin) %>%
  summarise(estimate = mean(freq_predict_cons),
            std_error = sd(freq_predict_cons)/sqrt(length(freq_predict_cons)))


plot_a <- for_plot_a %>%
  ggplot(aes(y = freq_predict_cons, x = pred_bin)) +
  geom_violin(alpha = 0.5, size = 0.7, bw = 0.03, fill = "#8c7b9d", color = "#655775") + 
  geom_pointrange(data = stat_data_1, aes(x = pred_bin,
                                          y = estimate,
                                          ymin = estimate - 1.96*std_error,
                                          ymax = estimate + 1.96*std_error),
                  color = "black", alpha = 0.7, size = 0.3) +
  #stat_summary(fun.data=data_summary, color = "black", alpha = 0.7) +
  theme_cowplot(16) + 
  theme(plot.title = element_text(hjust = 0, size = 16), 
        plot.subtitle = element_text(hjust = 0.5),
        panel.grid.major.y = element_line(color = "grey92", size=0.5),
        legend.position = "none") +
  scale_y_continuous(
    name = "Accuracy",
    limits = c(0, 1.0),
    breaks = seq(from = 0, to = 1.0, by = 0.1),
    expand = c(0, 0)) +
  scale_x_discrete(
    name = "CNN confidence")

plot_a

# Now, looking at wt predictions.

for_plot_b <- matches %>%
  select(c(gene, position, pred_bin, match_predict_wt )) %>%
  group_by(gene, pred_bin) %>%
  summarise(freq_predict_wt = sum(match_predict_wt, na.rm = TRUE)/sum(!is.na(match_predict_wt)))

stat_data_2 <- for_plot_b %>%
  select(-gene) %>%
  group_by(pred_bin) %>%
  summarise(estimate = mean(freq_predict_wt),
            std_error = sd(freq_predict_wt)/sqrt(length(freq_predict_wt)))

plot_b <- for_plot_b %>%
  ggplot(aes(y = freq_predict_wt, x = pred_bin)) +
  geom_violin(alpha = 0.5, size = 0.7, bw = 0.03, fill = "#8c7b9d", color = "#655775") + 
  geom_pointrange(data = stat_data_2, aes(x = pred_bin,
                                        y = estimate,
                                        ymin = estimate - 1.96*std_error,
                                        ymax = estimate + 1.96*std_error),
                  color = "black", alpha = 0.7, size = 0.3) +
  #stat_summary(fun.data=data_summary, color = "black", alpha = 0.7) +
  theme_cowplot(16) + 
  theme(plot.title = element_text(hjust = 0, size = 16), 
        plot.subtitle = element_text(hjust = 0.5),
        panel.grid.major.y = element_line(color = "grey92", size=0.5),
        legend.position = "none") +
  scale_y_continuous(
    name = "Accuracy",
    limits = c(0, 1.0),
    breaks = seq(from = 0, to = 1.0, by = 0.1),
    expand = c(0, 0)) +
  scale_x_discrete(
    name = "CNN confidence")

plot_b

acc_vs_conf <- plot_grid(plot_b, plot_a, nrow = 2, align = "h", labels = c('a', 'b'))
ggsave(filename = paste0("./analysis/figures/acc_vs_conf.png"), plot = acc_vs_conf, width = 8, height = 8)


ggsave(filename = "./analysis/figures/cons_acc_vs_pred.png", plot = plot_a, width = 7, height = 5)
ggsave(filename = "./analysis/figures/wt_acc_vs_pred.png", plot = plot_b, width = 7, height = 5)

#==========================================================================================
# Dodging the data from plots a and b into one plot
#==========================================================================================

freqs <- matches %>%
  select(c(gene, position, pred_bin, match_predict_wt, match_predict_cons)) %>%
  group_by(gene, pred_bin) %>%
  summarise(freq_predict_wt = sum(match_predict_wt, na.rm = TRUE)/sum(!is.na(match_predict_wt)),
            freq_predict_cons = sum(match_predict_cons, na.rm = TRUE)/sum(!is.na(match_predict_cons)))

matches_long <- freqs %>%
  pivot_longer(cols = c(freq_predict_cons, freq_predict_wt), names_to = "match_group", values_to = "freq") 
  
matches_long_2 <- matches_long %>%
  mutate(prediction = ifelse(match_group == "freq_predict_cons", "consensus", "wild type"))

custom_fills <- c("#9cbfe2", "#dbb070")
custom_colors <- c("#28598a", "#8f6424")

plot_c <- matches_long_2 %>%
  ggplot(aes(y = freq, x = pred_bin, fill = prediction, color = prediction)) +
  geom_sina(alpha = 0.6, size = 0.7, position = position_dodge(width = 0.6)) + 
  stat_summary(fun.data=data_summary, color = "black", alpha = 0.7, position = position_dodge(width = 0.6)) +
  theme_cowplot(14) + 
  theme(plot.title = element_text(hjust = 0, size = 14), 
        plot.subtitle = element_text(hjust = 0.5),
        panel.grid.major.y = element_line(color = "grey92", size=0.5),
        legend.position = "top",
        legend.direction = "vertical") +
  scale_fill_manual(values = custom_fills) +
  scale_color_manual(values = custom_colors) +
  scale_y_continuous(
    name = "Accuracy",
    limits = c(0.0, 1.0),
    breaks = seq(0, 1.0, by = 0.1),
    expand = c(0, 0)) +
  scale_x_discrete(
    name = "Predicted Probability")


plot_c

ggsave(filename = "./analysis/figures/sina_CNN_conf.png", plot = plot_c, width = 9, height = 5)


# boxplots of positions per protein within each prediction prob bin.

for_box_plot <- matches %>%
  select(gene, position, pred_bin) %>%
  group_by(gene, pred_bin) %>%
  count()
  
  
box_plot <- for_box_plot %>%
  ggplot(aes(x = factor(pred_bin), y = n)) +
  geom_boxplot(fill = "grey86") +
  scale_x_discrete(
    name = "Predicted probability") + 
  scale_y_continuous(
    name = "Number of positions per protein",
    limits = c(0, 130),
    breaks = seq(0, 130, by = 10),
    expand = c(0, 0)) +
  theme_cowplot(14)+
  theme(
    axis.text = element_text(color = "black", size = 14),
    strip.text.x = element_text(size = 16),
    panel.grid.major.y = element_line(color = "grey92", size=0.5)
  )

box_plot

ggsave(filename = "./analysis/figures/positions_per_bin.png", plot = box_plot, width = 9, height = 5)


#===============================================================================================
# Looking how predicted amino acid distributions relate to CNN confidence 
#===============================================================================================

fills <- c("#990008", "#0a2575", "#b35900", "#1a6600", "#5c0679", "#9e9e2e")

calc_class <- function(x) {
  aliphatic = c("M", "L", "I", "V", "A")
  small_polar = c("C", "S", "T", "N", "Q")
  negative = c("D", "E")
  positive = c("R", "K")
  aromatic = c("H", "Y", "W", "F")
  unique = c("P", "G")
  
  if (x %in% aliphatic) {
    return("aliphatic")
  }
  if (x %in% small_polar) {
    return("small_polar")
  }
  if (x %in% negative) {
    return("negative")
  }
  if (x %in% positive) {
    return("positive")
  }
  if (x %in% aromatic) {
    return("aromatic")
  }
  if (x %in% unique) {
    return("unique")
  }
  return("not found")
}

get_pred_bin <- function(x) {
  
  if (x > 0.80 & x <= 1.0) {
    return("Predicted at 80-100% confidence")
  }
  else{
    return(NA)
  }
}

# set working directory to: "Desktop/Natural_var_project/"
# loading data
cnn_data <- read.csv(file = "./data/PSICOV_box_20/output/cnn_wt_max_freq.csv", header=TRUE, sep=",")
natural_data <- read.csv(file = "./data/PSICOV_box_20/output/natural_max_freq_files/natural_max_freq_all.csv", header=TRUE, sep=",")

joined_data <- rbind(x = cnn_data, y = natural_data)

joined_data_trimmed <- joined_data %>%
  filter(!gene %in% c('1dbx', '1eaz', '1fvg', '1k7j', '1kq6', '1kw4', '1lpy', '1ne2', '1ny1', '1pko', '1rw1', '1vhu', '1w0h', '1wkc', '2tps'))

wide <- joined_data_trimmed %>%
  select(-c(aa_class, class_freq)) %>%
  pivot_wider(names_from = group, values_from = c(freq, aa)) 

wide_new <- wide %>%
  na.omit() %>%
  mutate(pred_bin = map_chr(freq_predicted, get_pred_bin)) %>%
  na.omit()

new_data <- wide_new %>%
  select(c(aa_predicted, pred_bin))

# finds the count of each aa acid per bin:
aa_counts <- new_data %>%
  group_by(pred_bin) %>%
  count(aa_predicted) %>%
  mutate(aa_count = n) %>%
  select(-n) %>%
  ungroup()

# getting aa fractions in the wild type structures (# alanines/ # total sites)
wt_fract <- wide %>%
  select(gene, position, aa_wt) %>%
  na.omit() %>%
  group_by(aa_wt) %>%
  count() %>%
  ungroup() %>%
  #mutate(sum = sum(n))
  #mutate(fract = n/sum) %>%
  #select(-c(n, sum)) %>%
  rename(aa_predicted = aa_wt) #renaming to an incorrect name so that we can do a join later

# normalizing the aa count per bin: 
# multiplying the aa count per bin by the total fraction of that amino acid in the wt structures:
aa_counts_norm <- aa_counts %>%
  left_join(wt_fract) %>%
  mutate(norm_count = aa_count/n) %>%
  select(-c(aa_count, n))
  

# now I need to add up all aa within each (normaized!!!) bin to get bin totals and append this to the aa_counts
for_barplot <- aa_counts_norm %>%
  group_by(pred_bin) %>%
  mutate(bin_count = sum(norm_count)) %>%
  ungroup()

# adding the class labels
with_classes <- for_barplot %>%
  mutate(class = map_chr(aa_predicted, calc_class))

for_barplot_2 <- with_classes %>%
  mutate(freq = norm_count/bin_count) %>%
  select(-c(norm_count, bin_count)) %>%
  mutate(aa_predicted = fct_rev(fct_relevel(aa_predicted, "G", "L", "P", "I", "W", "F", "C", "V", "D", "A", "M", "T", "S", "Y", "E", "N", "R", "H", "K", "Q"))) %>%
  mutate(class = fct_relevel(class, "aliphatic", "small_polar", "negative", "positive", "aromatic", "unique"))


# figuring out the order for the first facet:
order <- for_barplot_2 %>%
  filter(pred_bin == "Predicted at 80-100% confidence")

plot_e <- for_barplot_2 %>%
  ggplot(aes(x = freq, y = aa_predicted, fill = class)) +
  geom_col(alpha = 0.75) +
  #facet_wrap(vars(fct_relevel(pred_bin, "Predicted at 80-100% confidence", "Predicted at 0-20% confidence")), ncol = 1) +
  scale_fill_manual(
    values = fills,
    labels = c("aliphatic", "small polar", "negative", "positive", "aromatic", "unique")) +
  scale_x_continuous(
    name = "Frequency",
    limits = c(0.0, 0.165),
    breaks = seq(0.0, 0.16, by = 0.04),
    expand = c(0, 0)) + 
  scale_y_discrete(
    name = "Predicted amino acid",
    expand = c(0.03, 0.03)) + 
  theme_cowplot(16) +
  theme(
    axis.text = element_text(color = "black", size = 14),
    strip.text.x = element_text(size = 16),
    panel.grid.major.x = element_line(color = "grey92", size=0.5),
    panel.grid.minor.x = element_line(color = "grey92", size=0.5),
    panel.spacing = unit(2, "lines"),
    legend.position = "none")

plot_e

ggsave(filename = "./analysis/figures/Pred_aa_freq.png", plot = plot_e, width = 8, height = 8)

#===============================================================================================
#  training data amino acid distributions (b)
#===============================================================================================
calc_class <- function(x) {
  aliphatic = c("M", "L", "I", "V", "A")
  small_polar = c("C", "S", "T", "N", "Q")
  negative = c("D", "E")
  positive = c("R", "K")
  aromatic = c("H", "Y", "W", "F")
  unique = c("P", "G")
  
  if (x %in% aliphatic) {
    return("aliphatic")
  }
  if (x %in% small_polar) {
    return("small_polar")
  }
  if (x %in% negative) {
    return("negative")
  }
  if (x %in% positive) {
    return("positive")
  }
  if (x %in% aromatic) {
    return("aromatic")
  }
  if (x %in% unique) {
    return("unique")
  }
  return("not found")
}




training_data <- read.csv(file = "./data/training_data.csv", header=TRUE, sep=",")

freqs <- training_data %>%
  mutate(sum = sum(count)) %>%
  mutate(freq = count/sum) %>%
  mutate(class = map_chr(aa, calc_class)) %>%
  mutate(aa = fct_rev(fct_relevel(aa, "G", "L", "P", "I", "W", "F", "C", "V", "D", "A", "M", "T", "S", "Y", "E", "N", "R", "H", "K", "Q"))) %>%
  mutate(class = fct_relevel(class, "aliphatic", "small_polar", "negative", "positive", "aromatic", "unique"))


plot_train <- freqs %>%
  ggplot(aes(x = freq, y = aa, fill = class)) +
  geom_col(alpha = 0.75) +
  #facet_wrap(vars(fct_relevel(nat_bin, "Natural Frequency of 80-100%", "Natural Frequency of 0-20%")), ncol = 1) +
  scale_fill_manual(
    values = fills,
    labels = c("aliphatic", "small polar", "negative", "positive", "aromatic", "unique")) +
  scale_x_continuous(
    name = "Frequency",
    limits = c(0.0, 0.165),
    breaks = seq(0.0, 0.16, by = 0.04),
    expand = c(0, 0)) + 
  scale_y_discrete(
    name = "Amino acid in training data",
    expand = c(0.03, 0.03)) + 
  theme_cowplot(16) +
  theme(
    axis.text = element_text(color = "black", size = 14),
    strip.text.x = element_text(size = 16),
    panel.grid.major.x = element_line(color = "grey92", size=0.5),
    panel.grid.minor.x = element_line(color = "grey92", size=0.5),
    panel.spacing = unit(2, "lines"))

plot_train

ggsave(filename = "./analysis/figures/train_aa_freq.png", plot = plot_train, width = 8, height = 8)


#===============================================================================================
#  wild type amino acid distributions within their MSA frequency bins
#===============================================================================================

calc_class <- function(x) {
  aliphatic = c("M", "L", "I", "V", "A")
  small_polar = c("C", "S", "T", "N", "Q")
  negative = c("D", "E")
  positive = c("R", "K")
  aromatic = c("H", "Y", "W", "F")
  unique = c("P", "G")
  
  if (x %in% aliphatic) {
    return("aliphatic")
  }
  if (x %in% small_polar) {
    return("small_polar")
  }
  if (x %in% negative) {
    return("negative")
  }
  if (x %in% positive) {
    return("positive")
  }
  if (x %in% aromatic) {
    return("aromatic")
  }
  if (x %in% unique) {
    return("unique")
  }
  return("not found")
}

get_pred_bin <- function(x) {
  
  if (x > 0.80 & x <= 1.0) {
    return("Natural Frequency of 80-100%")
  }
  else{
    return(NA)
  }
}

# set working directory to: "Desktop/Natural_var_project/"
# loading data
cnn_data <- read.csv(file = "./data/PSICOV_box_20/output/cnn_wt_max_freq.csv", header=TRUE, sep=",")
natural_data <- read.csv(file = "./data/PSICOV_box_20/output/natural_max_freq_files/natural_max_freq_all.csv", header=TRUE, sep=",")

joined_data <- rbind(x = cnn_data, y = natural_data)

joined_data_trimmed <- joined_data %>%
  filter(!gene %in% c('1dbx', '1eaz', '1fvg', '1k7j', '1kq6', '1kw4', '1lpy', '1ne2', '1ny1', '1pko', '1rw1', '1vhu', '1w0h', '1wkc', '2tps'))

wide <- joined_data_trimmed %>%
  select(-c(aa_class, class_freq)) %>%
  pivot_wider(names_from = group, values_from = c(freq, aa)) 

wide_new <- wide %>%
  na.omit() %>%
  mutate(nat_bin = map_chr(freq_natural_wt, get_pred_bin)) %>%
  na.omit()

new_data <- wide_new %>%
  select(c(aa_wt, nat_bin))

# finds the count of each aa acid per bin:
aa_counts <- new_data %>%
  group_by(nat_bin) %>%
  count(aa_wt) %>%
  mutate(aa_count = n) %>%
  select(-n) %>%
  ungroup()

# getting aa fractions in the wild type structures (# alanines/ # total sites)
wt_fract <- wide %>%
  select(gene, position, aa_wt) %>%
  na.omit() %>%
  group_by(aa_wt) %>%
  count() %>%
  ungroup() 
  #mutate(sum = sum(n))
  #mutate(fract = n/sum) %>%
  #select(-c(n, sum)) %>%

# normalizing the aa count per bin: 
# multiplying the aa count per bin by the total fraction of that amino acid in the wt structures:
aa_counts_norm <- aa_counts %>%
  left_join(wt_fract) %>%
  mutate(norm_count = aa_count/n) %>%
  select(-c(aa_count, n))

# now I need to add up all aa within each (normaized!!!) bin to get bin totals and append this to the aa_counts
for_barplot <- aa_counts_norm %>%
  group_by(nat_bin) %>%
  mutate(bin_count = sum(norm_count)) %>%
  ungroup()

# adding the class labels
with_classes <- for_barplot %>%
  mutate(class = map_chr(aa_wt, calc_class))

for_barplot_3 <- with_classes %>%
  mutate(freq = norm_count/bin_count) %>%
  select(-c(norm_count, bin_count)) %>%
  mutate(aa_wt = fct_rev(fct_relevel(aa_wt, "G", "L", "P", "I", "W", "F", "C", "V", "D", "A", "M", "T", "S", "Y", "E", "N", "R", "H", "K", "Q"))) %>%
  mutate(class = fct_relevel(class, "aliphatic", "small_polar", "negative", "positive", "aromatic", "unique"))


# figuring out the order for the first facet:
order <- for_barplot_3 %>%
  filter(nat_bin == "Natural Frequency of 80-100%")

plot_g <- for_barplot_3 %>%
  ggplot(aes(x = freq, y = aa_wt, fill = class)) +
  geom_col(alpha = 0.75) +
  #facet_wrap(vars(fct_relevel(nat_bin, "Natural Frequency of 80-100%", "Natural Frequency of 0-20%")), ncol = 1) +
  scale_fill_manual(
    values = fills,
    labels = c("aliphatic", "small polar", "negative", "positive", "aromatic", "unique")) +
  scale_x_continuous(
    name = "Frequency",
    limits = c(0.0, 0.165),
    breaks = seq(0.0, 0.16, by = 0.04),
    expand = c(0, 0)) + 
  scale_y_discrete(
    name = "Wild type amino acid",
    expand = c(0.03, 0.03)) + 
  theme_cowplot(16) +
  theme(
    axis.text = element_text(color = "black", size = 14),
    strip.text.x = element_text(size = 16),
    panel.grid.major.x = element_line(color = "grey92", size=0.5),
    panel.grid.minor.x = element_line(color = "grey92", size=0.5),
    panel.spacing = unit(2, "lines"),
    legend.position = "none")

plot_g

ggsave(filename = "./analysis/figures/nat_aa_freq.png", plot = plot_g, width = 11, height = 8)

figure_final <- plot_grid(plot_e, plot_g, plot_train, nrow = 1, align = "h", labels = c('a', 'b', 'c'), rel_widths = c(1, 1, 1.5))

ggsave(filename = paste("./analysis/figures/aa_dist_CNN_conf.png"), plot = figure_final, width = 11, height = 9)

#===============================================================================================
#  consensus amino acid distributions within their MSA frequency bins (normalized by natural abundance)
#===============================================================================================

calc_class <- function(x) {
  aliphatic = c("M", "L", "I", "V", "A")
  small_polar = c("C", "S", "T", "N", "Q")
  negative = c("D", "E")
  positive = c("R", "K")
  aromatic = c("H", "Y", "W", "F")
  unique = c("P", "G")
  
  if (x %in% aliphatic) {
    return("aliphatic")
  }
  if (x %in% small_polar) {
    return("small_polar")
  }
  if (x %in% negative) {
    return("negative")
  }
  if (x %in% positive) {
    return("positive")
  }
  if (x %in% aromatic) {
    return("aromatic")
  }
  if (x %in% unique) {
    return("unique")
  }
  return("not found")
}

get_pred_bin <- function(x) {
  
  if (x > 0.80 & x <= 1.0) {
    return("Natural Frequency of 80-100%")
  }
  else{
    return(NA)
  }
}

# set working directory to: "Desktop/Natural_var_project/"
# loading data
cnn_data <- read.csv(file = "./data/PSICOV_box_20/output/cnn_wt_max_freq.csv", header=TRUE, sep=",")
natural_data <- read.csv(file = "./data/PSICOV_box_20/output/natural_max_freq_files/natural_max_freq_all.csv", header=TRUE, sep=",")

joined_data <- rbind(x = cnn_data, y = natural_data)

joined_data_trimmed <- joined_data %>%
  filter(!gene %in% c('1dbx', '1eaz', '1fvg', '1k7j', '1kq6', '1kw4', '1lpy', '1ne2', '1ny1', '1pko', '1rw1', '1vhu', '1w0h', '1wkc', '2tps'))

wide <- joined_data_trimmed %>%
  select(-c(aa_class, class_freq)) %>%
  pivot_wider(names_from = group, values_from = c(freq, aa)) 

wide_new <- wide %>%
  na.omit() %>%
  mutate(nat_bin = map_chr(freq_natural_max, get_pred_bin)) %>%
  na.omit()

new_data <- wide_new %>%
  select(aa_natural_max)

# finds the count of each aa acid per bin:
aa_counts <- new_data %>%
  count(aa_natural_max) %>%
  mutate(aa_count = n) %>%
  select(-n) %>%
  ungroup()

# getting aa fractions in the MSA sequences (# alanines/ # total sites)
natural_abundance <- read.csv(file = "./output/output_PSICOV/natural_abundance.csv", header=TRUE, sep=",")

natural_abundance <- natural_abundance %>%
  rename(aa_natural_max = aa) # renaming to incorrect name to join later on.
  

# normalizing the aa count per bin: 
# multiplying the aa count per bin by the total fraction of that amino acid in the wt structures:
aa_counts_norm <- aa_counts %>%
  left_join(natural_abundance) %>%
  mutate(norm_count = aa_count/count) %>%
  select(-c(aa_count, count))

# now I need to add up all aa within each (normaized!!!) bin to get bin totals and append this to the aa_counts
for_barplot <- aa_counts_norm %>%
  mutate(total_count = sum(norm_count)) %>%
  ungroup()

# adding the class labels
with_classes <- for_barplot %>%
  mutate(class = map_chr(aa_natural_max, calc_class))

for_barplot_4 <- with_classes %>%
  mutate(freq = norm_count/total_count) %>%
  select(-c(norm_count, total_count)) %>%
  mutate(aa_natural_max = fct_rev(fct_relevel(aa_natural_max, "G", "L", "P", "I", "W", "F", "C", "V", "D", "A", "M", "T", "S", "Y", "E", "N", "R", "H", "K", "Q"))) %>%
  mutate(class = fct_relevel(class, "aliphatic", "small_polar", "negative", "positive", "aromatic", "unique"))

plot_h <- for_barplot_4 %>%
  ggplot(aes(x = freq, y = aa_natural_max, fill = class)) +
  geom_col(alpha = 0.75) +
  scale_fill_manual(
    values = fills,
    labels = c("aliphatic", "small polar", "negative", "positive", "aromatic", "unique")) +
  scale_x_continuous(
    name = "Frequency",
    limits = c(0.0, 0.165),
    breaks = seq(0.0, 0.16, by = 0.04),
    expand = c(0, 0)) + 
  scale_y_discrete(
    name = "Natural amino acid",
    expand = c(0.03, 0.03)) + 
  theme_cowplot(16) +
  theme(
    axis.text = element_text(color = "black", size = 14),
    strip.text.x = element_text(size = 16),
    panel.grid.major.x = element_line(color = "grey92", size=0.5),
    panel.grid.minor.x = element_line(color = "grey92", size=0.5),
    panel.spacing = unit(2, "lines"),
    legend.position = "none")

plot_h

ggsave(filename = "./analysis/figures/cons_aa_dist.png", plot = plot_h, width = 11, height = 8)

figure_final <- plot_grid(plot_e, plot_h, plot_train, nrow = 1, align = "h", labels = c('a', 'b', 'c'), rel_widths = c(1, 1, 1.5))

ggsave(filename = paste("./analysis/figures/aa_dist_CNN_conf_natural.png"), plot = figure_final, width = 11, height = 9)

#=====================================================================================
# Natural N-eff distribution as a function of CNN confidence bins
#=====================================================================================


# set working directory to: "Desktop/Natural_var_project/"
# loading data
cnn_data <- read.csv(file = "./data/PSICOV_box_20/output/cnn_wt_max_freq.csv", header=TRUE, sep=",")

joined_data_trimmed <- cnn_data %>%
  filter(!gene %in% c('1dbx', '1eaz', '1fvg', '1k7j', '1kq6', '1kw4', '1lpy', '1ne2', '1ny1', '1pko', '1rw1', '1vhu', '1w0h', '1wkc', '2tps'))

wide <- joined_data_trimmed %>%
  select(-c(aa_class, class_freq)) %>%
  pivot_wider(names_from = group, values_from = c(freq, aa)) %>%
  select(c(gene, position, freq_predicted))


get_pred_bin <- function(x) {
  
  if (x > 0 & x <= 0.2) {
    return("(0-0.2]")
  }
  else if (x > 0.2 & x <= 0.4) {
    return("(0.2-0.4]")
  }
  else if (x > 0.4 & x <= 0.6) {
    return("(0.4-0.6]")
  }
  else if (x > 0.6 & x <= 0.8) {
    return("(0.6-0.8]")
  }
  else if (x > 0.8 & x <= 1.0) {
    return("(0.8-1.0]")
  }
}


wide_new <- wide %>%
  na.omit() %>%
  mutate(pred_bin = map_chr(freq_predicted, get_pred_bin))

conf_data <- wide_new %>%
  select(c(gene, position, pred_bin)) 
  
#loading the n-eff values:
natural_var <- read.csv(file="./output/output_PSICOV/stats_align_all.csv", header=TRUE, sep=",")

natural_var_clean <- natural_var %>%
  select(c(gene, position, n_eff)) 

joined <- inner_join(conf_data, natural_var_clean)


stat_data <- joined %>%
  select(-c(position, gene)) %>%
  group_by(pred_bin) %>%
  summarise(estimate = mean(n_eff),
            std_error = sd(n_eff)/sqrt(length(n_eff)))


plot_n_eff <- joined %>%
  ggplot(aes(y = n_eff, x = pred_bin)) +
  geom_violin(alpha = 0.6, size = 0.4, bw = 0.8, fill = "#99a88a", color = "#4d5841") +
  geom_pointrange(data = stat_data, aes(x = pred_bin,
                                        y = estimate,
                                        ymin = estimate - 1.96*std_error,
                                        ymax = estimate + 1.96*std_error),
                  color = "black", alpha = 0.7, size = 0.2) +
  #stat_summary(fun.data=data_summary, color = "black", alpha = 0.7) +
  theme_cowplot(16) + 
  theme(plot.title = element_text(hjust = 0, size = 16), 
        plot.subtitle = element_text(hjust = 0.5),
        panel.grid.major.y = element_line(color = "grey92", size=0.5),
        legend.position = "none") +
  scale_y_continuous(
    name = "Natural variation (n-eff)",
    limits = c(0, 20),
    breaks = seq(0, 20, by = 2),
    expand = c(0, 0)) +
  scale_x_discrete(
    name = "CNN confidence")

plot_n_eff


ggsave(filename = "./analysis/figures/n_eff_vs_pred.png", plot = plot_n_eff, width = 7, height = 4)

#===========================================================================================
# correlating n_eff natural vs. predicted probability for all positions and box sizes 
#=============================================================================================
# set working directory to: "Desktop/Natural_var_project/"
# loading data
cnn_data_12 <- read.csv(file = "./data/PSICOV_box_12/output/cnn_wt_max_freq.csv", header=TRUE, sep=",")
cnn_data_20 <- read.csv(file = "./data/PSICOV_box_20/output/cnn_wt_max_freq.csv", header=TRUE, sep=",")
cnn_data_30 <- read.csv(file = "./data/PSICOV_box_30/output/cnn_wt_max_freq.csv", header=TRUE, sep=",")
cnn_data_40 <- read.csv(file = "./data/PSICOV_box_40/output/cnn_wt_max_freq.csv", header=TRUE, sep=",")

#loading the natural n-eff values:
natural_var <- read.csv(file = "./output/output_PSICOV/stats_align_files/stats_align_60.csv", header=TRUE, sep=",")

natural_var <- natural_var %>%
  select(position, gene, n_eff) %>%
  mutate(group = "natural")

#box_size 12
cnn_var_12 <- cnn_var_12 %>%
  select(position, gene, n_eff, n_eff_class) %>%
  mutate(group = "predicted")

joined_12 <- rbind(natural_var, cnn_var_12) %>%
  mutate(box_size = "12") 









joined_data_trimmed <- cnn_data %>%
  filter(!gene %in% c('1dbx', '1eaz', '1fvg', '1k7j', '1kq6', '1kw4', '1lpy', '1ne2', '1ny1', '1pko', '1rw1', '1vhu', '1w0h', '1wkc', '2tps'))

wide <- joined_data_trimmed %>%
  select(-c(aa_class, class_freq)) %>%
  pivot_wider(names_from = group, values_from = c(freq, aa)) %>%
  select(c(gene, position, freq_predicted))


get_pred_bin <- function(x) {
  
  if (x > 0 & x <= 0.2) {
    return("(0-0.2]")
  }
  else if (x > 0.2 & x <= 0.4) {
    return("(0.2-0.4]")
  }
  else if (x > 0.4 & x <= 0.6) {
    return("(0.4-0.6]")
  }
  else if (x > 0.6 & x <= 0.8) {
    return("(0.6-0.8]")
  }
  else if (x > 0.8 & x <= 1.0) {
    return("(0.8-1.0]")
  }
}


wide_new <- wide %>%
  na.omit() %>%
  mutate(pred_bin = map_chr(freq_predicted, get_pred_bin))

conf_data <- wide_new %>%
  select(c(gene, position, pred_bin)) 

#loading the n-eff values:
natural_var <- read.csv(file = "./output/output_PSICOV/stats_align_files/stats_align_60.csv", header=TRUE, sep=",")

natural_var_clean <- natural_var %>%
  select(c(gene, position, n_eff)) 

joined <- inner_join(conf_data, natural_var_clean)


stat_data <- joined %>%
  select(-c(position, gene)) %>%
  group_by(pred_bin) %>%
  summarise(estimate = mean(n_eff),
            std_error = sd(n_eff)/sqrt(length(n_eff)))


plot_n_eff <- joined %>%
  ggplot(aes(y = n_eff, x = pred_bin)) +
  geom_violin(alpha = 0.6, size = 0.7, bw = 0.3, fill = "#99a88a", color = "#4d5841") +
  geom_pointrange(data = stat_data, aes(x = pred_bin,
                                        y = estimate,
                                        ymin = estimate - 1.96*std_error,
                                        ymax = estimate + 1.96*std_error),
                  color = "black", alpha = 0.7, size = 0.2) +
  #stat_summary(fun.data=data_summary, color = "black", alpha = 0.7) +
  theme_cowplot(16) + 
  theme(plot.title = element_text(hjust = 0, size = 16), 
        plot.subtitle = element_text(hjust = 0.5),
        panel.grid.major.y = element_line(color = "grey92", size=0.5),
        legend.position = "none") +
  scale_y_continuous(
    name = "Natural variation (n-eff)",
    limits = c(0, 20),
    breaks = seq(0, 20, by = 2),
    expand = c(0, 0)) +
  scale_x_discrete(
    name = "CNN confidence")

plot_n_eff


ggsave(filename = "./analysis/figures/n_eff_vs_pred.png", plot = plot_n_eff, width = 7, height = 4)

#=====================================================================================
# Consensus frequency as a function of CNN confidence bins
#=====================================================================================

# set working directory to: "Desktop/Natural_var_project/"
# loading data
cnn_data <- read.csv(file = "./data/PSICOV_box_20/output/cnn_wt_max_freq.csv", header=TRUE, sep=",")
natural_data <- read.csv(file = "./data/PSICOV_box_20/output/natural_max_freq_files/natural_max_freq_all.csv", header=TRUE, sep=",")

joined_data <- rbind(x = cnn_data, y = natural_data)

joined_data_trimmed <- joined_data %>%
  filter(!gene %in% c('1dbx', '1eaz', '1fvg', '1k7j', '1kq6', '1kw4', '1lpy', '1ne2', '1ny1', '1pko', '1rw1', '1vhu', '1w0h', '1wkc', '2tps'))

wide <- joined_data_trimmed %>%
  select(-c(aa_class, class_freq)) %>%
  pivot_wider(names_from = group, values_from = c(freq, aa)) 


get_pred_bin <- function(x) {
  
  if (x > 0 & x <= 0.2) {
    return("(0-0.2]")
  }
  else if (x > 0.2 & x <= 0.4) {
    return("(0.2-0.4]")
  }
  else if (x > 0.4 & x <= 0.6) {
    return("(0.4-0.6]")
  }
  else if (x > 0.6 & x <= 0.8) {
    return("(0.6-0.8]")
  }
  else if (x > 0.8 & x <= 1.0) {
    return("(0.8-1.0]")
  }
}


wide_new <- wide %>%
  na.omit() %>%
  mutate(pred_bin = map_chr(freq_predicted, get_pred_bin))


conf_data <- wide_new %>%
  select(c(gene, position, pred_bin, freq_natural_max)) 

data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}


plot_cons_freq <- conf_data %>%
  ggplot(aes(y = freq_natural_max, x = pred_bin)) +
  geom_violin(alpha = 0.6, size = 0.7, bw = 0.05, fill = "#d2a92d", color = "#8a7228") + 
  stat_summary(fun.data=data_summary, color = "black", alpha = 0.7) +
  theme_cowplot(16) + 
  theme(plot.title = element_text(hjust = 0, size = 16), 
        plot.subtitle = element_text(hjust = 0.5),
        panel.grid.major.y = element_line(color = "grey92", size=0.5),
        legend.position = "none") +
  scale_y_continuous(
    name = "Consensus Frequency",
    limits = c(0, 1.0),
    breaks = seq(0, 1.0, by = 0.2),
    expand = c(0, 0)) +
  scale_x_discrete(
    name = "CNN confidence \n (predicted probability)")

plot_cons_freq


ggsave(filename = "./analysis/figures/cons_freq_vs_pred.png", plot = plot_cons_freq, width = 7, height = 4)

#=====================================================================================
# Frequency of the predicted in MSA as a function of CNN confidence bins (mispredictions only)
#=====================================================================================

# set working directory to: "Desktop/Natural_var_project/"
# loading data
cnn_data <- read.csv(file = "./data/PSICOV_box_20/output/cnn_wt_max_freq.csv", header=TRUE, sep=",")
natural_data <- read.csv(file = "./data/PSICOV_box_20/output/natural_max_freq_files/natural_max_freq_all.csv", header=TRUE, sep=",")
natural_var <- read.csv(file="./output/output_PSICOV/stats_align_all.csv", header=TRUE, sep=",")


joined_data <- rbind(x = cnn_data, y = natural_data)

joined_data_trimmed <- joined_data %>%
  filter(!gene %in% c('1dbx', '1eaz', '1fvg', '1k7j', '1kq6', '1kw4', '1lpy', '1ne2', '1ny1', '1pko', '1rw1', '1vhu', '1w0h', '1wkc', '2tps'))

wide <- joined_data_trimmed %>%
  select(-c(aa_class, class_freq)) %>%
  pivot_wider(names_from = group, values_from = c(freq, aa)) 

joined_data2 <- inner_join(wide, natural_var)

joined_data_clean <- joined_data2 %>%
  select(c(gene, position, freq_predicted, aa_predicted, aa_wt, q_A:q_V))

mismatches <- joined_data_clean %>%
  filter(aa_predicted != aa_wt) %>%
  select(-aa_wt)

# renaming the natural frequencies
mismatches2 <- mismatches %>%
  pivot_longer(cols = c(q_A:q_V), names_to = "nat_aa", values_to = "nat_freq") %>%
  mutate(nat_aa = substr(nat_aa, 3, 3)) %>%
  filter(aa_predicted == nat_aa)
  
  
get_pred_bin <- function(x) {
  
  if (x > 0 & x <= 0.2) {
    return("(0-0.2]")
  }
  else if (x > 0.2 & x <= 0.4) {
    return("(0.2-0.4]")
  }
  else if (x > 0.4 & x <= 0.6) {
    return("(0.4-0.6]")
  }
  else if (x > 0.6 & x <= 0.8) {
    return("(0.6-0.8]")
  }
  else if (x > 0.8 & x <= 1.0) {
    return("(0.8-1.0]")
  }
}


for_plot <- mismatches2 %>%
  mutate(pred_bin = map_chr(freq_predicted, get_pred_bin)) %>%
  select(c(gene, position, nat_freq, pred_bin))



data_summary <- function(x) {
  m <- mean(x)
  ymin <- y-(sd(x)/sqrt(length(x)))
  ymax <- y+(sd(x)/sqrt(length(x)))
  return(c(y=m,ymin=ymin,ymax=ymax))
}

stat_data <- for_plot %>%
  select(-c(position, gene)) %>%
  group_by(pred_bin) %>%
  summarise(estimate = mean(nat_freq),
         std_error = sd(nat_freq)/sqrt(length(nat_freq)))
  


plot_nat_conf <- for_plot %>%
  ggplot(aes(y = nat_freq, x = pred_bin)) +
  geom_violin(alpha = 0.6, size = 0.4, bw = 0.02, fill = "#99a88a", color = "#4d5841") + 
  geom_pointrange(data = stat_data, aes(x = pred_bin,
                      y = estimate,
                      ymin = estimate - 1.96*std_error,
                      ymax = estimate + 1.96*std_error),
                  color = "black", alpha = 0.7, size = 0.4) +
  #stat_summary(fun.data=data_summary, color = "black", alpha = 0.7) +
  theme_cowplot(14) + 
  theme(plot.title = element_text(hjust = 0, size = 14), 
        plot.subtitle = element_text(hjust = 0.5),
        panel.grid.major.y = element_line(color = "grey92", size=0.5),
        panel.grid.minor.y = element_line(color = "grey92", size=0.5),
        legend.position = "none") +
  scale_y_continuous(
    name = "Amino acid frequency",
    limits = c(0, 1.0),
    breaks = seq(0, 1.0, by = 0.2),
    expand = c(0, 0)) +
  scale_x_discrete(
    name = "CNN confidence")

plot_nat_conf

#using standard_error
#ggsave(filename = "./analysis/figures/std_error.png", plot = plot_nat_conf, width = 8, height = 5)


ggsave(filename = "./analysis/figures/nat_freq_vs_pred.png", plot = plot_nat_conf, width = 8.5, height = 4.5)

#now I need a quick boxplot of the number of positions per group:

for_bar <- for_plot %>%
  group_by(pred_bin) %>%
  count()


bar_plot <- for_bar %>%
  ggplot(aes(x = factor(pred_bin), y = n)) +
  geom_col(fill = "#dbc7bd", color = "#cdaa98", alpha = 0.8) +
  scale_x_discrete(
    name = "Predicted probability") + 
  scale_y_continuous(
    name = "Number of positions per bin",
    limits = c(0, 3200),
    breaks = seq(0, 3000, by = 500),
    expand = c(0, 0)) +
  theme_cowplot(14)+
  theme(
    axis.text = element_text(color = "black", size = 14),
    strip.text.x = element_text(size = 16),
    panel.grid.major.y = element_line(color = "grey92", size=0.5),
    panel.grid.minor.y = element_line(color = "grey92", size=0.5)
  )

bar_plot

ggsave(filename = "./analysis/figures/positions_per_CNN_mispred.png", plot = bar_plot, width = 9, height = 5)

# now doing a correlation test between the two variables. 
cor<- mismatches2 %>%
  select(c(gene, freq_predicted, nat_freq)) %>%
  na.omit() %>%
  summarise(cor = cor(freq_predicted, nat_freq))

# fitting a linear model to data (getting R^2 and p-values)
lm_summary <- lm(mismatches2$nat_freq ~ mismatches2$freq_predicted)
summary(lm_summary)


#=======================================================================================
# Performing a paired t-test for mispredictions (between confidence bins)
#=======================================================================================


for_t_test <- for_plot %>%
  select(-c(gene, position)) %>%
  group_by(pred_bin) %>%
  mutate(row = row_number()) %>%
  ungroup() %>%
  pivot_wider(names_from = pred_bin, values_from = nat_freq) %>%
  select(-row)



results_1 <- t.test(for_t_test$`(0-0.2]`, for_t_test$`(0.2-0.4]`, paired = FALSE, alternative = "two.sided")
results_1
#p-value = 4.421e-06 (there is a significant difference)

results_2 <- t.test(for_t_test$`(0.2-0.4]`, for_t_test$`(0.4-0.6]`, paired = FALSE, alternative = "two.sided")
results_2
#p-value = 3.286e-12 (there is a significant difference)

results_3 <- t.test(for_t_test$`(0.4-0.6]`, for_t_test$`(0.6-0.8]`, paired = FALSE, alternative = "two.sided")
results_3
#p-value = 1.228e-05 (there is a significant difference)

results_4 <- t.test(for_t_test$`(0.6-0.8]`, for_t_test$`(0.8-1.0]`, paired = FALSE, alternative = "two.sided")
results_4
#p-value = 0.01287 (there is a significant difference)

results_5 <- t.test(for_t_test$`(0-0.2]`, for_t_test$`(0.8-1.0]`, paired = FALSE, alternative = "two.sided")
results_5
#p-value = p-value < 2.2e-16 (there is a significant difference)


    