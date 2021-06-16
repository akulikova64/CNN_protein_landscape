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


data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}


plot_a <- for_plot_a %>%
  ggplot(aes(y = freq_predict_cons, x = pred_bin)) +
  geom_violin(alpha = 0.6, size = 0.7, fill = "#d2a92d", color = "#8a7228") + 
  stat_summary(fun.data=data_summary, color = "black", alpha = 0.7) +
  theme_cowplot(16) + 
  theme(plot.title = element_text(hjust = 0, size = 16), 
        plot.subtitle = element_text(hjust = 0.5),
        panel.grid.major.y = element_line(color = "grey92", size=0.5),
        legend.position = "none") +
  scale_y_continuous(
    name = "Consensus Accuracy",
    limits = c(0.0, 1.0),
    breaks = seq(0, 1.0, by = 0.2),
    expand = c(0, 0)) +
  scale_x_discrete(
    name = "Predicted Probability")

plot_a

# Now, looking at wt predictions.

for_plot_b <- matches %>%
  select(c(gene, position, pred_bin, match_predict_wt )) %>%
  group_by(gene, pred_bin) %>%
  summarise(freq_predict_wt = sum(match_predict_wt, na.rm = TRUE)/sum(!is.na(match_predict_wt)))


plot_b <- for_plot_b %>%
  ggplot(aes(y = freq_predict_wt, x = pred_bin)) +
  geom_violin(alpha = 0.6, size = 0.7, fill = "#d2a92d", color = "#8a7228") + 
  stat_summary(fun.data=data_summary, color = "black", alpha = 0.7) +
  theme_cowplot(16) + 
  theme(plot.title = element_text(hjust = 0, size = 16), 
        plot.subtitle = element_text(hjust = 0.5),
        panel.grid.major.y = element_line(color = "grey92", size=0.5),
        legend.position = "none") +
  scale_y_continuous(
    name = "Wild Type Accuracy",
    limits = c(0.0, 1.0),
    breaks = seq(0, 1.0, by = 0.2),
    expand = c(0, 0)) +
  scale_x_discrete(
    name = "Predicted Probability")

plot_b


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
# Looking how natural/predicted amino acid distributions relate to CNN confidence 
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
  
  if (x > 0 & x <= 0.25) {
    return("Predicted at 0-25% confidence")
  }
  else if (x > 0.25 & x <= 0.50) {
    return("Predicted at 25-50% confidence")
  }
  else if (x > 0.50 & x <= 0.75) {
    return("Predicted at 50-75% confidence")
  }
  else if (x > 0.75 & x <= 1.0) {
    return("Predicted at 75-100% confidence")
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
  mutate(pred_bin = map_chr(freq_predicted, get_pred_bin))

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
  mutate(sum = sum(n)) %>%
  mutate(fract = n/sum) %>%
  select(-c(n, sum)) %>%
  rename(aa_predicted = aa_wt) #renaming to an incorrect name so that we can do a join later

# normalizing the aa count per bin: 
# multiplying the aa count per bin by the total fraction of that amino acid in the wt structures:
aa_counts_norm <- aa_counts %>%
  left_join( wt_fract) %>%
  mutate(norm_count = aa_count*fract) %>%
  select(-c(aa_count, fract))
  

# what is this?
aa_per_bin <- with_classes %>%
  group_by(pred_bin) %>%
  count() %>%
  mutate(bin_count = n) %>%
  select(-n) %>%
  ungroup()

for_barplot <- left_join(aa_counts, aa_per_bin)

with_classes <- for_barplot %>%
  mutate(class = map_chr(aa_predicted, calc_class))


for_barplot_2 <- with_classes %>%
  mutate(freq = aa_count/bin_count) %>%
  select(-c(aa_count, bin_count))

plot_e <- for_barplot_2 %>%
  ggplot(aes(x = freq, y = fct_reorder(aa_predicted, freq), fill = class)) +
  geom_col() +
  facet_wrap(vars(fct_relevel(pred_bin, "Predicted at 75-100% confidence", "Predicted at 50-75% confidence", "Predicted at 25-50% confidence", "Predicted at 0-25% confidence"))) +
  scale_x_continuous(
    name = "Amino acid frequency \n (Normalized by aa abundance)",
    limits = c(0.0, 0.2),
    breaks = seq(0.0, 0.2, by = 0.04),
    expand = c(0, 0)) + 
  scale_y_discrete(
    name = "Predicted amino acid",
    expand = c(0.03, 0.03)) + 
  theme_cowplot(14) +
  theme(
    axis.text = element_text(color = "black", size = 14),
    strip.text.x = element_text(size = 16),
    panel.grid.major.x = element_line(color = "grey92", size=0.5),
    panel.grid.minor.x = element_line(color = "grey92", size=0.5),
    panel.spacing = unit(2, "lines"))

plot_e

ggsave(filename = "./analysis/figures/Pred_aa_freq.png", plot = plot_e, width = 11, height = 8)








