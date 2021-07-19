library(tidyverse)
library(cowplot)
library(ggforce)

# misprediction analysis (what amino acids are more likely to be mispredicted)

# set working directory to: "Desktop/Natural_var_project/"
# loading data
cnn_data <- read.csv(file = "./data/PSICOV_box_20/output/cnn_wt_max_freq.csv", header=TRUE, sep=",")
natural_data <- read.csv(file = "./data/PSICOV_box_20/output/natural_max_freq_files/natural_max_freq_all.csv", header=TRUE, sep=",")

joined_data <- rbind(x = cnn_data, y = natural_data)

joined_data_trimmed <- joined_data %>%
  filter(!gene %in% c('1dbx', '1eaz', '1fvg', '1k7j', '1kq6', '1kw4', '1lpy', '1ne2', '1ny1', '1pko', '1rw1', '1vhu', '1w0h', '1wkc', '2tps'))

wide <- joined_data_trimmed %>%
  select(-c(aa_class, class_freq)) %>%
  pivot_wider(names_from = group, values_from = c(freq, aa)) %>%
  select(c(gene, position, aa_predicted, aa_wt))


mismatches <- wide %>%
  filter(aa_predicted != aa_wt) %>%
  select(aa_wt)

#===============================================================================================
#  mispredicted amino acid distributions
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


# finds the count of each aa acid per bin:
aa_counts <- mismatches %>%
  count(aa_wt) %>%
  mutate(aa_count = n) %>%
  select(-n) %>%
  ungroup()

# getting aa fractions in the wild type structures (# alanines/ # total sites)
# wt_fract <- wide %>%
#   select(gene, position, aa_wt) %>%
#   na.omit() %>%
#   group_by(aa_wt) %>%
#   count() %>%
#   ungroup() 


# normalizing the aa count per bin: 
# multiplying the aa count per bin by the total fraction of that amino acid in the wt structures:
# aa_counts_norm <- aa_counts %>%
#   left_join(wt_fract) %>%
#   mutate(norm_count = aa_count/n) %>%
#   select(-c(aa_count, n))

# now I need to add up all aa within each (normaized!!!) bin to get bin totals and append this to the aa_counts
for_barplot <- aa_counts %>%
  mutate(total_count = sum(aa_count)) %>%
  ungroup()

# adding the class labels
with_classes <- for_barplot %>%
  mutate(class = map_chr(aa_wt, calc_class))

for_barplot_1 <- with_classes %>%
  mutate(freq = aa_count/total_count) %>%
  select(-c(aa_count, total_count)) %>%
  mutate(aa_wt = fct_rev(fct_relevel(aa_wt, "K", "E", "S", "A", "Q", "R", "N", "T", "V", "D", "L", "I", "H", "M", "F", "Y", "P", "C", "W", "G"))) %>%
  mutate(class = fct_relevel(class, "aliphatic", "small_polar", "negative", "positive", "aromatic", "unique"))

fills <- c("#990008", "#0a2575", "#b35900", "#1a6600", "#5c0679", "#9e9e2e")

plot_a <- for_barplot_1 %>%
  ggplot(aes(x = freq, y = aa_wt, fill = class)) +
  geom_col(alpha = 0.75) +
  scale_fill_manual(
    values = fills,
    labels = c("aliphatic", "small polar", "negative", "positive", "aromatic", "unique")) +
  scale_x_continuous(
    name = "Misprediction frequency",
    limits = c(0.0, 0.11),
    breaks = seq(0.0, 0.10, by = 0.02),
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

plot_a

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

#filtering for correct predictions only:
# wide <- wide %>%
#   filter(aa_predicted == aa_wt)

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


for_barplot <- aa_counts %>%
  group_by(pred_bin) %>%
  mutate(bin_count = sum(aa_count)) %>%
  ungroup()

# adding the class labels
with_classes <- for_barplot %>%
  mutate(class = map_chr(aa_predicted, calc_class))

for_barplot_2 <- with_classes %>%
  mutate(freq = aa_count/bin_count) %>%
  select(-c(aa_count, bin_count)) %>%
  mutate(aa_predicted = fct_rev(fct_relevel(aa_predicted, "K", "E", "S", "A", "Q", "R", "N", "T", "V", "D", "L", "I", "H", "M", "F", "Y", "P", "C", "W", "G"))) %>%
  mutate(class = fct_relevel(class, "aliphatic", "small_polar", "negative", "positive", "aromatic", "unique"))


# figuring out the order for the first facet:
order <- for_barplot_2 %>%
  filter(pred_bin == "Predicted at 80-100% confidence")


plot_b <- for_barplot_2 %>%
  ggplot(aes(x = freq, y = aa_predicted, fill = class)) +
  geom_col(alpha = 0.75) +
  scale_fill_manual(
    values = fills,
    labels = c("aliphatic", "small polar", "negative", "positive", "aromatic", "unique")) +
  scale_x_continuous(
    name = "Confident \n prediction frequency",
    limits = c(0.0, 0.22),
    breaks = seq(0.0, 0.20, by = 0.05),
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
    panel.spacing = unit(2, "lines"))

plot_b


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



figure_final <- plot_grid(plot_a, plot_b, nrow = 1, align = "h", labels = c('a', 'b'), rel_widths = c(1, 1.5))

ggsave(filename = paste("./analysis/figures/mispredic_unnorm.png"), plot = figure_final, width = 9, height = 9)

#---------------------------------------------------------------------------------
# getting odds ratio of mispredictions over correct predictions. 
#---------------------------------------------------------------------------------

#data from the natural amino acid distributions:
mispred_data <- for_barplot_1 %>%
  select(c(aa_wt, class, freq)) %>%
  rename(freq_mispred = freq,
         aa = aa_wt)

# nat_data <- for_barplot_3 %>%
#   select(c(aa_natural_max, freq)) %>%
#   rename(aa = aa_natural_max,
#          freq_nat = freq)

#correct predictions:

wide <- joined_data_trimmed %>%
  select(-c(aa_class, class_freq)) %>%
  pivot_wider(names_from = group, values_from = c(freq, aa)) 

wide <- wide %>%
  filter(aa_predicted == aa_wt)

# finds the count of each aa acid per bin:
aa_counts <- wide %>%
  select(aa_predicted) %>%
  count(aa_predicted) %>%
  mutate(aa_count = n) %>%
  select(-n) %>%
  ungroup()

sums <- aa_counts %>%
  mutate(bin_count = sum(aa_count)) %>%
  ungroup()

# adding the class labels
with_classes <- sums %>%
  mutate(class = map_chr(aa_predicted, calc_class))

freqs <- with_classes %>%
  mutate(freq = aa_count/bin_count)

pred_data <- freqs %>%
  select(c(aa_predicted, freq)) %>%
  rename(freq_pred = freq,
         aa = aa_predicted)

#getting odds ratio
for_odds <- inner_join(pred_data, mispred_data)

for_odds_2 <- for_odds %>%
  mutate(ratio = freq_mispred/freq_pred) %>%
  mutate(aa = fct_reorder(aa, ratio)) %>%
  mutate(class = fct_relevel(class, "aliphatic", "small_polar", "negative", "positive", "aromatic", "unique"))

plot_odds <- for_odds_2 %>%
  ggplot(aes(x = ratio, y = aa, fill = class)) +
  geom_col(alpha = 0.75) +
  scale_fill_manual(
    values = fills,
    labels = c("aliphatic", "small polar", "negative", "positive", "aromatic", "unique")) +
  scale_x_log10(
    name = "Mispredicted frequency over \n correctly predicted frequency",
    #limits = c(0, 3.5),
    #breaks = seq(0, 3, by = 0.5),
    expand = c(0, 0)) + 
  scale_y_discrete(
    name = "Amino acid",
    expand = c(0.03, 0.03)) + 
  geom_hline(yintercept = 10.5, linetype = "dashed", color = "grey20", alpha = 0.8, size = 0.85) +
  theme_cowplot(16) +
  theme(
    axis.text = element_text(color = "black", size = 14),
    strip.text.x = element_text(size = 16),
    panel.grid.major.x = element_line(color = "grey92", size=0.5),
    panel.grid.minor.x = element_line(color = "grey92", size=0.5),
    panel.spacing = unit(2, "lines"))

plot_odds

ggsave(filename = paste("./analysis/figures/aa_odds_mispred.png"), plot = plot_odds, width = 6.5, height = 9)

polar <- c("C", "S", "T", "N", "Q", "D", "E", "R", "K", "Y", "H")
non_polar <- c("G", "P", "I", "A", "V", "L", "F", "M", "W")

for_odds_3 <- for_odds_2 %>%
  mutate(hydro = ifelse(aa %in% polar, "polar", "non-polar"))
  
plot_odds_3 <- for_odds_3 %>%
  ggplot(aes(x = ratio, y = aa, fill = hydro)) +
  geom_col(alpha = 0.75) +
  scale_fill_manual(
    values = fills) +
    #labels = c("aliphatic", "small polar", "negative", "positive", "aromatic", "unique")) +
  scale_x_log10(
    name = "Mispredicted frequency over \n correctly predicted frequency",
    #limits = c(0, 3.5),
    #breaks = seq(0, 3, by = 0.5),
    expand = c(0, 0)) + 
  scale_y_discrete(
    name = "Amino acid",
    expand = c(0.03, 0.03)) + 
  geom_hline(yintercept = 10.5, linetype = "dashed", color = "grey20", alpha = 0.8, size = 0.85) +
  theme_cowplot(16) +
  theme(
    axis.text = element_text(color = "black", size = 14),
    strip.text.x = element_text(size = 16),
    panel.grid.major.x = element_line(color = "grey92", size=0.5),
    panel.grid.minor.x = element_line(color = "grey92", size=0.5),
    panel.spacing = unit(2, "lines"))

plot_odds_3

ggsave(filename = paste("./analysis/figures/aa_odds_mispred_hydro.png"), plot = plot_odds_3, width = 6.5, height = 9)



