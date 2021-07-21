# this program bins by solvent accessibility and calculates amino acid distributions for each bin

library(tidyverse)
library(cowplot)
library(broom)

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

get_SA_bin <- function(x) {
  
  if (x == 0.0) {
    return("SA = 0")
  }
  if (x > 0.0 & x <= 5.0) {
    return("SA = (0-5]")
  }
  if (x > 5.0 & x <= 10.0) {
    return("SA = (5-10]")
  }
  if (x > 10.0 & x <= 15.0) {
    return("SA = (10-15]")
  }
  if (x > 15.0 & x <= 25.0) {
    return("SA = (15-25]")
  }
  if (x > 25.0 & x <= 50.0) {
    return("SA = (25-50]")
  }
  if (x > 50.0) {
    return("SA > 50")
  }
}


# set working directory to: "Desktop/Natural_var_project/"
# loading data
sa_data <- read.csv(file = "./output/output_PSICOV/solvent_accessibility.csv", header=TRUE, sep=",")

#cnn_data <- read.csv(file = "./data/PSICOV_box_20/output/cnn_wt_max_freq.csv", header=TRUE, sep=",")
#natural_data <- read.csv(file = "./data/PSICOV_box_20/output/natural_max_freq_files/natural_max_freq_all.csv", header=TRUE, sep=",")


# SASA abs_total distribution
dens_1 <- sa_data %>%
  ggplot(aes(x = SASA_abs_total)) +
  geom_density()
dens_1

count_1 <- sa_data %>%
  ggplot(aes(x = SASA_abs_total)) +
  geom_histogram()
count_1


# SASA rel_total distribution
dens_2 <- sa_data %>%
  ggplot(aes(x = SASA_rel_total)) +
  geom_density()
dens_2

count_2 <- sa_data %>%
  ggplot(aes(x = SASA_rel_total)) +
  geom_histogram()
count_2


# get an extra column with the SA bin
binned <- sa_data %>%
  mutate(SA_bin = map_chr(SASA_rel_total, get_SA_bin)) %>%
  na.omit()

#count positions per bin
counts <- binned %>%
  group_by(SA_bin) %>%
  count()

#------------------------------------------------------------------------------------
# finds the count of each aa acid per bin:
aa_counts <- binned %>%
  group_by(SA_bin) %>%
  count(aa) %>%
  mutate(aa_count = n) %>%
  select(-n) %>%
  ungroup()

# now I need to add up all aa within each bin to get bin totals and append this to the aa_counts
with_sum <- aa_counts %>%
  group_by(SA_bin) %>%
  mutate(bin_count = sum(aa_count)) %>%
  ungroup()

# adding the class labels
with_classes <- with_sum %>%
  mutate(class = map_chr(aa, calc_class))  

for_barplot <- with_classes %>%
  mutate(freq = aa_count/bin_count) %>%
  select(-c(aa_count, bin_count)) %>%
  #mutate(aa = fct_rev(fct_relevel(aa, "L", "V", "A", "I", "G", "F", "T", "S", "Y", "D", "P", "M", "N", "C", "E", "H", "Q", "W", "R", "K"))) %>%
  mutate(aa = fct_rev(fct_relevel(aa, "K", "E", "G", "R", "D", "P", "N", "Q", "A", "S", "T", "L", "V", "H", "M", "Y", "F", "I", "W", "C"))) %>%
  mutate(class = fct_relevel(class, "aliphatic", "small_polar", "negative", "positive", "aromatic", "unique"))


# figuring out the order for the first facet:
order <- for_barplot  %>%
  filter(SA_bin == "SA > 50")


plot_a <- for_barplot %>%
  ggplot(aes(x = freq, y = aa, fill = class)) +
  geom_col(alpha = 0.75) +
  facet_wrap(vars(fct_rev(fct_relevel(SA_bin, "SA = 0", "SA = (0-5]", "SA = (5-10]", "SA = (10-15]", "SA = (15-25]", "SA = (25-50]", "SA > 50"))), ncol = 7) +
  scale_fill_manual(
    values = fills,
    labels = c("aliphatic", "small polar", "negative", "positive", "aromatic", "unique")) +
  scale_x_continuous(
    name = "Frequency",
    limits = c(0.0, 0.16),
    breaks = seq(0.0, 0.16, by = 0.08),
    expand = c(0, 0)) + 
  scale_y_discrete(
    name = "Wild type amino acid",
    expand = c(0.03, 0.03)) + 
  theme_cowplot(16) +
  theme(
    axis.text = element_text(color = "black", size = 12),
    strip.text.x = element_text(size = 16),
    panel.grid.major.x = element_line(color = "grey92", size=0.5),
    panel.grid.minor.x = element_line(color = "grey92", size=0.5),
    panel.spacing = unit(2, "lines"),
    legend.position = "top")

plot_a

ggsave(filename = paste("./analysis/figures/SA_7_bins_charged.png"), plot = plot_a, width = 12, height = 6)

#supplementary boxplot with SA distributions per bin. 

plot_b <- binned %>%
  ggplot(aes(x = fct_relevel(SA_bin, "SA = 0", "SA = (0-5]", "SA = (5-10]", "SA = (10-15]", "SA = (15-25]", "SA = (25-50]", "SA > 50"))) +
  geom_bar(fill = "#988981", color = "#70635c", alpha = 0.8) +
  scale_x_discrete(
    name = "") +
  scale_y_continuous(
    name = "Count",
    limits = c(0, 9000),
    breaks = seq(0, 9000, by = 2000),
    expand = c(0, 0)) +
  geom_text(stat = 'count', aes(label = ..count..), vjust = - 0.5) +
  theme_cowplot(14)+
  theme(
    axis.text = element_text(color = "black", size = 14),
    strip.text.x = element_text(size = 16),
    panel.grid.major.y = element_line(color = "grey92", size=0.5),
    panel.grid.minor.y = element_line(color = "grey92", size=0.5)
  )
plot_b

ggsave(filename = paste("./analysis/figures/SA_bins_barplot.png"), plot = plot_b, width = 10, height = 6)


#-------------------------------------------------------------------------------
# Two extreames: High and low solvent accessibility
#-------------------------------------------------------------------------------


get_SA_extremes <- function(x) {
  
  if (x == 0.0) {
    return("SA = 0")
  }
  if (x > 25.0) {
    return("SA > 25")
  }
  else {
    return(NA)
  }
}

# get an extra column with the SA bin
binned2 <- sa_data %>%
  mutate(SA_bin = map_chr(SASA_rel_total, get_SA_extremes)) %>%
  na.omit()

#count positions per bin
counts2 <- binned2 %>%
  group_by(SA_bin) %>%
  count()

#------------------------------------------------------------------------------------
# finds the count of each aa acid per bin:
aa_counts <- binned2 %>%
  group_by(SA_bin) %>%
  count(aa) %>%
  mutate(aa_count = n) %>%
  select(-n) %>%
  ungroup()

# now I need to add up all aa within each bin to get bin totals and append this to the aa_counts
with_sum <- aa_counts %>%
  group_by(SA_bin) %>%
  mutate(bin_count = sum(aa_count)) %>%
  ungroup()

# adding the class labels
with_classes <- with_sum %>%
  mutate(class = map_chr(aa, calc_class))  

for_barplot <- with_classes %>%
  mutate(freq = aa_count/bin_count) %>%
  select(-c(aa_count, bin_count)) %>%
  # order base on SA = 0
  mutate(aa = fct_rev(fct_relevel(aa, "L", "V", "A", "I", "G", "F", "T", "S", "Y", "D", "P", "M", "N", "C", "E", "H", "Q", "W", "R", "K"))) %>%
  # order base on SA > 50
  #mutate(aa = fct_rev(fct_relevel(aa, "K", "E", "G", "R", "D", "P", "N", "Q", "A", "S", "T", "L", "V", "H", "M", "Y", "F", "I", "W", "C"))) %>%
  mutate(class = fct_relevel(class, "aliphatic", "small_polar", "negative", "positive", "aromatic", "unique"))


# figuring out the order for the first facet:
order <- for_barplot  %>%
  filter(SA_bin == "SA = 0")


plot_c <- for_barplot %>%
  ggplot(aes(x = freq, y = aa, fill = class)) +
  geom_col(alpha = 0.75) +
  facet_wrap(vars(fct_relevel(SA_bin, "SA = 0", "SA > 25")), ncol = 2) +
  scale_fill_manual(
    values = fills,
    labels = c("aliphatic", "small polar", "negative", "positive", "aromatic", "unique")) +
  scale_x_continuous(
    name = "Frequency",
    limits = c(0.0, 0.16),
    breaks = seq(0.0, 0.16, by = 0.04),
    expand = c(0, 0)) + 
  scale_y_discrete(
    name = "Wild type amino acid",
    expand = c(0.03, 0.03)) + 
  theme_cowplot(16) +
  theme(
    axis.text = element_text(color = "black", size = 12),
    strip.text.x = element_text(size = 16),
    panel.grid.major.x = element_line(color = "grey92", size=0.5),
    panel.grid.minor.x = element_line(color = "grey92", size=0.5),
    panel.spacing = unit(2, "lines"))

plot_c

ggsave(filename = paste("./analysis/figures/SA_extreams.png"), plot = plot_c, width = 9, height = 6)

#===============================================================================================
#  ODDS RATIO amino acid distributions within their MSA frequency bins
#===============================================================================================

#data from the natural amino acid distributions:
low_data <- for_barplot %>%
  filter(SA_bin == "SA = 0") %>%
  select(c(aa, class, freq)) %>%
  rename(freq_SA_zero = freq)

high_data <- for_barplot %>%
  filter(SA_bin == "SA > 25") %>%
  select(c(aa, class, freq)) %>%
  rename(freq_SA_high = freq)

for_odds <- inner_join(low_data, high_data)

for_odds_2 <- for_odds %>%
  mutate(ratio = freq_SA_high/freq_SA_zero) %>%
  mutate(aa = fct_reorder(aa, ratio)) %>%
  mutate(class = fct_relevel(class, "aliphatic", "small_polar", "negative", "positive", "aromatic", "unique"))

plot_odds <- for_odds_2 %>%
  ggplot(aes(x = ratio, y = aa, fill = class)) +
  geom_col(alpha = 0.75) +
  scale_fill_manual(
    values = fills,
    labels = c("aliphatic", "small polar", "negative", "positive", "aromatic", "unique")) +
  scale_x_log10(
    name = "High SA (>25) over \n Low SA ( = 0)",
    #limits = c(0, 3.5),
    #breaks = seq(0, 3, by = 0.5),
    expand = c(0, 0)) + 
  scale_y_discrete(
    name = "Amino acid",
    expand = c(0.03, 0.03)) + 
  geom_hline(yintercept = 11.5, linetype = "dashed", color = "grey20", alpha = 0.8, size = 0.85) +
  theme_cowplot(16) +
  theme(
    axis.text = element_text(color = "black", size = 14),
    strip.text.x = element_text(size = 16),
    panel.grid.major.x = element_line(color = "grey92", size=0.5),
    panel.grid.minor.x = element_line(color = "grey92", size=0.5),
    panel.spacing = unit(2, "lines"))

plot_odds

ggsave(filename = paste("./analysis/figures/aa_dist_odds.png"), plot = plot_odds, width = 6.5, height = 9)




