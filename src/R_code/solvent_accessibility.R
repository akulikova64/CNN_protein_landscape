# this program bins by solvent accessibility and calculates amino acid distributions for each bin

library(tidyverse)
library(cowplot)
library(broom)
library(ggplot2)
library(colorspace)

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
    return("0")
  }
  if (x > 0.0 & x <= 0.1) {
    return("(0-0.1]")
  }
  if (x > 0.1 & x <= 0.2) {
    return("(0.1-0.2]")
  }
   if (x > 0.2 & x <= 0.3) {
     return("(0.2-0.3]")
   }
  if (x > 0.3 & x <= 0.4) {
    return("(0.3-0.4]")
  }
  if (x > 0.4 & x <= 0.5) {
    return("(0.4-0.5]")
  }
  if (x > 0.5 & x <= 0.6) {
    return("(0.5-0.6]")
  }
  if (x > 0.6 & x <= 0.7) {
    return("(0.6-0.7]")
  }
  if (x > 0.7) {
    return("> 0.7")
  }
  else {
    return(NA)
  }
}


# set working directory to: "Desktop/Natural_var_project/"
# loading data
sa_data <- read.csv(file = "./output/output_PSICOV/SASA_scores.csv", header=TRUE, sep=",")

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
  scale_x_continuous(
    name = "RSA",
    limits = c(0.0, 1.67),
    breaks = seq(0.0, 1.6, by = 0.2),
    expand = c(0,0)) +
  scale_y_continuous(
    expand = c(0,0)) +
  geom_density(alpha = 0.6, size = 0.4, bw = 0.02, fill = "#99a88a", color = "#4d5841") +
  theme_cowplot()
dens_2
ggsave(filename = paste("./analysis/figures/RSA_density.png"), plot = dens_2, width = 10, height = 4.5)


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
  mutate(aa = fct_rev(fct_relevel(aa, "V", "A", "L", "I", "G", "F", "T", "S", "C", "M", "P", "Y", "N", "W", "D", "E", "H", "Q", "K", "R"))) %>%
  #mutate(aa = fct_rev(fct_relevel(aa, "L", "V", "A", "I", "G", "F", "T", "S", "Y", "D", "P", "M", "N", "C", "E", "H", "Q", "W", "R", "K"))) %>%
  #mutate(aa = fct_rev(fct_relevel(aa, "K", "E", "G", "R", "D", "P", "N", "Q", "A", "S", "T", "L", "V", "H", "M", "Y", "F", "I", "W", "C"))) %>%
  mutate(class = fct_relevel(class, "aliphatic", "small_polar", "negative", "positive", "aromatic", "unique"))


# figuring out the order for the first facet:
order <- for_barplot  %>%
  filter(SA_bin == "SA = 0")


plot_a <- for_barplot %>%
  ggplot(aes(x = freq, y = aa, fill = class)) +
  geom_col(alpha = 0.75) +
  facet_wrap(vars(fct_relevel(SA_bin, "0","(0-0.1]","(0.1-0.2]","(0.2-0.3]","(0.3-0.4]","(0.4-0.5]","(0.5-0.6]","(0.6-0.7]","> 0.7")), ncol = 10) +
  scale_fill_manual(
    values = fills,
    labels = c("aliphatic", "small polar", "negative", "positive", "aromatic", "unique")) +
  scale_x_continuous(
    name = "Frequency",
    limits = c(0.0, 0.20),
    breaks = seq(0.0, 0.20, by = 0.1),
    expand = c(0, 0)) + 
  scale_y_discrete(
    name = "Wild type amino acid",
    expand = c(0.03, 0.03)) + 
  theme_cowplot(16) +
  theme(
    axis.text = element_text(color = "black", size = 12),
    strip.text.x = element_text(size = 12),
    panel.grid.major.x = element_line(color = "grey92", size=0.5),
    panel.grid.minor.x = element_line(color = "grey92", size=0.5),
    panel.spacing = unit(2, "lines"),
    legend.position = "top")

plot_a

ggsave(filename = paste("./analysis/figures/SA_10_bins_aliph.png"), plot = plot_a, width = 15, height = 5)

#ggsave(filename = paste("./analysis/figures/SA_7_bins_charged.png"), plot = plot_a, width = 12, height = 6)

#supplementary boxplot with SA distributions per bin. 

plot_b <- binned %>%
  ggplot(aes(x = fct_relevel(SA_bin, "0","(0-0.1]","(0.1-0.2]","(0.2-0.3]","(0.3-0.4]","(0.4-0.5]","(0.5-0.6]","(0.6-0.7]","> 0.7"))) +
  geom_bar(fill = "#988981", color = "#70635c", alpha = 0.8) +
  scale_x_discrete(
    name = "Relative Solvent Accessibiity") +
  scale_y_continuous(
    name = "Count",
    limits = c(0, 6000),
    breaks = seq(0, 5000, by = 1000),
    expand = c(0, 0)) +
  geom_text(stat = 'count', aes(label = ..count..), vjust = - 0.5) +
  theme_cowplot(14)+
  theme(
    axis.text = element_text(color = "black", size = 14),
    strip.text.x = element_text(size = 14),
    panel.grid.major.y = element_line(color = "grey92", size=0.5),
    panel.grid.minor.y = element_line(color = "grey92", size=0.5)
  )
plot_b

ggsave(filename = paste("./analysis/figures/SA_bins_barplot2.png"), plot = plot_b, width = 15, height = 5)


#-------------------------------------------------------------------------------
# Two extreames: High and low solvent accessibility
#-------------------------------------------------------------------------------


get_SA_extremes <- function(x) {
  
  if (x == 0.0) {
    return("0")
  }
  if (x > 0.7) {
    return("> 0.7")
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
  mutate(aa = fct_rev(fct_relevel(aa, "V", "A", "L", "I", "G", "F", "T", "S", "C", "M", "P", "Y", "N", "W", "D", "E", "H", "Q", "K", "R"))) %>%
  #mutate(aa = fct_rev(fct_relevel(aa, "L", "V", "A", "I", "G", "F", "T", "S", "Y", "D", "P", "M", "N", "C", "E", "H", "Q", "W", "R", "K"))) %>%
  # order base on SA > 50
  #mutate(aa = fct_rev(fct_relevel(aa, "K", "E", "G", "R", "D", "P", "N", "Q", "A", "S", "T", "L", "V", "H", "M", "Y", "F", "I", "W", "C"))) %>%
  mutate(class = fct_relevel(class, "aliphatic", "small_polar", "negative", "positive", "aromatic", "unique"))


# figuring out the order for the first facet:
order <- for_barplot  %>%
  filter(SA_bin == "0")


plot_c <- for_barplot %>%
  ggplot(aes(x = freq, y = aa, fill = class)) +
  geom_col(alpha = 0.75) +
  facet_wrap(vars(fct_relevel(SA_bin, "0", "> 0.7")), ncol = 2) +
  scale_fill_manual(
    values = fills,
    labels = c("aliphatic", "small polar", "negative", "positive", "aromatic", "unique")) +
  scale_x_continuous(
    name = "Frequency",
    limits = c(0.0, 0.20),
    breaks = seq(0.0, 0.20, by = 0.05),
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

ggsave(filename = paste("./analysis/figures/SA_extreams2.png"), plot = plot_c, width = 9, height = 6)

#===============================================================================================
#  ODDS RATIO amino acid distributions within their MSA frequency bins
#===============================================================================================

#data from the natural amino acid distributions:
low_data <- for_barplot %>%
  filter(SA_bin == "0") %>%
  select(c(aa, class, freq)) %>%
  rename(freq_SA_zero = freq)

high_data <- for_barplot %>%
  filter(SA_bin == "> 0.7") %>%
  select(c(aa, class, freq)) %>%
  rename(freq_SA_high = freq)

for_odds <- inner_join(low_data, high_data)

for_odds_2 <- for_odds %>%
  mutate(ratio = freq_SA_zero/freq_SA_high) %>%
  mutate(aa = fct_reorder(aa, ratio)) %>%
  mutate(class = fct_relevel(class, "aliphatic", "small_polar", "negative", "positive", "aromatic", "unique"))

plot_odds <- for_odds_2 %>%
  ggplot(aes(x = ratio, y = aa, fill = class)) +
  geom_col(alpha = 0.75) +
  scale_fill_manual(
    values = fills,
    labels = c("aliphatic", "small polar", "negative", "positive", "aromatic", "unique")) +
  scale_x_log10(
    name = "Low SA (=0) over \n High SA (>0.7)",
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

ggsave(filename = paste("./analysis/figures/aa_dist_SA_odds2.png"), plot = plot_odds, width = 8, height = 9)

#-----------------------------------------------------------------------------------
# dot plot (and error) of accuracy as a function of SA bins (use the 10 from above)
#-----------------------------------------------------------------------------------

# loading in data:

sa_data <- read.csv(file = "./output/output_PSICOV/SASA_scores.csv", header=TRUE, sep=",")
cnn_data <- read.csv(file = "./data/PSICOV_box_20/output/cnn_wt_max_freq.csv", header=TRUE, sep=",")

sa_data_new <- sa_data %>%
  select(-SASA_abs_total) %>%
  mutate(group = "wt",
         position = as.numeric(position)) %>%
  rename(aa_wt_sa = aa)
  

cnn_data_new <- cnn_data %>%
  select(-c(aa_class, class_freq)) %>%
  pivot_wider(names_from = group, values_from = c(freq, aa)) %>%
  select(-freq_wt)


joined_data <- inner_join(x = cnn_data_new, y = sa_data_new)

#bin each posotion by rsa
binned3 <- joined_data %>%
  mutate(SA_bin = map_chr(SASA_rel_total, get_SA_bin)) %>%
  na.omit()

#get the frequency of correct matches:
matches <- binned3 %>%
  mutate(match_predict_wt = aa_predicted == aa_wt) 

for_dot_plot <- matches %>%
  group_by(SA_bin) %>%
  summarise(freq_predict_wt = sum(match_predict_wt, na.rm = TRUE)/sum(!is.na(match_predict_wt)))


dot_plot <- for_dot_plot %>%
  ggplot(aes(x = fct_relevel(SA_bin, "0","(0-0.1]", "(0.1-0.2]", "(0.2-0.3]","(0.3-0.4]","(0.4-0.5]","(0.5-0.6]","(0.6-0.7]","> 0.7"), y = freq_predict_wt)) +
  geom_point(size = 2, color = "#2873bd") +
  scale_x_discrete(
    name = "Relative Solvent Accessibiity",
    expand = c(0.05, 0.05)) + 
  scale_y_continuous(
    name = "Accuracy \n (predicting wt)",
    limits = c(0.2, 1.0),
    breaks = seq(0.2, 1.0, by = 0.2),
    expand = c(0, 0)) + 
  theme_cowplot(16) +
  theme(
    axis.text = element_text(color = "black", size = 14),
    strip.text.x = element_text(size = 16),
    panel.grid.major.x = element_line(color = "grey92", size=0.5),
    panel.grid.minor.x = element_line(color = "grey92", size=0.5),
    panel.grid.major.y = element_line(color = "grey92", size=0.5),
    panel.grid.minor.y = element_line(color = "grey92", size=0.5),
    panel.spacing = unit(2, "lines"))

dot_plot

ggsave(filename = paste("./analysis/figures/dot_plot.png"), plot = dot_plot, width = 10, height = 5)

#-------------------------------------------------------------------------------------
# let's correlate the pred. prob (CNN conf) by Solvent Accessibility (SA)
#-------------------------------------------------------------------------------------

scatter_plot <- joined_data %>%
  ggplot(aes(x = SASA_rel_total, y = freq_predicted)) +
  #geom_pointrange() +
  geom_hex(bins = 30) +
  scale_x_continuous(
    name = "Relative Solvent Accessibiity (Å^2)",
    #limits = c(0.0, 1.7),
    #breaks = seq(0.0, 1.0, by = 0.2),
    expand = c(0, 0)) +
  scale_y_continuous(
    name = "CNN confidence",
    limits = c(0.0, 1.0),
    breaks = seq(0.0, 1.0, by = 0.1),
    expand = c(0, 0)) + 
  scale_fill_binned_sequential(palette = "Teal", limits = c(0, 50)) +
  theme_cowplot(16) +
  theme(
    axis.text = element_text(color = "black", size = 14),
    strip.text.x = element_text(size = 16),
    panel.grid.major.x = element_line(color = "grey92", size=0.5),
    panel.grid.minor.x = element_line(color = "grey92", size=0.5),
    panel.grid.major.y = element_line(color = "grey92", size=0.5),
    panel.grid.minor.y = element_line(color = "grey92", size=0.5),
    panel.spacing = unit(2, "lines"))

scatter_plot

ggsave(filename = paste("./analysis/figures/CNN_conf_vs_SA.png"), plot = scatter_plot, width = 8, height = 6)


#-------------------------------------------------------------------------------------
# violin plots of the pred. prob (CNN conf) by Solvent Accessibility (SA)
#-------------------------------------------------------------------------------------
data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}

violin_plot <- binned3 %>%
  ggplot(aes(x = fct_relevel(SA_bin, "0","(0-0.025]","(0.025-0.1]","(0.1-0.2]","(0.2-0.5]","(0.3-0.4]","(0.4-0.5]","(0.5-0.6]","(0.6-0.7]","> 0.7"), y = freq_predicted)) +
  geom_violin(alpha = 0.6, size = 0.4, bw = 0.02, fill = "#99a88a", color = "#4d5841") + 
  stat_summary(fun.data=data_summary, color = "black", alpha = 0.7) +
  scale_x_discrete(
    name = "Relative Solvent Accessibiity",
    #limits = c(0.0, 1.7),
    #breaks = seq(0.0, 1.0, by = 0.2),
    expand = c(0, 0)) +
  scale_y_continuous(
    name = "CNN confidence",
    limits = c(0.0, 1.0),
    breaks = seq(0.0, 1.0, by = 0.2),
    expand = c(0, 0)) + 
  scale_fill_binned_sequential(palette = "Teal", limits = c(0, 50)) +
  theme_cowplot(16) +
  theme(
    axis.text = element_text(color = "black", size = 14),
    strip.text.x = element_text(size = 16),
    panel.grid.major.y = element_line(color = "grey92", size=0.5),
    panel.grid.minor.y = element_line(color = "grey92", size=0.5),
    panel.spacing = unit(2, "lines"))

violin_plot

ggsave(filename = paste("./analysis/figures/CNN_conf_vs_SA_bins.png"), plot = violin_plot, width = 12, height = 5)

#-------------------------------------------------------------------------------------
# aa distribution plots comparing the high density areas in 2D histogram.
#-------------------------------------------------------------------------------------

high_conf <- binned3 %>%
  filter(SASA_rel_total > 0.4 & SASA_rel_total < 0.8) %>%
  filter(freq_predicted > 0.9) %>%
  select(c(gene, position, aa_wt)) %>%
  mutate(group = "high conf.")

low_conf <- binned3 %>%
  filter(SASA_rel_total > 0.4 & SASA_rel_total < 0.8) %>%
  filter(freq_predicted < 0.4 & freq_predicted > 0.15) %>%
  select(c(gene, position, aa_wt)) %>%
  mutate(group = "low conf.")

dist_data <- rbind(low_conf, high_conf)

# finds the count of each aa acid per bin:
aa_counts <- dist_data %>%
  group_by(group) %>%
  count(aa_wt) %>%
  mutate(aa_count = n) %>%
  select(-n) %>%
  ungroup()

# now I need to add up all aa within each bin to get bin totals and append this to the aa_counts
with_sum <- aa_counts %>%
  group_by(group) %>%
  mutate(bin_count = sum(aa_count)) %>%
  ungroup()

# adding the class labels
with_classes <- with_sum %>%
  mutate(class = map_chr(aa_wt, calc_class))  

for_barplot <- with_classes %>%
  mutate(freq = aa_count/bin_count) %>%
  select(-c(aa_count, bin_count)) %>%
  # order base on SA = 0
  mutate(aa_wt = fct_rev(fct_relevel(aa_wt, "G", "P",  "D", "L", "E", "K", "A", "R", "V", "I", "N", "S", "T", "F", "C", "H", "Y", "Q", "W", "M"))) %>%
  #mutate(aa = fct_rev(fct_relevel(aa, "L", "V", "A", "I", "G", "F", "T", "S", "Y", "D", "P", "M", "N", "C", "E", "H", "Q", "W", "R", "K"))) %>%
  # order base on SA > 50
  #mutate(aa = fct_rev(fct_relevel(aa, "K", "E", "G", "R", "D", "P", "N", "Q", "A", "S", "T", "L", "V", "H", "M", "Y", "F", "I", "W", "C"))) %>%
  mutate(class = fct_relevel(class, "aliphatic", "small_polar", "negative", "positive", "aromatic", "unique"))


# figuring out the order for the first facet:
order <- for_barplot  %>%
  filter(group == "high conf.")


plot_high_low <- for_barplot %>%
  ggplot(aes(x = freq, y = aa_wt, fill = class)) +
  geom_col(alpha = 0.75) +
  facet_wrap(vars(fct_relevel(group, "high conf.", "low conf.")), ncol = 2) +
  scale_fill_manual(
    values = fills,
    labels = c("aliphatic", "small polar", "negative", "positive", "aromatic", "unique")) +
  scale_x_continuous(
    name = "Frequency",
    limits = c(0.0, 0.50),
    breaks = seq(0.0, 0.50, by = 0.10),
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

plot_high_low

ggsave(filename = paste("./analysis/figures/SA_vs_conf.png"), plot = plot_high_low, width = 9, height = 6)

#===============================================================================================
#  ODDS RATIO of amino acid distributions (high over low conf)
#===============================================================================================

#
low_data <- for_barplot %>%
  filter(group == "low conf.") %>%
  select(c(aa_wt, class, freq)) %>%
  rename(freq_low = freq)

high_data <- for_barplot %>%
  filter(group == "high conf.") %>%
  select(c(aa_wt, class, freq)) %>%
  rename(freq_high = freq)

for_odds <- inner_join(low_data, high_data)

for_odds_2 <- for_odds %>%
  mutate(ratio = freq_high/freq_low) %>%
  mutate(aa_wt = fct_reorder(aa_wt, ratio)) %>%
  mutate(class = fct_relevel(class, "aliphatic", "small_polar", "negative", "positive", "aromatic", "unique"))

plot_odds <- for_odds_2 %>%
  ggplot(aes(x = ratio, y = aa_wt, fill = class)) +
  geom_col(alpha = 0.75) +
  scale_fill_manual(
    values = fills,
    labels = c("aliphatic", "small polar", "negative", "positive", "aromatic", "unique")) +
  scale_x_log10(
    name = "High conf. (>0.9) over \n low conf. (>0.15, <0.4)",
    #limits = c(0, 3.5),
    #breaks = seq(0, 3, by = 0.5),
    expand = c(0, 0)) + 
  scale_y_discrete(
    name = "Wild type amino acid",
    expand = c(0.03, 0.03)) + 
  geom_hline(yintercept = 15.5, linetype = "dashed", color = "grey20", alpha = 0.8, size = 0.85) +
  theme_cowplot(16) +
  theme(
    axis.text = element_text(color = "black", size = 14),
    strip.text.x = element_text(size = 16),
    panel.grid.major.x = element_line(color = "grey92", size=0.5),
    panel.grid.minor.x = element_line(color = "grey92", size=0.5),
    panel.spacing = unit(2, "lines"))

plot_odds

ggsave(filename = paste("./analysis/figures/SA_vs_conf_odds.png"), plot = plot_odds, width = 8, height = 9)






