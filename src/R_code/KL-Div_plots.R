library(tidyverse)
library(cowplot)
library(ggforce)


#========================================================================================================
# here I will bin positions based on KL-divergence and find the aa distributions within each bin. 
#========================================================================================================

# function receives two distributions x and y
get_D_KL_2 <- function(y, x) {
  sum(ifelse(x == 0 | y == 0, 0, x*log(x/y)))
}

# function receives two distributions y (natural) and x (predicted)
get_D_KL <- function(y, x) {
  sum <- 0
  for (i in 1:20) { 
    y[i] = y[i] + 1e-20 # the KL-div statistic can't be calculated when y is as zero
    x[i] = x[i] + 1e-20
    if (x[i] != 0 & y[i] != 0) {
      sum = sum + x[i] * log(x[i]/y[i])
    }
    if (x[i] != 0 & y[i] == 0) {
      sum = sum + x[i] * log(x[i]/y[i])
    }
  }
  return(sum)
}

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

cnn_var <- read.csv(file = paste0("./data/PSICOV_box_20/output/stats_cnn.csv"), header=TRUE, sep=",")
natural_var <- read.csv(file="./output/output_PSICOV/stats_align_all.csv", header=TRUE, sep=",")

cnn_var2 <- cnn_var %>%
  nest(q_cnn = q_A:q_V) %>%
  select(-c(q_aliphatic:n_eff_class, entropy, n_eff))

natural_var2 <- natural_var %>%
  nest(q_nat = q_A:q_V) %>%
  select(-c(q_aliphatic:n_eff_class, entropy, n_eff)) 

joined_data <- inner_join(cnn_var2, natural_var2)

joined_data_trimmed <- joined_data %>%
  filter(!gene %in% c('1dbx', '1eaz', '1fvg', '1k7j', '1kq6', '1kw4', '1lpy', '1ne2', '1ny1', '1pko', '1rw1', '1vhu', '1w0h', '1wkc', '2tps'))

with_d_kl <- joined_data_trimmed %>%
  rowwise() %>%
  mutate(D_KL = get_D_KL(as.numeric(q_nat), as.numeric(q_cnn))) %>%
  select(-c(q_nat, q_cnn))

threshold = "500"

high_div <- with_d_kl %>%
  arrange(desc(D_KL)) %>%
  #filter(D_KL >= 7.660969) %>% # top 1000
  filter(D_KL >= 18.57776) %>% #top 500 
  mutate(group = paste(threshold,"predictions \n least similar to nature"))

low_div <- with_d_kl %>%
  arrange(D_KL) %>%
  #filter(D_KL <= 0.1774506 ) %>% # bottom 1000
  filter(D_KL <= 0.08089885 ) %>% # bottom 500
  mutate(group = paste(threshold,"predictions \n most similar to nature"))

final_df <- rbind(high_div, low_div)  

#Now lets also add the predicted amino acid to this dataframe. 

cnn_data <- read.csv(file = "./data/PSICOV_box_20/output/cnn_wt_max_freq.csv", header=TRUE, sep=",")

predicted <- cnn_data %>%
  select(c(gene, group, position, aa)) %>%
  pivot_wider(names_from = group, values_from = aa) 

final_df_2 <- inner_join(final_df, predicted) %>%
  select(-wt_aa)


#==========================================================================================
#now lets start calculating the counts/ frequencies:
#==========================================================================================

# finds the count of each aa acid per bin:
aa_counts <- final_df_2 %>%
  group_by(group) %>%
  count(predicted) %>%
  mutate(aa_count = n) %>%
  select(-n) %>%
  ungroup()

# check if some amino acids are missing from the group, if so add them with count 0.
#aa_list = c("G", "C", "P", "W", "F", "D", "L", "A", "I", "V", "Y", "E", "T", "N", "R", "S", "H", "K", "Q", "M")
#most similar does not have "M" when threshold is 500!
new_row <- c(paste(threshold,"predictions \n most similar to nature"), "M", 0)
aa_counts <- rbind(aa_counts, new_row)
aa_counts$aa_count <- as.numeric(aa_counts$aa_count)  


# getting aa fractions in the wild type structures (# alanines/ # total sites) to use for normalization 
wt_fract <- predicted %>%
  select(gene, position, wt) %>%
  na.omit() %>%
  group_by(wt) %>%
  count() %>%
  ungroup() %>%
  rename(predicted = wt) #renaming to an incorrect name so that we can do a join later

# normalizing the aa count per bin: 
# multiplying the aa count per bin by the total fraction of that amino acid in the wt structures:
aa_counts_norm <- aa_counts %>%
  left_join(wt_fract) %>%
  mutate(norm_count = aa_count/n) %>%
  select(-c(aa_count, n))


# now I need to add up all aa within each (normaized!!!) bin to get bin totals and append this to the aa_counts
for_barplot <- aa_counts_norm %>%
  group_by(group) %>%
  mutate(group_count = sum(norm_count)) %>%
  ungroup()

# adding the class labels
with_classes <- for_barplot %>%
  mutate(class = map_chr(predicted, calc_class))

# edit from here
for_barplot_2 <- with_classes %>%
  mutate(freq = norm_count/group_count) %>%
  select(-c(norm_count, group_count)) %>%
  #mutate(predicted = fct_rev(fct_relevel(predicted, "G", "C", "P", "W", "F", "D", "L", "A", "I", "V", "Y", "E", "T", "N", "R", "S", "H", "K", "Q", "M"))) %>%
  mutate(predicted = fct_rev(fct_relevel(predicted, "C", "G", "W", "P", "F", "L", "D", "A", "I", "T", "V", "Y", "H", "S", "N", "E", "Q", "R", "K", "M"))) %>%
  mutate(class = fct_relevel(class, "aliphatic", "small_polar", "negative", "positive", "aromatic", "unique"))

# figuring out the order for the first facet:
order <- for_barplot_2 %>%
  filter(group == paste(threshold,"predictions \n most similar to nature"))

plot_a <- for_barplot_2 %>%
  ggplot(aes(x = freq, y = predicted, fill = class)) +
  geom_col(alpha = 0.7) +
  facet_wrap(vars(fct_relevel(group, paste(threshold,"predictions \n most similar to nature"), paste(threshold,"predictions \n least similar to nature"))), nrow = 2) +
  scale_fill_manual(
    values = c("#99004d", "#001a66", "#994d00", "#1a6600", "#330066", "#9e9e2e")) +
  scale_x_continuous(
    name = "Frequency within divergence group \n (normalized by wt aa abundance)",
    limits = c(0.0, 0.32),
    breaks = seq(0.0, 0.30, by = 0.05),
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
    panel.spacing = unit(2, "lines"),
    legend.position = "none")

plot_a

ggsave(filename = "./analysis/figures/D_KL_aa_dist.png", plot = plot_a, width = 11, height = 8)  
  
 
#==========================================================================================
# plotting the wt as a function of divergence
#==========================================================================================


wt <- cnn_data %>%
  select(c(gene, group, position, aa)) %>%
  pivot_wider(names_from = group, values_from = aa) 

final_df_2 <- inner_join(final_df, wt) %>%
  select(-c(predicted, wt_aa))


# finds the count of each aa acid per bin:
aa_counts <- final_df_2 %>%
  group_by(group) %>%
  count(wt) %>%
  mutate(aa_count = n) %>%
  select(-n) %>%
  ungroup()

# check for missing aa's and add them with count 0
new_row <- c(paste(threshold,"predictions \n most similar to nature"), "K", 0)
aa_counts <- rbind(aa_counts, new_row)
new_row <- c(paste(threshold,"predictions \n most similar to nature"), "M", 0)
aa_counts <- rbind(aa_counts, new_row)
aa_counts$aa_count <- as.numeric(aa_counts$aa_count)  

# getting aa fractions in the wild type structures (# alanines/ # total sites) to use for normalization 
wt_fract <- wt %>%
  select(gene, position, wt) %>%
  na.omit() %>%
  group_by(wt) %>%
  count() %>%
  ungroup() 

# normalizing the aa count per bin: 
# multiplying the aa count per bin by the total fraction of that amino acid in the wt structures:
aa_counts_norm <- aa_counts %>%
  left_join(wt_fract) %>%
  mutate(norm_count = aa_count/n) %>%
  select(-c(aa_count, n))


# now I need to add up all aa within each (normaized!!!) bin to get bin totals and append this to the aa_counts
for_barplot <- aa_counts_norm %>%
  group_by(group) %>%
  mutate(group_count = sum(norm_count)) %>%
  ungroup()

# adding the class labels
with_classes <- for_barplot %>%
  mutate(class = map_chr(wt, calc_class))

# edit from here
for_barplot_3 <- with_classes %>%
  mutate(freq = norm_count/group_count) %>%
  select(-c(norm_count, group_count)) %>%
  #mutate(wt = fct_rev(fct_relevel(wt, "G", "C", "P", "W", "F", "D", "L", "A", "I", "V", "Y", "E", "T", "N", "R", "S", "H", "K", "Q", "M"))) %>%
  mutate(wt = fct_rev(fct_relevel(wt, "C", "G", "W", "P", "F", "L", "D", "A", "I", "T", "V", "Y", "H", "S", "N", "E", "Q", "R", "K", "M"))) %>%
  mutate(class = fct_relevel(class, "aliphatic", "small_polar", "negative", "positive", "aromatic", "unique"))


# figuring out the order for the first facet:
#order <- for_barplot_2 %>%
  #filter(group == "Predictions \n most similar to nature")

plot_b <- for_barplot_3 %>%
  ggplot(aes(x = freq, y = wt, fill = class)) +
  geom_col(alpha = 0.7) +
  facet_wrap(vars(fct_relevel(group, paste(threshold,"predictions \n most similar to nature"), paste(threshold,"predictions \n least similar to nature"))), nrow = 2) +
  scale_fill_manual(
    values = c("#99004d", "#001a66", "#994d00", "#1a6600", "#330066", "#9e9e2e"),
    labels = c("aliphatic", "small polar", "negative", "positive", "aromatic", "unique")) +
  scale_x_continuous(
    name = "Frequency within divergence group \n (normalized by wt aa abundance)",
    limits = c(0.0, 0.32),
    breaks = seq(0.0, 0.30, by = 0.05),
    expand = c(0, 0)) + 
  scale_y_discrete(
    name = "Wild type amino acid",
    expand = c(0.03, 0.03)) + 
  theme_cowplot(14) +
  theme(
    axis.text = element_text(color = "black", size = 14),
    strip.text.x = element_text(size = 16),
    panel.grid.major.x = element_line(color = "grey92", size=0.5),
    panel.grid.minor.x = element_line(color = "grey92", size=0.5),
    panel.spacing = unit(2, "lines"))

plot_b

ggsave(filename = "./analysis/figures/D_KL_wt_dist.png", plot = plot_b, width = 11, height = 8)  

  
figure_final <- plot_grid(plot_a, plot_b, nrow = 1, align = "h", labels = c('a', 'b'), rel_widths = c(1, 1.25))

ggsave(filename = paste("./analysis/figures/figure_KL_div_",threshold,".png"), plot = figure_final, width = 10.5, height = 9)


#========================================================================================================================
# comparing n-eff distributions (between n_eff_cnn and n_eff_nat)
#========================================================================================================================
get_neff_group <- function(x, y) { #x is pred, y is nat
  if (round((x - y), 1) == 0 & round(y, 1) != 1) {
    return("pred=nat \n nat not 1")
  }
  if (round((x - y), 1) == 0 & round(y, 1) == 1) {
    return("pred=nat \n nat = 1")
  }
  if (x > y) {
    return("pred>nat")
  }
  if (x < y) {
    return("pred<nat")
  }
}


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

cnn_var <- read.csv(file = paste0("./data/PSICOV_box_20/output/stats_cnn.csv"), header=TRUE, sep=",")
natural_var <- read.csv(file="./output/output_PSICOV/stats_align_all.csv", header=TRUE, sep=",")

cnn_var2 <- cnn_var %>%
  select(-c(q_A:q_V)) %>%
  select(-c(q_aliphatic:n_eff_class, entropy)) %>%
  rename(n_eff_cnn = n_eff)

natural_var2 <- natural_var %>%
  select(-c(q_A:q_V)) %>%
  select(-c(q_aliphatic:n_eff_class, entropy)) %>%
  rename(n_eff_nat = n_eff)

joined_data <- inner_join(cnn_var2, natural_var2)

joined_data_trimmed <- joined_data %>%
  filter(!gene %in% c('1dbx', '1eaz', '1fvg', '1k7j', '1kq6', '1kw4', '1lpy', '1ne2', '1ny1', '1pko', '1rw1', '1vhu', '1w0h', '1wkc', '2tps'))


#Now lets also add the predicted amino acid to this dataframe. 

cnn_data <- read.csv(file = "./data/PSICOV_box_20/output/cnn_wt_max_freq.csv", header=TRUE, sep=",")

predicted <- cnn_data %>%
  select(c(gene, group, position, aa)) %>%
  pivot_wider(names_from = group, values_from = aa) 

final <- inner_join(joined_data_trimmed, predicted) %>%
  select(-wt_aa) %>%
  mutate(diff = round((n_eff_cnn - n_eff_nat), 3))

final_2 <- final %>%
  rowwise() %>%
  mutate(n_eff_group = get_neff_group(as.numeric(n_eff_cnn), as.numeric(n_eff_nat))) %>%
  mutate(facet_group = ifelse(n_eff_group == "pred=nat \n nat not 1" | n_eff_group == "pred=nat \n nat = 1", "pred = nat", "pred ≠ nat"))

counts <- final_2 %>%
  group_by(n_eff_group) %>%
  summarise(count = n()) %>%
  ungroup() %>%
  mutate(sum = sum(count)) %>%
  mutate(freq = count/sum) %>%
  select(c(n_eff_group, freq)) %>%
  mutate(facet_group = ifelse(n_eff_group == "pred=nat \n nat not 1" | n_eff_group == "pred=nat \n nat = 1", "pred = nat", "pred ≠ nat"))


summary <- final_2 %>%
  group_by(n_eff_group) %>%
  summarise(count = n())
summary

plot_a <- counts %>%
  filter(n_eff_group == "pred<nat" | n_eff_group == "pred>nat") %>%
  ggplot(aes(x=n_eff_group, y = freq)) +
  geom_col() 
plot_a

plot_b <- counts %>%
  filter(n_eff_group == "pred=nat \n nat = 1" | n_eff_group == "pred=nat \n nat not 1") %>%
  ggplot(aes(x=n_eff_group, y = freq)) +
  geom_col() 
plot_b

plot_c <- final_2 %>%
  ggplot(aes(x = fct_relevel(n_eff_group,"pred<nat", "pred>nat", "pred=nat \n nat = 1", "pred=nat \n nat not 1"))) +
  geom_bar(aes(y = ..prop.., group = 1), fill = "#988981", color = "#70635c", alpha = 0.8) +
  scale_x_discrete(
    name = "Variation group (n-eff)",
    labels = c("predicted < natural",
               "predicted > natural",
               "predicted = natural \n (n-eff = 1)",
               "predicted = natural \n (n-eff ≠ 1)")) + 
  scale_y_continuous(
    name = "Proportion",
    limits = c(0, 0.8),
    breaks = seq(0, 0.8, by = 0.1),
    expand = c(0, 0)) +
  theme_cowplot(14)+
  theme(
    axis.text = element_text(color = "black", size = 14),
    strip.text.x = element_text(size = 16),
    panel.grid.major.y = element_line(color = "grey92", size=0.5),
    panel.grid.minor.y = element_line(color = "grey92", size=0.5)
  )
plot_c

ggsave(filename = paste("./analysis/figures/n_eff_bins.png"), plot = plot_c, width = 8.5, height = 7)


#===============================================================================
# Next, we want to add in KL-divergence (make violin plots)
#===============================================================================

# function receives two distributions x and y
get_D_KL_2 <- function(y, x) {
  sum(ifelse(x == 0 | y == 0, 0, x*log(x/y)))
}

# function receives two distributions y (natural) and x (predicted)
get_D_KL <- function(y, x) {
  sum <- 0
  for (i in 1:20) { 
    y[i] = y[i] + 1e-20 # the KL-div statistic can't be calculated when y is zero
    x[i] = x[i] + 1e-20
    if (x[i] != 0 & y[i] != 0) {
      sum = sum + x[i] * log(x[i]/y[i])
    }
    if (x[i] != 0 & y[i] == 0) {
      sum = sum + x[i] * log(x[i]/y[i])
    }
  }
  return(sum)
}

get_neff_group <- function(x, y) { #x is pred, y is nat
  if (round((x - y), 1) == 0 & round(y, 1) != 1) {
    return("pred=nat \n nat not 1")
  }
  if (round((x - y), 1) == 0 & round(y, 1) == 1) {
    return("pred=nat \n nat = 1")
  }
  if (x > y) {
    return("pred>nat")
  }
  if (x < y) {
    return("pred<nat")
  }
}


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

cnn_var <- read.csv(file = paste0("./data/PSICOV_box_20/output/stats_cnn.csv"), header=TRUE, sep=",")
natural_var <- read.csv(file="./output/output_PSICOV/stats_align_all.csv", header=TRUE, sep=",")

cnn_var2 <- cnn_var %>%
  nest(q_cnn = q_A:q_V) %>%
  select(-c(q_aliphatic:n_eff_class, entropy)) %>%
  rename(n_eff_cnn = n_eff)

natural_var2 <- natural_var %>%
  nest(q_nat = q_A:q_V) %>%
  select(-c(q_aliphatic:n_eff_class, entropy)) %>%
  rename(n_eff_nat = n_eff)

joined_data <- inner_join(cnn_var2, natural_var2)

joined_data_trimmed <- joined_data %>%
  filter(!gene %in% c('1dbx', '1eaz', '1fvg', '1k7j', '1kq6', '1kw4', '1lpy', '1ne2', '1ny1', '1pko', '1rw1', '1vhu', '1w0h', '1wkc', '2tps'))

with_d_kl <- joined_data_trimmed %>%
  rowwise() %>%
  mutate(D_KL = get_D_KL(as.numeric(q_nat), as.numeric(q_cnn))) %>%
  select(-c(q_nat, q_cnn))

#Now lets also add the predicted amino acid to this dataframe. 

cnn_data <- read.csv(file = "./data/PSICOV_box_20/output/cnn_wt_max_freq.csv", header=TRUE, sep=",")

predicted <- cnn_data %>%
  select(c(gene, group, position, aa)) %>%
  pivot_wider(names_from = group, values_from = aa) 

final <- inner_join(with_d_kl, predicted) %>%
  select(-wt_aa) %>%
  mutate(diff = round((n_eff_cnn - n_eff_nat), 3))

final_2 <- final %>%
  rowwise() %>%
  mutate(n_eff_group = get_neff_group(as.numeric(n_eff_cnn), as.numeric(n_eff_nat))) %>%
  mutate(facet_group = ifelse(n_eff_group == "pred=nat \n nat not 1" | n_eff_group == "pred=nat \n nat = 1", "pred = nat", "pred ≠ nat"))

# counts <- final_2 %>%
#   group_by(n_eff_group) %>%
#   summarise(count = n()) %>%
#   ungroup() %>%
#   mutate(sum = sum(count)) %>%
#   mutate(freq = count/sum) %>%
#   select(c(n_eff_group, freq)) %>%
#   mutate(facet_group = ifelse(n_eff_group == "pred=nat \n nat not 1" | n_eff_group == "pred=nat \n nat = 1", "pred = nat", "pred ≠ nat"))


# summary <- final_2 %>%
#   group_by(n_eff_group) %>%
#   summarise(count = n())
# summary
# 
# plot_a <- counts %>%
#   filter(n_eff_group == "pred<nat" | n_eff_group == "pred>nat") %>%
#   ggplot(aes(x=n_eff_group, y = freq)) +
#   geom_col() 
# plot_a
# 
# plot_b <- counts %>%
#   filter(n_eff_group == "pred=nat \n nat = 1" | n_eff_group == "pred=nat \n nat not 1") %>%
#   ggplot(aes(x=n_eff_group, y = freq)) +
#   geom_col() 
# plot_b

plot_d <- final_2 %>%
  ggplot(aes(y = D_KL, x = fct_relevel(n_eff_group,"pred<nat", "pred>nat", "pred=nat \n nat = 1", "pred=nat \n nat not 1"))) +
  geom_sina(fill = "#988981", color = "#70635c", alpha = 0.7, size = 0.5, bw = 0.3) +
  scale_x_discrete(
    name = "Variation group (n-eff)",
    labels = c("predicted < natural",
               "predicted > natural",
               "predicted = natural \n (n-eff = 1)",
               "predicted = natural \n (n-eff ≠ 1)")) + 
  scale_y_continuous(
    name = "KL-Divergence",
    limits = c(0, 10),
    breaks = seq(0, 10, by = 1),
    expand = c(0, 0)) +
  theme_cowplot(14)+
  theme(
    axis.text = element_text(color = "black", size = 14),
    strip.text.x = element_text(size = 16),
    panel.grid.major.y = element_line(color = "grey92", size=0.5),
    panel.grid.minor.y = element_line(color = "grey92", size=0.5)
  )
plot_d

ggsave(filename = paste("./analysis/figures/n_eff_vs_KL.png"), plot = plot_d, width = 8.5, height = 7)

#Let's look at the fourth group: "pred=nat \n nat not 1"
for_observation <- final_2 %>%
  filter(n_eff_group == "pred=nat \n nat not 1") %>%
  select(c(gene, position, wt, predicted, D_KL, n_eff_group, n_eff_cnn, n_eff_nat))

#===============================================================================
# Let's get 3 KL groups
#===============================================================================

get_div_group <- function(n_eff_cnn, n_eff_nat, KL_div) { 
  if (round(n_eff_cnn, 1) == 1.0 & round(n_eff_nat, 1) == 1.0 & round(KL_div, 1) == 0) {
    return("n-eff = 1 \n KL-Div ~ 0")
  }
  if (round(n_eff_cnn, 1) > 1.5 & round(n_eff_nat, 1) > 1.5 & round(KL_div, 1) < 0.5) {
    return("n-eff > 1.5 \n KL-Div < 0.5")
  }
  if (round(KL_div, 1) > 5) {
    return("KL-Div > 5")
  }
  return(NA)
}

final_2

new_data <- final_2 %>%
  rowwise() %>%
  mutate(div_group = get_div_group(as.numeric(n_eff_cnn), as.numeric(n_eff_nat), as.numeric(D_KL))) %>%
  na.omit()

count_div_group <- new_data %>%
  group_by(div_group) %>%
  count()
count_div_group


#-------------------------------------------------------------------------------
#ok, now lets get the aa distributions for each group. 

# finds the count of each aa acid per bin:
aa_counts <- new_data %>%
  select(c(gene, position, wt, div_group)) %>%
  group_by(div_group) %>%
  count(wt) %>%
  mutate(aa_count = n) %>%
  select(-n) %>%
  ungroup()

aa <- c('A', 'R', 'N',  'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V')

all_aa <- tibble(wt = rep(aa, times = 3), 
                 count = rep(0, times = 60),
                 div_group = c(rep("KL-Div > 5", times = 20), rep("n-eff = 1 \n KL-Div ~ 0", times = 20), rep("n-eff > 1.5 \n KL-Div < 0.5", 20)))

with_missing_aa <- aa_counts %>%
  right_join(all_aa) %>%
  mutate(aa_count = ifelse(is.na(aa_count), 0, aa_count)) %>%
  select(-count)


# now I need to add up all aa within each (normaized!!!) bin to get bin totals and append this to the aa_counts
for_barplot <- with_missing_aa %>%
  group_by(div_group) %>%
  mutate(bin_count = sum(aa_count)) %>%
  ungroup()

# adding the class labels
with_classes <- for_barplot %>%
  mutate(class = map_chr(wt, calc_class))

for_barplot_2 <- with_classes %>%
  mutate(freq = aa_count/bin_count) %>%
  select(-c(aa_count, bin_count)) %>%
  mutate(wt = fct_rev(fct_relevel(wt, "E", "K", "A", "D", "V", "S", "T", "R", "I", "Q", "P", "N", "Y", "L", "F", "H", "M", "W", "C", "G"))) %>%
  mutate(class = fct_relevel(class, "aliphatic", "small_polar", "negative", "positive", "aromatic", "unique"))

# figuring out the order for the first facet:
order <- for_barplot_2 %>%
  filter(div_group == "n-eff > 1.5 \n KL-Div < 0.5")

fills <- c("#990008", "#0a2575", "#b35900", "#1a6600", "#5c0679", "#9e9e2e")

aa_dist_plot <- for_barplot_2 %>%
  ggplot(aes(x = freq, y = wt, fill = class)) +
  geom_col(alpha = 0.75) +
  facet_wrap(vars(fct_relevel(div_group, "n-eff > 1.5 \n KL-Div < 0.5", "KL-Div > 5", "n-eff = 1 \n KL-Div ~ 0")), ncol = 3) +
  scale_fill_manual(
    values = fills,
    labels = c("aliphatic", "small polar", "negative", "positive", "aromatic", "unique")) +
  scale_x_continuous(
    name = "Frequency",
    limits = c(0.0, 0.65),
    breaks = seq(0.0, 0.60, by = 0.1),
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
    panel.spacing = unit(2, "lines"))

aa_dist_plot

ggsave(filename = "./analysis/figures/aa_dist_kl_bins.png", plot = aa_dist_plot, width = 10, height = 6)

aa_dist_plot_2 <- for_barplot_2 %>%
  filter(div_group != "n-eff = 1 \n KL-Div ~ 0") %>%
  ggplot(aes(x = freq, y = wt, fill = class)) +
  geom_col(alpha = 0.75) +
  facet_wrap(vars(fct_relevel(div_group, "n-eff > 1.5 \n KL-Div < 0.5", "KL-Div > 5")), ncol = 2) +
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
    axis.text = element_text(color = "black", size = 14),
    strip.text.x = element_text(size = 16),
    panel.grid.major.x = element_line(color = "grey92", size=0.5),
    panel.grid.minor.x = element_line(color = "grey92", size=0.5),
    panel.spacing = unit(2, "lines"))

aa_dist_plot_2

ggsave(filename = "./analysis/figures/aa_dist_kl_bins_2.png", plot = aa_dist_plot_2, width = 10, height = 6)

#---------------------------------------------------------------------------------------------------------------------
# ok, I now need to make these same bar charts, only with secondary structure, instead of amino acids.

sec_struc <- read.csv(file = "./output/output_PSICOV/second_struc.csv", header=TRUE, sep=",")

sec_struc_counts <- sec_struc %>%
  group_by(second_struc) %>%
  count()

with_struc <- new_data %>%
  inner_join(sec_struc) %>%
  mutate(second_struc = ifelse(second_struc == "Coil", "Random coil", second_struc),
         second_struc = ifelse(second_struc == "Strand", "Beta strand", second_struc),
         second_struc = ifelse(second_struc == "AlphaHelix", "Alpha helix", second_struc),
         second_struc = ifelse(second_struc == "310Helix" | second_struc == "PiHelix", "Noncanon. helix", second_struc))
  

# finds the count of each structure per bin:
struc_counts <- with_struc %>%
  select(c(gene, position, second_struc, div_group)) %>%
  group_by(div_group) %>%
  count(second_struc) %>%
  mutate(struc_count = n) %>%
  select(-n) %>%
  ungroup()


# now I need to add up all aa within each (normaized!!!) bin to get bin totals and append this to the aa_counts
for_barplot <- struc_counts %>%
  group_by(div_group) %>%
  mutate(bin_count = sum(struc_count)) %>%
  ungroup()

# quick barplot with counts:
counts_table <- tibble(bin = c("n-eff > 1.5 \n KL-Div < 0.5", "KL-Div > 5", "n-eff = 1 \n KL-Div ~ 0"), 
                       count = c(1471, 1393, 136))
counts_plot <- counts_table %>%
  ggplot(aes(x = fct_relevel(bin, "n-eff > 1.5 \n KL-Div < 0.5", "KL-Div > 5", "n-eff = 1 \n KL-Div ~ 0"), y = count)) +
  geom_col(fill = "#988981", color = "#70635c", alpha = 0.8) +
  geom_text(aes(label = count, vjust = -0.25)) +
  scale_x_discrete(
    name = "") + 
  scale_y_continuous(
    name = "Count",
    limits = c(0, 1550),
    breaks = seq(0, 1400, by = 200),
    expand = c(0, 0)) +
  theme_cowplot(14)+
  theme(
    axis.text = element_text(color = "black", size = 14),
    strip.text.x = element_text(size = 16),
    panel.grid.major.y = element_line(color = "grey92", size=0.5),
    panel.grid.minor.y = element_line(color = "grey92", size=0.5)
  )
  
counts_plot

ggsave(filename = "./analysis/figures/kl_bin_counts.png", plot = counts_plot, width = 5.5, height = 4)


for_barplot_2 <- for_barplot %>%
  mutate(freq = struc_count/bin_count) %>%
  select(-c(struc_count, bin_count)) %>%
  mutate(second_struc = fct_rev(fct_relevel(second_struc, "Alpha helix", "Beta strand", "Turn", "Random coil", "Noncanon. helix", "Bridge")))
  
# figuring out the order for the first facet:
order <- for_barplot_2 %>%
  filter(div_group == "n-eff > 1.5 \n KL-Div < 0.5")

struc_fills = c("#bf8040", "#ac5396", "#70adc2", "#748f3d", "#cc5933", "#7070c2")

struc_dist_plot <- for_barplot_2 %>%
  ggplot(aes(x = freq, y = second_struc, fill = second_struc)) +
  geom_col(alpha = 0.9) +
  facet_wrap(vars(fct_relevel(div_group, "n-eff > 1.5 \n KL-Div < 0.5", "KL-Div > 5", "n-eff = 1 \n KL-Div ~ 0")), ncol = 3) +
  scale_fill_manual(
    values = struc_fills) +
  scale_x_continuous(
    name = "Frequency",
    limits = c(0.0, 0.45),
    breaks = seq(0.0, 0.40, by = 0.1),
    expand = c(0, 0)) + 
  scale_y_discrete(
    name = "Secondary structure",
    expand = c(0.03, 0.03)) + 
  theme_cowplot(16) +
  theme(
    axis.text = element_text(color = "black", size = 14),
    strip.text.x = element_text(size = 16),
    panel.grid.major.x = element_line(color = "grey92", size=0.5),
    panel.grid.minor.x = element_line(color = "grey92", size=0.5),
    panel.spacing = unit(2, "lines"),
    legend.position = "none")

struc_dist_plot

ggsave(filename = "./analysis/figures/struc_dist_kl_bins.png", plot = struc_dist_plot, width = 10, height = 5)

#====================================================================================
# Let's now look at the distribution of mispredicted sites across secondary structure
#====================================================================================

#Now lets also add the predicted amino acid to this dataframe. 

cnn_data <- read.csv(file = "./data/PSICOV_box_20/output/cnn_wt_max_freq.csv", header=TRUE, sep=",")

predicted <- cnn_data %>%
  select(c(gene, group, position, aa)) %>%
  pivot_wider(names_from = group, values_from = aa)

#get the wt freq:

predicted_2 <- cnn_data %>%
  select(c(gene, group, position, freq)) %>%
  pivot_wider(names_from = group, values_from = freq) %>%
  select(gene, position, wt) %>%
  rename(pred_wt = wt)

predicted_final <- inner_join(predicted, predicted_2)

# run code at top of script to get with_d_kls
final <- inner_join(with_d_kl, predicted_final) %>%
  select(-wt_aa) %>%
  mutate(diff = round((n_eff_cnn - n_eff_nat), 3))

final_2 <- final %>%
  rowwise() %>%
  mutate(n_eff_group = get_neff_group(as.numeric(n_eff_cnn), as.numeric(n_eff_nat))) %>%
  mutate(facet_group = ifelse(n_eff_group == "pred=nat \n nat not 1" | n_eff_group == "pred=nat \n nat = 1", "pred = nat", "pred ≠ nat"))


mispred <- final_2 %>%
  filter(predicted != wt & pred_wt < 0.1)


with_struc1 <- mispred %>%
  inner_join(sec_struc) %>%
  mutate(second_struc = ifelse(second_struc == "Coil", "Random coil", second_struc),
         second_struc = ifelse(second_struc == "Strand", "Beta strand", second_struc),
         second_struc = ifelse(second_struc == "AlphaHelix", "Alpha helix", second_struc),
         second_struc = ifelse(second_struc == "310Helix" | second_struc == "PiHelix", "Noncanon. helix", second_struc))

# finds the count of each aa acid per bin:
struc_counts1 <- with_struc1 %>%
  select(c(gene, position, second_struc)) %>%
  count(second_struc) %>%
  mutate(struc_count = n) %>%
  select(-n) 


# now I need to add up all aa within each (normaized!!!) bin to get bin totals and append this to the aa_counts
for_barplot_1 <- struc_counts1 %>%
  ungroup() %>%
  mutate(total_count = sum(struc_count)) 


for_barplot_2 <- for_barplot_1 %>%
  mutate(freq = struc_count/total_count) %>%
  select(-c(struc_count, total_count)) %>%
  mutate(second_struc = fct_rev(fct_relevel(second_struc, "Alpha helix", "Random coil", "Turn", "Beta strand", "Noncanon. helix", "Bridge")))

struc_fills = c("#bf8040", "#ac5396", "#70adc2", "#748f3d", "#cc5933", "#7070c2")

struc_mispr_plot <- for_barplot_2 %>%
  ggplot(aes(x = freq, y = second_struc, fill = second_struc)) +
  geom_col(alpha = 0.9) +
  scale_fill_manual(
    values = struc_fills) +
  scale_x_continuous(
    name = "Frequency of secondary structure \n at mispredicted positions",
    limits = c(0.0, 0.34),
    breaks = seq(0.0, 0.35, by = 0.05),
    expand = c(0, 0)) + 
  scale_y_discrete(
    name = "Secondary structure",
    expand = c(0.03, 0.03)) + 
  theme_cowplot(16) +
  theme(
    axis.text = element_text(color = "black", size = 14),
    strip.text.x = element_text(size = 16),
    panel.grid.major.x = element_line(color = "grey92", size=0.5),
    panel.grid.minor.x = element_line(color = "grey92", size=0.5),
    panel.spacing = unit(2, "lines"),
    legend.position = "none")

struc_mispr_plot

ggsave(filename = "./analysis/figures/struc_dist_mispr_engineer.png", plot = struc_mispr_plot, width = 6, height = 5)

#--------------------------------------------------------------------------------------
#now, lets looks at the distribution of correct predictions

cor_pred <- final_2 %>%
  filter(predicted == wt)

with_struc2 <- cor_pred %>%
  inner_join(sec_struc) %>%
  mutate(second_struc = ifelse(second_struc == "Coil", "Random coil", second_struc),
         second_struc = ifelse(second_struc == "Strand", "Beta strand", second_struc),
         second_struc = ifelse(second_struc == "AlphaHelix", "Alpha helix", second_struc),
         second_struc = ifelse(second_struc == "310Helix" | second_struc == "PiHelix", "Noncanon. helix", second_struc))

# finds the count of each aa acid per bin:
struc_counts2 <- with_struc2 %>%
  select(c(gene, position, second_struc)) %>%
  count(second_struc) %>%
  mutate(struc_count = n) %>%
  select(-n) 


# now I need to add up all aa within each (normaized!!!) bin to get bin totals and append this to the aa_counts
for_barplot3 <- struc_counts2 %>%
  ungroup() %>%
  mutate(total_count = sum(struc_count)) 


for_barplot_3 <- for_barplot3 %>%
  mutate(freq = struc_count/total_count) %>%
  select(-c(struc_count, total_count)) %>%
  mutate(second_struc = fct_rev(fct_relevel(second_struc, "Alpha helix", "Random coil", "Turn", "Beta strand", "Noncanon. helix", "Bridge")))

struc_fills = c("#bf8040", "#ac5396", "#70adc2", "#748f3d", "#cc5933", "#7070c2")

cor_pred_plot <- for_barplot_3 %>%
  ggplot(aes(x = freq, y = second_struc, fill = second_struc)) +
  geom_col(alpha = 0.9) +
  scale_fill_manual(
    values = struc_fills) +
  scale_x_continuous(
    name = "Frequency of secondary structure \n where wild type was predicted",
    limits = c(0.0, 0.34),
    breaks = seq(0.0, 0.35, by = 0.05),
    expand = c(0, 0)) + 
  scale_y_discrete(
    name = "Secondary structure",
    expand = c(0.03, 0.03)) + 
  theme_cowplot(16) +
  theme(
    axis.text = element_text(color = "black", size = 14),
    strip.text.x = element_text(size = 16),
    panel.grid.major.x = element_line(color = "grey92", size=0.5),
    panel.grid.minor.x = element_line(color = "grey92", size=0.5),
    panel.spacing = unit(2, "lines"),
    legend.position = "none")

cor_pred_plot

ggsave(filename = "./analysis/figures/struc_cor_predictions.png", plot = cor_pred_plot, width = 6, height = 5)

figure_final <- plot_grid(struc_mispr_plot, cor_pred_plot, ncol = 2, align = "h", labels = c('a', 'b'), rel_widths = c(1, 1))

ggsave(filename = paste("./analysis/figures/mispr_to_cor_pred.png"), plot = figure_final, width = 13, height = 6)

#----------------------------------------------------------------------------------------
# getting the odds ratio of mispredictions to correct predictions (secondary structures)
mispredictions <- for_barplot_2 %>%
  rename(freq_mispr = freq)

cor_predictions <- for_barplot_3 %>%
  rename(freq_cor = freq)

for_odds <- inner_join(mispredictions, cor_predictions)

for_odds2 <- for_odds %>%
  mutate(odds = freq_mispr/freq_cor)

plot_odds <- for_odds2 %>%
  ggplot(aes(x = odds, y = second_struc, fill = second_struc)) +
  geom_col(alpha = 0.9) +
  geom_text(aes(label = round(odds, 2), hjust = -0.05)) +
  scale_fill_manual(
    values = struc_fills) +
  scale_x_log10(
    name = "Odds ratio (log scale)"
  ) +
  scale_y_discrete(
    name = "Secondary structure",
    expand = c(0.03, 0.03)) + 
  theme_cowplot(16) +
  theme(
    axis.text = element_text(color = "black", size = 14),
    strip.text.x = element_text(size = 16),
    panel.grid.major.x = element_line(color = "grey92", size=0.5),
    panel.grid.minor.x = element_line(color = "grey92", size=0.5),
    panel.spacing = unit(2, "lines"),
    legend.position = "none")

plot_odds

ggsave(filename = "./analysis/figures/struc_odds.png", plot = plot_odds, width = 9, height = 5)


#doing some tests with KL-div
get_D_KL <- function(y, x) {
  sum <- 0
  for (i in 1:20) { 
    y[i] = y[i] + 1e-20 # the KL-div statistic can't be calculated when y is as zero
    x[i] = x[i] + 1e-20
    if (x[i] != 0 & y[i] != 0) {
      sum = sum + x[i] * log(x[i]/y[i])
    }
    if (x[i] != 0 & y[i] == 0) {
      sum = sum + x[i] * log(x[i]/y[i])
    }
  }
  return(sum)
}


