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

# fix the code below
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

  
