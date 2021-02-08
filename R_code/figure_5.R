library(tidyverse)
library(cowplot)
library(ggforce)

# loading data
cnn_data <- read.csv(file = "./cnn_wt_max_freq.csv", header=TRUE, sep=",")
natural_data_20 <- read.csv(file = "./natural_max_freq_files/natural_max_freq_20.csv", header=TRUE, sep=",")
natural_data_40 <- read.csv(file = "./natural_max_freq_files/natural_max_freq_40.csv", header=TRUE, sep=",")
natural_data_60 <- read.csv(file = "./natural_max_freq_files/natural_max_freq_60.csv", header=TRUE, sep=",")
natural_data_80 <- read.csv(file = "./natural_max_freq_files/natural_max_freq_80.csv", header=TRUE, sep=",")
natural_data_100 <- read.csv(file = "./natural_max_freq_files/natural_max_freq_100.csv", header=TRUE, sep=",")


joined_data_20 <- rbind(x = cnn_data, y = natural_data_20)
joined_data_40 <- rbind(x = cnn_data, y = natural_data_40)
joined_data_60 <- rbind(x = cnn_data, y = natural_data_60)
joined_data_80 <- rbind(x = cnn_data, y = natural_data_80)
joined_data_100 <- rbind(x = cnn_data, y = natural_data_100)

joined_data_20 <- joined_data_20 %>%
  mutate(perc_sim = "(0-20%]") 

joined_data_40 <- joined_data_40 %>%
  mutate(perc_sim = "(20-40%]")

joined_data_60 <- joined_data_60 %>%
  mutate(perc_sim = "(40-60%]")

joined_data_80 <- joined_data_80 %>%
  mutate(perc_sim = "(60-80%]")

joined_data_100 <- joined_data_100 %>%
  mutate(perc_sim = "(80-100%]")

all_data <- rbind(joined_data_20, joined_data_40, joined_data_60, joined_data_80, joined_data_100)

match_consensus <- all_data %>%
  pivot_wider(names_from = group, values_from = c(aa, freq, aa_class, class_freq)) %>%
  mutate(match = aa_predicted == aa_natural_max) %>%
  select(gene, position, match, perc_sim)


stats_for_plot <- match_consensus %>%
  group_by(gene, match, perc_sim) %>%
  summarise(count = n()) %>%
  mutate(freq = count / sum(count)) %>%
  filter(match == TRUE) %>%
  mutate(group = "predicted aa = consensus")

stats_for_plot %>%
  ggplot(aes(y = freq, x = perc_sim)) +
  geom_violin(alpha = 0.5) + 
  geom_sina(size = 0.75) +
  theme_cowplot() + 
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5), legend.position = "none") +
  ylab("Accuracy") +
  xlab("") 












