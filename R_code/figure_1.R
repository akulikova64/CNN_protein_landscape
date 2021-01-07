library(tidyverse)
library(yardstick)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(cowplot)
library(tidyr)
library(sqldf)
library(sinaplot)
library(ggforce)

# code for figure 1. 
#==============================================================================
# a) Basic stats (% match across all genes)
#===============================================================================

# loading data
cnn_data <- read.csv(file = "./cnn_wt_max_freq.csv", header=TRUE, sep=",")
natural_data <- read.csv(file = "./natural_max_freq.csv", header=TRUE, sep=",")

joined_data <- rbind(x = cnn_data, y = natural_data)

# check if predicted aa is the one found in the wt structure
matches_1 <- joined_data %>%
  pivot_wider(names_from = group, values_from = c(aa, freq, aa_class, class_freq)) %>%
  mutate(match = aa_predicted == aa_wt) %>%
  select(gene, position, match)

# restart r
library(dplyr)

# selecting the data entries where the predicted amino acid matches the 
summary_stats_1 <- matches_1 %>%
  group_by(gene, match) %>%
  summarise(count = n()) %>%
  mutate(freq = count / sum(count))
  
stats <- summary_stats_1 %>%
  filter(match == TRUE) %>%
  mutate(group = "predicted aa = wt aa")

# now, we need to add accuracy within each class (predicted class == wt class)
matches_2 <- joined_data %>%
  select(c(gene, position, group, aa_class, class_freq)) %>%
  pivot_wider(names_from = group, values_from = c(aa_class, class_freq)) %>%
  mutate(match = aa_class_predicted == aa_class_wt) %>%
  select(gene, position, match)

# selecting the data entries where the predicted amino acid class matches the wt aa class.
library(dplyr)
stats_2 <- matches_2 %>%
  group_by(gene, match) %>%
  summarise(count = n()) %>%
  mutate(freq = count / sum(count)) %>%
  filter(match == TRUE) %>%
  mutate(group = "predicted class = wt class")

library(tidyverse)
library(yardstick)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(cowplot)
library(tidyr)
library(sqldf)
library(sinaplot)
library(ggforce)

joined_single_and_class <- rbind(stats, stats_2)

custom_colors <- c("#9875bd", "#ecb613")


a <- joined_single_and_class %>%
  ggplot(aes(x = freq, fill = group)) +
  geom_density(alpha = 0.5) + 
  scale_colour_manual(values = custom_colors, aesthetics = c("colour", "fill")) +
  #ggtitle(label = "CNN Accuracy by Protein") +
  theme_cowplot() + 
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) +
  ylab("Count") +
  xlab("Accuracy (predicting wt)") +
  coord_cartesian(xlim = c(0.06, 0.9), ylim = c(0, 8)) +
  scale_x_continuous(breaks = seq(0.1, 0.9, by = 0.1)) +
  scale_y_continuous(breaks = seq(0, 8, by = 2))

stats %>%
  summarise(mean = mean(freq))
# mean = 0.713

stats_2 %>%
  summarise(mean = mean(freq))
# mean = 0.799

#===================================================================================
# b) check if predicted aa is the most frequent aa in the alignment (consensus)
#===================================================================================
matches_3 <- joined_data %>%
  pivot_wider(names_from = group, values_from = c(aa, freq, aa_class, class_freq)) %>%
  mutate(match = aa_predicted == aa_natural_max) %>%
  select(gene, position, match)

# restart r
library(dplyr)
stats_3 <- matches_3 %>%
  group_by(gene, match) %>%
  summarise(count = n()) %>%
  mutate(freq = count / sum(count)) %>%
  filter(match == TRUE) %>%
  mutate(group = "predicted aa = consensus")

library(tidyr)

matches_4 <- joined_data %>%
  select(c(gene, position, group, aa_class, class_freq)) %>%
  pivot_wider(names_from = group, values_from = c(aa_class,class_freq)) %>%
  mutate(match = aa_class_predicted == aa_class_natural_max) %>%
  select(gene, position, match)
  
stats_4 <- matches_4 %>%
  group_by(gene, match) %>%
  summarise(count = n()) %>%
  mutate(freq = count / sum(count)) %>%
  filter(match == TRUE) %>%
  mutate(group = "predicted class = consensus")

library(tidyverse)
library(yardstick)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(cowplot)
library(tidyr)
library(sqldf)
library(sinaplot)
library(ggforce)

joined_consensus <- rbind(stats_3, stats_4)

b <- joined_consensus %>%
  ggplot(aes(x = freq, fill = group)) +
  geom_density(alpha = 0.5) + 
  scale_colour_manual(values = custom_colors, aesthetics = c("colour", "fill")) +
  #ggtitle(label = "CNN Predictions Compared to \n Alignment Consensus") +
  theme_cowplot() + 
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) +
  ylab("Count") +
  xlab("Accuracy (predicting alignment consensus)") +
  coord_cartesian(xlim = c(0.06, 0.9), ylim = c(0, 8)) +
  scale_x_continuous(breaks = seq(0.1, 0.9, by = 0.1)) +
  scale_y_continuous(breaks = seq(0, 8, by = 2))

summary_stats_2 %>%
  summarise(mean = mean(freq))
# mean = 0.331

#==============================================================
# making and saving the plot
#==============================================================

plot_grid(a, b, nrow = 2, align="v", labels = c('A', 'B'))
