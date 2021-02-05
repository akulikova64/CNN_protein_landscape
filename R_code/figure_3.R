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
# figure 3. 

# heat maps. 

# loading data
cnn_data <- read.csv(file = "./cnn_wt_max_freq.csv", header=TRUE, sep=",")
natural_data <- read.csv(file = "./natural_max_freq.csv", header=TRUE, sep=",")

joined_data <- rbind(x = cnn_data, y = natural_data)

for_heat <- joined_data %>%
  select(group, aa, position, gene) %>%
  pivot_wider(names_from = group, values_from = aa) %>%
  group_by(wt, predicted) %>%
  summarise(count = n())

for_heat <- na.omit(for_heat)

custom_colors <- c("#9875bd", "#ecb613")

a <- for_heat %>%
  ggplot(aes(x = wt, y = predicted, fill = count)) +
  geom_tile() +
  scale_fill_gradient(low = "#e5ddee", high = "#3e2b55") +
  xlab("WT Residue") +
  ylab("Predicted Residue") +
  theme(panel.background = element_blank()) +
  theme_cowplot()

for_heat_2 <- joined_data %>%
  select(group, aa_class, position, gene) %>%
  pivot_wider(names_from = group, values_from = aa_class) %>%
  group_by(wt, predicted) %>%
  summarise(count = n())

for_heat_2 <- na.omit(for_heat_2)

a2 <- for_heat_2 %>%
  ggplot(aes(x = wt, y = predicted, fill = count)) +
  geom_tile() +
  scale_fill_gradient(low = "#fdf8e8", high = "#463606") +
  xlab("WT Residue Class") +
  ylab("Predicted Residue Class") +
  theme(panel.background = element_blank()) +
  theme_cowplot()

for_heat_3 <- joined_data %>%
  select(group, aa, position, gene) %>%
  pivot_wider(names_from = group, values_from = aa) %>%
  group_by(natural_max, predicted) %>%
  summarise(count = n())

for_heat_3 <- na.omit(for_heat_3)

b <- for_heat_3 %>%
  ggplot(aes(x = natural_max, y = predicted, fill = count)) +
  geom_tile() +
  scale_fill_gradient(low = "#e5ddee", high = "#3e2b55") +
  xlab("Consensus Residue") +
  ylab("Predicted Residue") +
  theme(panel.background = element_blank()) +
  theme_cowplot()

for_heat_4 <- joined_data %>%
  select(group, aa_class, position, gene) %>%
  pivot_wider(names_from = group, values_from = aa_class) %>%
  group_by(natural_max, predicted) %>%
  summarise(count = n())

for_heat_4 <- na.omit(for_heat_4)

b2 <- for_heat_4 %>%
  ggplot(aes(x = natural_max, y = predicted, fill = count)) +
  geom_tile() +
  scale_fill_gradient(low = "#fdf8e8", high = "#463606") +
  xlab("Consensus Residue Class") +
  ylab("Predicted Residue Class") +
  theme(panel.background = element_blank()) +
  theme_cowplot()

#==============================================================
# making and saving the plot
#==============================================================

figure_3a <- plot_grid(a, a2, nrow = 1, align="h", labels = c('A', 'B'))
figure_3b <- plot_grid(b, b2, nrow = 1, align="h", labels = c('C', 'D'))

ggsave(filename = "figure_3a.png", plot = figure_3a, width = 14, height = 5.5)
ggsave(filename = "figure_3b.png", plot = figure_3b, width = 14, height = 5.5)


figure_3 <- plot_grid(a, a2, b, b2, nrow = 2, align="hv", labels = c('A', 'B', 'C', 'D'))
ggsave(filename = "figure_3.png", plot = figure_3, width = 14, height = 9)








