library(tidyverse)
library(cowplot)

# figure 3. 

# heat maps. 

# loading data
cnn_data <- read.csv(file = "./cnn_wt_max_freq.csv", header=TRUE, sep=",")
natural_data <- read.csv(file = "./natural_max_freq.csv", header=TRUE, sep=",")

joined_data <- rbind(x = cnn_data, y = natural_data)

for_heat <- joined_data %>%
  select(group, aa, position, gene) %>%
  pivot_wider(names_from = group, values_from = aa) 

for_heat_sums <- for_heat %>%
  group_by(wt, predicted) %>%
  summarise(count = n())
  
for_heat_sums2 <- for_heat_sums %>%
  group_by(wt) %>%
  mutate(
    freq = count/sum(count),
    freq = ifelse(wt == predicted, NA, freq)) %>%
  ungroup()
  
for_heat_sums2 <- na.omit(for_heat_sums2)

custom_colors <- c("#9875bd", "#ecb613")

for_heatplot_final <- for_heat_sums2 %>%
  mutate(
    predicted = fct_relevel(predicted, "G","A","V","M","I","L","S","C","N","T","Q","D","E","H","K","R","F","Y","W","P"),
    wt = fct_relevel(wt, "G","A","V","M","I","L","S","C","N","T","Q","D","E","H","K","R","F","Y","W","P")) 
 
#add aa class label: 
calc_class <- function(x) {
  
  aliphatic = c("G", "A", "V", "M", "I", "L")
  aromatic = c("F", "Y", "W")
  polar = c("S", "C", "N", "T", "Q")
  negative = c("D", "E")
  positive = c("H", "K", "R")
    
  if (x %in% aliphatic) {
    return("aliphatic")
  }
  if (x %in% aromatic) {
    return("aromatic")
  }
  if (x %in% polar) {
    return("polar")
  }
  if (x %in% negative) {
    return("negative")
  }
  if (x %in% positive) {
    return("positive")
  }
  if (x == "P") {
    return("proline")
  }
  
}

for_heatplot_with_classes <- for_heatplot_final %>%
  mutate(wt_class = map_chr(wt, calc_class))
  
plot_a1 <- for_heatplot_with_classes %>% 
  ggplot(aes(x = wt, y = predicted, alpha = freq, fill = wt_class)) +
  geom_tile() + 
  #scale_fill_gradient(low = "#e5ddee", high = "#3e2b55") +
  scale_alpha_continuous(guide = guide_legend(order = 2)) +
  scale_fill_manual(values = c("red", "green", "blue", "yellow", "magenta", "cyan"), guide = guide_legend(order = 1) ) +
  xlab("WT Residue") +
  ylab("Predicted Residue") +
  theme(panel.background = element_blank()) +
  theme_cowplot()

plot_a1

# by classes:

for_heat_2 <- joined_data %>%
  select(group, aa_class, position, gene) %>%
  pivot_wider(names_from = group, values_from = aa_class) 

for_heat_sums <- for_heat_2 %>%
  group_by(wt, predicted) %>%
  summarise(count = n())

for_heat_sums2 <- for_heat_sums %>%
  group_by(wt) %>%
  mutate(freq = count/sum(count),
         freq = ifelse(wt == predicted, NA, freq)) %>%
  ungroup()


for_heat_sums2 <- na.omit(for_heat_sums2)

for_plot_a2 <- for_heat_sums2 %>%
  mutate(
    wt = fct_relevel(wt, "aliphatic", "polar", "negative", "positive", "aromatic", "proline"),
    predicted = fct_relevel(predicted, "aliphatic", "polar", "negative", "positive", "aromatic", "proline")
  )

a2 <- for_plot_a2 %>%
  ggplot(aes(x = wt, y = predicted, alpha = freq)) +
  geom_tile(fill = "#463606") +
  #scale_fill_gradient(low = "#fdf8e8", high = "#463606") +
  scale_alpha_continuous() +
  xlab("WT Residue Class") +
  ylab("Predicted Residue Class") +
  theme(panel.background = element_blank()) +
  theme_cowplot()

a2

# how heat maps of predicting the consensus

for_heat_3 <- joined_data %>%
  select(group, aa, position, gene) %>%
  pivot_wider(names_from = group, values_from = aa) 

for_heat_sums <- for_heat_3 %>%
  group_by(natural_max, predicted) %>%
  summarise(count = n())

for_heat_sums2 <- for_heat_sums %>%
  group_by(natural_max) %>%
  mutate(freq = count/sum(count))

for_heat_3 <- na.omit(for_heat_sums2)

plot_b <- for_heat_3 %>%
  ggplot(aes(x = natural_max, y = predicted, fill = freq)) +
  geom_tile() +
  scale_fill_gradient(low = "#e5ddee", high = "#3e2b55") +
  xlab("Consensus Residue") +
  ylab("Predicted Residue") +
  theme(panel.background = element_blank()) +
  theme_cowplot()

plot_b

# predicting consensus class heat map:
for_heat_4 <- joined_data %>%
  select(group, aa_class, position, gene) %>%
  pivot_wider(names_from = group, values_from = aa_class) 

for_heat_sums <- for_heat_4 %>%
  group_by(natural_max, predicted) %>%
  summarise(count = n())

for_heat_sums2 <- for_heat_sums %>%
  group_by(natural_max) %>%
  mutate(freq = count/sum(count))

for_heat_4 <- na.omit(for_heat_sums2)

plot_b2 <- for_heat_4 %>%
  ggplot(aes(x = natural_max, y = predicted, fill = freq)) +
  geom_tile() +
  scale_fill_gradient(low = "#fdf8e8", high = "#463606") +
  xlab("Consensus Residue Class") +
  ylab("Predicted Residue Class") +
  theme(panel.background = element_blank()) +
  theme_cowplot()

plot_b2
#==============================================================
# making and saving the plot
#==============================================================

figure_3a <- plot_grid(plot_1a, a2, nrow = 1, align="h", labels = c('A', 'B'))
figure_3b <- plot_grid(plot_b, plot_b2, nrow = 1, align="h", labels = c('C', 'D'))

ggsave(filename = "../../analysis/figures/figure_3a_new.png", plot = figure_3a, width = 14, height = 5.5)
ggsave(filename = "../../analysis/figures/figure_3b_new.png", plot = figure_3b, width = 14, height = 5.5)


figure_3 <- plot_grid(a, a2, b, b2, nrow = 2, align="hv", labels = c('A', 'B', 'C', 'D'))
ggsave(filename = "figure_3.png", plot = figure_3, width = 14, height = 9)


#tests ...

table <- tibble(y = c(5, 3, 1, 3, 7), x = c(5, 3, 4, 3, 1))

table2 <- table %>%
  count(x, y) %>%
  mutate(freq = n/sum(n))

# it works 





