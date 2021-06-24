library(tidyverse)
library(cowplot)

# figure 3. 

# heat maps. 

font_size <- 10

# set working directory to: "Desktop/Natural_var_project/"
# loading data
cnn_data <- read.csv(file = "./data/PSICOV_box_20/output/cnn_wt_max_freq.csv", header=TRUE, sep=",")
natural_data <- read.csv(file = "./data/PSICOV_box_20/output/natural_max_freq_files/natural_max_freq_all.csv", header=TRUE, sep=",")

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
    freq = count/sum(count)) %>%
    #freq = ifelse(wt == predicted, NA, freq)) %>%
  ungroup()
  
for_heat_sums2 <- na.omit(for_heat_sums2)

custom_colors <- c("#9875bd", "#ecb613")

for_heatplot_final <- for_heat_sums2 %>%
  mutate(
    predicted = fct_rev(fct_relevel(predicted, "M","L","I","V","A","C","S","T","N","Q","D","E","R","K","H","Y","F","W","P","G")),
    wt = fct_relevel(wt, "M","L","I","V","A","C","S","T","N","Q","D","E","R","K","H","Y","F","W","P","G")) 
 
#add aa class label: 
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

for_heatplot_with_classes <- for_heatplot_final %>%
  mutate(wt_class = map_chr(wt, calc_class)) %>%
  mutate(class = fct_relevel(wt_class, "aliphatic", "small_polar", "negative", "positive", "aromatic", "unique"))
  
gray_zone = tibble(
  x = c("G","A","V","M","I","L","S","C","N","T","Q","D","E","H","K","R","F","Y","W","P"), 
  y = c("G","A","V","M","I","L","S","C","N","T","Q","D","E","H","K","R","F","Y","W","P"), 
  value = rep(1, times=20))

plot_20_a1 <- ggplot() +
  geom_tile(data = for_heatplot_with_classes, aes(x = wt, y = predicted, alpha = freq, fill = class)) + 
  scale_alpha_continuous(
    guide = guide_legend(order = 2, reverse = TRUE),
    range = c(0.2, 1)) +
  scale_fill_manual(
    values = c("#991f00", "#001a66", "#994d00", "#1a6600", "#330066", "#9e9e2e"),
    guide = guide_legend(order = 1)) +
  #geom_tile(data = gray_zone, aes(x,y), fill = "grey57") +
  scale_x_discrete(
    name = "Wild type residue",
    expand = c(0,0)) +
  scale_y_discrete(
    name = "Predicted residue",
    expand = c(0,0)) +
  labs(fill = "Wild type class", alpha = "Frequency") +
  theme_cowplot(font_size) +
  theme(
    panel.background = element_blank(),
    axis.text = element_text(color = "black", size = font_size),
    legend.position = "none",
    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")) 


plot_20_a1


colors <-  c("#991f00", "#001a66", "#994d00", "#1a6600", "#330066", "#9e9e2e")
alphas <- c(0.2, 0.4, 0.6, 0.8, 1.0)
class_list <- c("aliphatic", "small_polar", "negative", "positive", "aromatic", "unique")

twoD_legend <- tibble(class = c(rep("aliphatic", times = 5), 
                                rep("small_polar", times = 5), 
                                rep("negative", times = 5), 
                                rep("positive", times = 5), 
                                rep("aromatic", times = 5), 
                                rep("unique", times = 5)), 
                      freq = c(rep(alphas, times = 6)))

twoD_legend <- twoD_legend %>%
  mutate(class = fct_relevel(class, "aliphatic", "small_polar", "negative", "positive", "aromatic", "unique"))

legend <- twoD_legend %>%
  ggplot(aes(x = class, y = factor(freq), fill = class, alpha = factor(freq))) +
  geom_tile() +
  scale_fill_manual(
    values = c("#991f00", "#001a66", "#994d00", "#1a6600", "#330066", "#9e9e2e"),
    guide = guide_legend(order = 1)) +
  scale_alpha_manual(
    guide = guide_legend(order = 2, reverse = TRUE),
    values = c(0.2, 0.4, 0.6, 0.8, 1.0)) +
  scale_x_discrete(
    name = "Amino acid class",
    position = "top",
    label = c("aliphatic", "small polar", "negative", "positive", "aromatic", "unique"),
    expand = c(0,0)) +
  scale_y_discrete(
    name = "Frequency \n",
    position = "right",
    labels = rev(c("1.0", "0.8", "0.6", "0.4", "0.2")),
    expand = c(0,0)) +
  theme_minimal() +
  theme(
    legend.position = "none",
    panel.grid = element_blank(),
    axis.title = element_text(
      size = font_size,
      color = "black"),
    axis.text.x = element_text(
      size = font_size,
      color = "black",
      angle = 45, 
      hjust = 0.08),
    axis.text.y = element_text(
      size = font_size,
      color = "black"),
    aspect.ratio = 5/5)

    
legend

ggsave(filename = "./analysis/figures/legend_fig_3.png", plot = legend, width = 1, height = 1)

# by classes:

for_heat_2 <- joined_data %>%
  select(group, aa_class, position, gene) %>%
  pivot_wider(names_from = group, values_from = aa_class) 

for_heat_sums <- for_heat_2 %>%
  group_by(wt, predicted) %>%
  summarise(count = n())

for_heat_sums2 <- for_heat_sums %>%
  group_by(wt) %>%
  mutate(freq = count/sum(count)) %>%
         #freq = ifelse(wt == predicted, NA, freq)) %>%
  ungroup()


for_heat_sums2 <- na.omit(for_heat_sums2)

for_plot_a2 <- for_heat_sums2 %>%
  mutate(
    wt = fct_relevel(wt, "aliphatic", "small_polar", "negative", "positive", "aromatic", "unique"),
    predicted = fct_rev(fct_relevel(predicted, "aliphatic", "small_polar", "negative", "positive", "aromatic", "unique"))
  )

gray_zone2 = tibble(
  x = c("aliphatic", "small_polar", "negative", "positive", "aromatic", "unique"), 
  y = c("aliphatic", "small_polar", "negative", "positive", "aromatic", "unique"), 
  value = rep(1, times=6))

plot_20_a2 <- ggplot() +
  geom_tile(data = for_plot_a2, aes(x = wt, y = predicted, alpha = freq, fill = wt)) +
  #geom_tile(data = gray_zone2, aes(x,y), fill = "grey57") +
  scale_alpha_continuous(
    guide = guide_legend(order = 2, reverse = TRUE),
    range = c(0.2, 1)) +
  scale_fill_manual(
    values = c("#991f00", "#001a66", "#994d00", "#1a6600", "#330066", "#9e9e2e"),
    guide = guide_legend(order = 1) ) +
  scale_x_discrete(
    name = "Wild type residue class",
    label = c("aliphatic", "small polar", "negative", "positive", "aromatic", "unique"),
    expand = c(0,0)) +
  scale_y_discrete(
    name = "Predicted residue class",
    label = c("aliphatic", "small polar", "negative", "positive", "aromatic", "unique"),
    expand = c(0,0)) +
  labs(fill = "Wild type class", alpha = "Frequency") +
  theme_cowplot(font_size) +
  theme(
    panel.background = element_blank(),
    axis.text = element_text(
      color = "black", 
      size = font_size),
    axis.text.x = element_text(
      angle = 45, 
      hjust = 1),
    legend.position = "none",
    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))

plot_20_a2

figure_3a <- plot_grid(plot_20_a1, plot_20_a2, nrow = 1, align = "h", labels = c('A', 'B'), label_x = 0)

ggsave(filename = "./analysis/figures/figure_3a.png", plot = figure_3a, width = 10, height = 4.5)


#===========================================================================================
# now heat maps of predicting the consensus
#===========================================================================================

for_heat_3 <- joined_data %>%
  select(group, aa_class, position, gene) %>%
  pivot_wider(names_from = group, values_from = aa_class) 

for_heat_sums <- for_heat_2 %>%
  group_by(natural_max, predicted) %>%
  summarise(count = n())

for_heat_sums2 <- for_heat_sums %>%
  group_by(natural_max) %>%
  mutate(freq = count/sum(count)) %>%
  ungroup()

for_heat_sums2 <- na.omit(for_heat_sums2)

for_plot_b2 <- for_heat_sums2 %>%
  mutate(
    natural_max = fct_relevel(natural_max, "aliphatic", "small_polar", "negative", "positive", "aromatic", "unique"),
    predicted = fct_rev(fct_relevel(predicted, "aliphatic", "small_polar", "negative", "positive", "aromatic", "unique"))
  )


plot_20_b2 <- for_plot_b2 %>%
  ggplot(aes(x = natural_max, y = predicted, alpha = freq, fill = natural_max)) +
  geom_tile() +
  scale_alpha_continuous(
    guide = guide_legend(order = 2, reverse = TRUE),
    range = c(0.2, 1)) +
  scale_fill_manual(
    values = c("#991f00", "#001a66", "#994d00", "#1a6600", "#330066", "#9e9e2e"),
    guide = 'none' ) +
  scale_x_discrete(
    name = "Consensus residue class",
    label = c("aliphatic", "small polar", "negative", "positive", "aromatic", "unique"),
    expand = c(0,0)) +
  scale_y_discrete(
    name = "Predicted residue class",
    label = c("aliphatic", "small polar", "negative", "positive", "aromatic", "unique"),
    expand = c(0,0)) +
  labs(fill = "Consensus \n class", alpha = "Frequency") +
  theme_cowplot(font_size) +
  theme(
    panel.background = element_blank(),
    axis.text = element_text(
      color = "black", 
      size = font_size),
    axis.text.x = element_text(
      angle = 45, 
      hjust = 1),
    legend.position = "none",
    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))

plot_20_b2


# predicting consensus heat map:
for_heat_4 <- joined_data %>%
  select(group, aa, position, gene) %>%
  pivot_wider(names_from = group, values_from = aa) 

for_heat_sums <- for_heat_4 %>%
  group_by(natural_max, predicted) %>%
  summarise(count = n())

for_heat_sums2 <- for_heat_sums %>%
  group_by(natural_max) %>%
  mutate(
    freq = count/sum(count)) %>%
    #freq = ifelse(natural_max == predicted, NA, freq)) %>%
  ungroup()

for_heat_4 <- na.omit(for_heat_sums2)

for_heatplot_final <- for_heat_4 %>%
  mutate(
    predicted = fct_rev(fct_relevel(predicted, "M","L","I","V","A","C","S","T","N","Q","D","E","R","K","H","Y","F","W","P","G")),
    wt = fct_relevel(natural_max, "M","L","I","V","A","C","S","T","N","Q","D","E","R","K","H","Y","F","W","P","G")) 

for_heatplot_with_classes <- for_heatplot_final %>%
  mutate(natural_max_class = map_chr(natural_max, calc_class)) %>%
  mutate(class = fct_relevel(natural_max_class, "aliphatic", "small_polar", "negative", "positive", "aromatic", "unique"))


plot_20_b1 <- for_heatplot_with_classes %>%
  ggplot(aes(
    x = fct_relevel(natural_max, "M","L","I","V","A","C","S","T","N","Q","D","E","R","K","H","Y","F","W","P","G"), 
    y = fct_rev(fct_relevel(predicted, "M","L","I","V","A","C","S","T","N","Q","D","E","R","K","H","Y","F","W","P","G")), 
                alpha = freq, 
                fill = class)) +
  geom_tile() + 
  scale_alpha_continuous(
    guide = guide_legend(order = 1, reverse = TRUE),
    range = c(0.2, 1)) +
  scale_fill_manual(
    values = c("#991f00", "#001a66", "#994d00", "#1a6600", "#330066", "#9e9e2e"),
    guide = 'none') +
  scale_x_discrete(
    name = "Consensus residue",
    expand = c(0,0)) +
  scale_y_discrete(
    name = "Predicted residue",
    expand = c(0,0)) +
  #labs(fill = "Consensus Class", alpha = "Frequency") +
  theme_cowplot(font_size) +
  theme(
    panel.background = element_blank(),
    axis.text = element_text(color = "black", size = font_size),
    legend.position = "none",
    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")) 


plot_20_b1

figure_3abcd <- plot_grid(plot_20_a1, plot_20_a2, plot_20_b1, plot_20_b2, rel_widths = c(1, 1, 1, 1), nrow = 2, axis = "btrl", labels = c('A', 'B', 'C', 'D'), label_x = 0)
ggsave(filename = "./analysis/figures/figure_3abcd.png", plot = figure_3abcd, width = 6.5, height = 6.5)


#figure_3_almost <- plot_grid(figure_3a, figure_3b, nrow = 2, align="h", labels = c('', ''))
#figure_3_almost

figure_3 <- plot_grid(figure_3abcd, legend, nrow = 1, rel_widths = c(3, 1), scale = c(1, 0.80), labels = c('', ''))

ggsave(filename = "./analysis/figures/figure_3_box_20.png", plot = figure_3, width = 9, height = 6.5)








#Now I want to use the 80-100% homology group:

cnn_data <- read.csv(file = "./cnn_wt_max_freq.csv", header=TRUE, sep=",")
natural_data <- read.csv(file = "./natural_max_freq_files/natural_max_freq_20.csv", header=TRUE, sep=",")

joined_data <- rbind(x = cnn_data, y = natural_data)

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

for_heatplot_final <- for_heat_4 %>%
  mutate(
    predicted = fct_rev(fct_relevel(predicted, "M","L","I","V","A","C","S","T","N","Q","D","E","R","K","H","Y","F","W","P","G")),
    natural_max = fct_relevel(natural_max, "M","L","I","V","A","C","S","T","N","Q","D","E","R","K","H","Y","F","W","P","G")) %>%
  mutate(natural_class = map_chr(natural_max, calc_class)) %>%
  mutate(natural_class = fct_relevel(natural_class, "aliphatic", "small_polar", "negative", "positive", "aromatic", "unique"))


plot_c <- for_heatplot_final %>%
  ggplot(aes(
    x = fct_relevel(natural_max, "M","L","I","V","A","C","S","T","N","Q","D","E","R","K","H","Y","F","W","P","G"), 
    y = fct_rev(fct_relevel(predicted, "M","L","I","V","A","C","S","T","N","Q","D","E","R","K","H","Y","F","W","P","G")), 
    alpha = freq, 
    fill = natural_class)) +
  geom_tile() + 
  scale_alpha_continuous(
    guide = guide_legend(order = 2, reverse = TRUE),
    range = c(0.2, 2)) +
  scale_fill_manual(
    values = c("#991f00", "#001a66", "#994d00", "#1a6600", "#330066", "#9e9e2e"),
    guide = guide_legend(order = 1) ) +
  scale_x_discrete(
    name = "Consensus Residue",
    expand = c(0,0)) +
  scale_y_discrete(
    name = "Predicted Residue",
    expand = c(0,0)) +
  labs(fill = "WT Class", alpha = "Frequency") +
  theme_cowplot(12) +
  theme(
    panel.background = element_blank(),
    axis.text = element_text(color = "black", size = 12))


plot_c






#==============================================================
# making and saving the plot
#==============================================================

figure_3a <- plot_grid(plot_1a, plot_a2, nrow = 1, align="h", labels = c('A', 'B'), label_x = 0)
figure_3b <- plot_grid(plot_b, plot_b2, nrow = 1, align="h", labels = c('C', 'D'), label_x = 0)

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





