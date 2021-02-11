library(tidyverse)
library(yardstick)
library(ggpubr)
library(cowplot)
library(sqldf)
library(sinaplot)
library(dplyr)
library(ggforce)

# reading csv files
natural_var <- read.csv(file="./stats_align_all.csv", header=TRUE, sep=",")
cnn_data <- read.csv(file = "./cnn_wt_max_freq.csv", header=TRUE, sep=",")

# joining the two data frames.

natural_var <- natural_var %>%
  nest(q_natural = c(q_H:q_C))

joined_data <- inner_join(x = natural_var, y = cnn_data, by = c('position', 'gene'))

# Figure 4
# wt prob (CNN) as a function of neff in alignment.

n_eff_data_1 <- joined_data %>%
  select(n_eff, gene, position)
names(n_eff_data_1) <- c('n_eff', 'gene', 'position')
n_eff_data_1$group <- 'n_eff_natural'

data_2 <- joined_data %>%
  filter(group == 'wt') %>%
  select(freq, gene, position)
data_2$group <- 'freq_predicted'

n_eff_data <- inner_join(x = n_eff_data_1, y = data_2, by = c('gene', 'position'))
n_eff_unique <- unique(n_eff_data)

n_eff_averaged <- n_eff_unique %>%
  group_by(gene, position, group.x, group.y) %>%
  summarise(n_eff = mean(n_eff), freq = mean(freq))

library(dplyr)
# finding all the corelations of all proteins
cor_1 <- n_eff_averaged %>%
  select(gene, n_eff, freq) %>%
  na.omit() %>%
  group_by(gene) %>%
  summarise(cor = cor(n_eff, freq)) %>%
  mutate(group = "by aa") %>%
  mutate(R2 = cor^2)

goi_1 = "1i71p"
goi_2 = "1rw7"
goi_3 = "1g2r"

#finding the linear regression
linear_reg_1 <- n_eff_averaged %>%
  filter(gene == goi_1) %>%
  select(n_eff, freq) %>%
  lm(formula = n_eff ~ freq)
summary(linear_reg_1)

slope_1 <- linear_reg_1$coefficients[2]
intercept_1 <- linear_reg_1$coefficients[1]

purple <- "#9875bd"
mustard <- "#ecb613"

# scatter plot
a <- n_eff_averaged %>%
  filter(gene == goi_1) %>%
  ggplot(aes(x = n_eff, y = freq)) +
  geom_point(color = purple) +
  ggtitle(label = goi_1) +
  #geom_abline(slope = slope_1, intercept = intercept_1, color = "grey65" ) +
  stat_cor(label.y = 1.23, label.x = 8) +
  stat_regline_equation(label.y = 1.1, label.x = 8) +
  theme_cowplot() +
  scale_colour_manual(values = purple, aesthetics = c("colour", "fill")) +
  theme(plot.title = element_text(hjust=0.5), plot.subtitle = element_text(hjust=0.5)) +
  ylab("WT probability (CNN)") +
  xlab("N-eff Natural") 
  #coord_cartesian(xlim = c(1, 16), ylim = c(0, 1.5)) +
  #scale_x_continuous(breaks = seq(2, 16, by = 2)) +
  #scale_y_continuous(breaks = seq(0, 1.4 , by = 0.2))

# second gene:
linear_reg_2 <- n_eff_averaged %>%
  filter(gene == goi_2) %>%
  select(n_eff, freq) %>%
  lm(formula = n_eff ~ freq)
summary(linear_reg_2)

slope_2 <- linear_reg_2$coefficients[2]
intercept_2 <- linear_reg_2$coefficients[1]

# scatter plot
b <- n_eff_averaged %>%
  filter(gene == goi_2) %>%
  ggplot(aes(x = n_eff, y = freq)) +
  geom_point(color = purple) +
  ggtitle(label = goi_2) +
  #geom_abline(slope = slope_2, intercept = intercept_2, color = "grey65" ) +
  stat_cor(label.y = 1.23, label.x = 8) +
  stat_regline_equation(label.y = 1.1, label.x = 8) +
  theme_cowplot() +
  scale_colour_manual(values = purple, aesthetics = c("colour", "fill")) +
  theme(plot.title = element_text(hjust=0.5), plot.subtitle = element_text(hjust=0.5)) +
  ylab("WT probability (CNN)") +
  xlab("N-eff Natural") +
  coord_cartesian(xlim = c(1, 16), ylim = c(0, 1.5)) +
  scale_x_continuous(breaks = seq(2, 16, by = 2)) +
  scale_y_continuous(breaks = seq(0, 1.4 , by = 0.2))

# third gene:
linear_reg_3 <- n_eff_averaged %>%
  filter(gene == goi_3) %>%
  select(n_eff, freq) %>%
  lm(formula = n_eff ~ freq)
summary(linear_reg_3)

slope_3 <- linear_reg_3$coefficients[2]
intercept_3 <- linear_reg_3$coefficients[1]

# scatter plot
c <- n_eff_averaged %>%
  filter(gene == goi_3) %>%
  ggplot(aes(x = n_eff, y = freq)) +
  geom_point(color = purple) +
  ggtitle(label = goi_3) +
  #geom_abline(slope = slope_3, intercept = intercept_3, color = "grey65" ) +
  stat_cor(label.y = 1.23, label.x = 8) +
  stat_regline_equation(label.y = 1.1, label.x = 8) +
  theme_cowplot() +
  scale_colour_manual(values = purple, aesthetics = c("colour", "fill")) +
  theme(plot.title = element_text(hjust=0.5), plot.subtitle = element_text(hjust=0.5)) +
  ylab("WT probability (CNN)") +
  xlab("N-eff Natural") +
  coord_cartesian(xlim = c(1, 16), ylim = c(0, 1.5)) +
  scale_x_continuous(breaks = seq(2, 16, by = 2)) +
  scale_y_continuous(breaks = seq(0, 1.4 , by = 0.2))

#============================================================================================
# comparing the effective number of predicted classes to wt class probability (CNN)
#============================================================================================


class_data <- joined_data %>%
  select(n_eff, gene, position)
names(class_data) <- c('n_eff_class', 'gene', 'position')
class_data$group <- 'n_eff_natural'

class_data_2 <- joined_data %>%
  filter(group == 'wt') %>%
  select(class_freq, gene, position)
class_data_2$group <- 'freq_predicted'

n_eff_data <- inner_join(x = class_data, y = class_data_2, by = c('gene', 'position'))
n_eff_unique <- unique(n_eff_data)

n_eff_averaged <- n_eff_unique %>%
  group_by(gene, position, group.x, group.y) %>%
  summarise(n_eff_class = mean(n_eff_class), class_freq = mean(class_freq))

# finding all the corelations of all proteins
cor_2 <- n_eff_averaged %>%
  select(gene, n_eff_class, class_freq) %>%
  na.omit() %>%
  group_by(gene) %>%
  summarise(cor = cor(n_eff_class, class_freq)) %>%
  mutate(group = "by class") %>%
  mutate(R2 = cor^2)

#finding the linear regression
linear_reg_4 <- n_eff_averaged %>%
  filter(gene == goi_1) %>%
  select(n_eff_class, class_freq) %>%
  lm(formula = n_eff_class ~ class_freq)
summary(linear_reg_4)

slope_4 <- linear_reg_4$coefficients[2]
intercept_4 <- linear_reg_4$coefficients[1]

# scatter plot
a2 <- n_eff_averaged %>%
  filter(gene == goi_1) %>%
  ggplot(aes(x = n_eff_class, y = class_freq)) +
  geom_point(color = mustard) +
  ggtitle(label = goi_1) +
  #geom_abline(slope = slope_4, intercept = intercept_4, color = "grey65" ) +
  stat_cor(label.y = 1.23, label.x = 8) +
  stat_regline_equation(label.y = 1.1, label.x = 8) +
  theme_cowplot() +
  scale_colour_manual(values = mustard, aesthetics = c("colour", "fill")) +
  theme(plot.title = element_text(hjust=0.5), plot.subtitle = element_text(hjust=0.5)) +
  ylab("WT Class probability (CNN)") +
  xlab("N-eff Class Natural") +
  coord_cartesian(xlim = c(1, 16), ylim = c(0, 1.5)) +
  scale_x_continuous(breaks = seq(2, 16, by = 2)) +
  scale_y_continuous(breaks = seq(0, 1.4 , by = 0.2))

# second gene
#finding the linear regression
linear_reg_5 <- n_eff_averaged %>%
  filter(gene == goi_2) %>%
  select(n_eff_class, class_freq) %>%
  lm(formula = n_eff_class ~ class_freq)
summary(linear_reg_5)

slope_5 <- linear_reg_5$coefficients[2]
intercept_5 <- linear_reg_5$coefficients[1]

# scatter plot
b2 <- n_eff_averaged %>%
  filter(gene == goi_2) %>%
  ggplot(aes(x = n_eff_class, y = class_freq)) +
  geom_point(color = mustard) +
  ggtitle(label = goi_2) +
  #geom_abline(slope = slope_5, intercept = intercept_5, color = "grey65" ) +
  stat_cor(label.y = 1.23, label.x = 8) +
  stat_regline_equation(label.y = 1.1, label.x = 8) +
  theme_cowplot() +
  scale_colour_manual(values = mustard, aesthetics = c("colour", "fill")) +
  theme(plot.title = element_text(hjust=0.5), plot.subtitle = element_text(hjust=0.5)) +
  ylab("WT Class probability (CNN)") +
  xlab("N-eff Class Natural") +
  coord_cartesian(xlim = c(1, 16), ylim = c(0, 1.5)) +
  scale_x_continuous(breaks = seq(2, 16, by = 2)) +
  scale_y_continuous(breaks = seq(0, 1.4 , by = 0.2))

# third gene
#finding the linear regression
linear_reg_6 <- n_eff_averaged %>%
  filter(gene == goi_3) %>%
  select(n_eff_class, class_freq) %>%
  lm(formula = n_eff_class ~ class_freq)
summary(linear_reg_6)

slope_6 <- linear_reg_6$coefficients[2]
intercept_6 <- linear_reg_6$coefficients[1]

# scatter plot
c2 <- n_eff_averaged %>%
  filter(gene == goi_3) %>%
  ggplot(aes(x = n_eff_class, y = class_freq)) +
  geom_point(color = mustard) +
  ggtitle(label = goi_3) +
  #geom_abline(slope = slope_6, intercept = intercept_6, color = "grey65" ) +
  stat_cor(label.y = 1.23, label.x = 8) +
  stat_regline_equation(label.y = 1.1, label.x = 8) +
  theme_cowplot() +
  scale_colour_manual(values = mustard, aesthetics = c("colour", "fill")) +
  theme(plot.title = element_text(hjust=0.5), plot.subtitle = element_text(hjust=0.5)) +
  ylab("WT Class probability (CNN)") +
  xlab("N-eff Class Natural") +
  coord_cartesian(xlim = c(1, 16), ylim = c(0, 1.5)) +
  scale_x_continuous(breaks = seq(2, 16, by = 2)) +
  scale_y_continuous(breaks = seq(0, 1.4 , by = 0.2))

#=========================================================================================================
# density plot comparing the distribution of correlation coefficients between the two groups (aa vs class)
#=========================================================================================================

correlation_coeffs <- rbind(cor_1, cor_2)

custom_colors <- c("#9875bd", "#ecb613")

d <- correlation_coeffs %>%
  ggplot(aes(x = group, y = cor, fill = group)) +
  geom_violin(alpha = 0.5) + geom_sina() +
  scale_colour_manual(values = custom_colors, aesthetics = c("colour", "fill")) +
  #ggtitle(label = "CNN Class Predictions Compared to \n Effective Number of Classes in Alignment") +
  theme_cowplot() + 
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5), legend.position = "none") +
  ylab("Correlation Coefficients") +
  xlab("")

#coord_cartesian(xlim = c(0.0, 0.5), ylim = c(0, 3)) +
#scale_x_continuous(breaks = seq(0, 0.6, by = 0.1)) +
#scale_y_continuous(breaks = seq(0, 3, by = 0.5))

#=======================================================================
# making and saving final plot
#=======================================================================
p1 <- plot_grid(a, a2, b, b2, c, c2, nrow = 3, align="hv", axis = "b", labels = c('A', 'B', 'C', 'D', 'E', 'F'))
p2 <- plot_grid(d, labels = "G")
ggsave(filename = "figure_4b.png", plot = p2, width = 8, height = 4)

figure_4 <- plot_grid(p1, p2, nrow = 2, align = "hv", axis = "b", rel_heights = c(3, 1))

ggsave(filename = "figure_4.png", plot = figure_4, width = 10, height = 12)



