library(tidyverse)
library(yardstick)
library(ggpubr)
library(cowplot)
library(sqldf)
library(sinaplot)
library(dplyr)
library(ggforce)

# code for figure 2

#=======================================================================
# Comparing n-eff predicted to n-eff natural (scatter plots)
#=======================================================================

# reading csv files
natural_var <- read.csv(file="./stats_align_all.csv", header=TRUE, sep=",")
cnn_var <- read.csv(file="./stats_cnn.csv", header=TRUE, sep=",")

# joining the two data frames.

natural_var <- natural_var %>%
  nest(q_natural = c(q_H:q_C))

cnn_var <- cnn_var %>%
  nest(q_cnn = c(q_H:q_C))

joined_data <- inner_join(x = natural_var, y = cnn_var, by = c('position', 'gene'))

# comparing n_eff between cnn and natural data.
n_eff_data_1 <- joined_data %>%
  select(n_eff.x, gene, position)
names(n_eff_data_1) <- c('n_eff', 'gene', 'position')
n_eff_data_1$group <- 'natural'

n_eff_data_2 <- joined_data %>%
  select(n_eff.y, gene, position)
names(n_eff_data_2) <- c('n_eff', 'gene', 'position')
n_eff_data_2$group <- 'predicted'

n_eff_data <- rbind(n_eff_data_1, n_eff_data_2)


# making df wider for plot
n_eff_data_wide <- n_eff_data %>%
  pivot_wider(names_from = group, values_from = n_eff)

# finding all the corelations of all proteins
library(dplyr)
cor_1 <- n_eff_data_wide %>%
  select(gene, natural, predicted) %>%
  na.omit() %>%
  group_by(gene) %>%
  summarise(cor = cor(natural, predicted)) %>%
  mutate(group = "by aa")

goi_1 = "1w7w"
goi_2 = "1znn"
goi_3 = "1ci0"

#finding the linear regression
linear_reg_1 <- n_eff_data_wide %>%
  filter(gene == goi_1) %>%
  select(predicted, natural) %>%
  lm(formula = predicted ~ natural)
summary(linear_reg)

slope_1 <- linear_reg_1$coefficients[2]
intercept_1 <- linear_reg_1$coefficients[1]

purple <- "#9875bd"
mustard <- "#ecb613"

# scatter plot
a <- n_eff_data_wide %>%
  filter(gene == goi_1) %>%
  ggplot(aes(x = natural, y = predicted)) +
  geom_point(color = purple) +
  ggtitle(label = goi) +
  geom_abline(slope = slope_1, intercept = intercept_1, color = "grey65" ) +
  stat_cor(label.y = 12, label.x = 8) +
  stat_regline_equation(label.y = 10.5, label.x = 8) +
  theme_cowplot() +
  scale_colour_manual(values = purple, aesthetics = c("colour", "fill")) +
  theme(plot.title = element_text(hjust=0.5), plot.subtitle = element_text(hjust=0.5)) +
  ylab("N-eff Predicted") +
  xlab("N-eff Natural") +
  coord_cartesian(xlim = c(1, 12.5), ylim = c(1, 12.5)) +
  scale_x_continuous(breaks = seq(2, 12, by = 2)) +
  scale_y_continuous(breaks = seq(2, 12, by = 2))

# second gene:
linear_reg_2 <- n_eff_data_wide %>%
  filter(gene == goi_2) %>%
  select(predicted, natural) %>%
  lm(formula = predicted ~ natural)
summary(linear_reg_2)

slope_2 <- linear_reg_2$coefficients[2]
intercept_2 <- linear_reg_2$coefficients[1]

# scatter plot
b <- n_eff_data_wide %>%
  filter(gene == goi_2) %>%
  ggplot(aes(x = natural, y = predicted)) +
  geom_point(color = purple) +
  ggtitle(label = goi_2) +
  geom_abline(slope = slope_2, intercept = intercept_2, color = "grey65" ) +
  stat_cor(label.y = 12, label.x = 8) +
  stat_regline_equation(label.y = 10.5, label.x = 8) +
  theme_cowplot() +
  scale_colour_manual(values = purple, aesthetics = c("colour", "fill")) +
  theme(plot.title = element_text(hjust=0.5), plot.subtitle = element_text(hjust=0.5)) +
  ylab("N-eff Predicted") +
  xlab("N-eff Natural")+
  coord_cartesian(xlim = c(1, 12.5), ylim = c(1, 12.5)) +
  scale_x_continuous(breaks = seq(2, 12, by = 2)) +
  scale_y_continuous(breaks = seq(2, 12, by = 2))

# third gene:
linear_reg_3 <- n_eff_data_wide %>%
  filter(gene == goi_3) %>%
  select(predicted, natural) %>%
  lm(formula = predicted ~ natural)
summary(linear_reg_3)

slope_3 <- linear_reg_3$coefficients[2]
intercept_3 <- linear_reg_3$coefficients[1]

# scatter plot
c <- n_eff_data_wide %>%
  filter(gene == goi_3) %>%
  ggplot(aes(x = natural, y = predicted)) +
  geom_point(color = purple) +
  ggtitle(label = goi_3) +
  geom_abline(slope = slope_3, intercept = intercept_3, color = "grey65" ) +
  stat_cor(label.y = 12, label.x = 8) +
  stat_regline_equation(label.y = 10.5, label.x = 8) +
  theme_cowplot() +
  scale_colour_manual(values = purple, aesthetics = c("colour", "fill")) +
  theme(plot.title = element_text(hjust=0.5), plot.subtitle = element_text(hjust=0.5)) +
  ylab("N-eff Predicted") +
  xlab("N-eff Natural") +
  coord_cartesian(xlim = c(1, 12.5), ylim = c(1, 12.5)) +
  scale_x_continuous(breaks = seq(2, 12, by = 2)) +
  scale_y_continuous(breaks = seq(2, 12, by = 2))

#============================================================================================
# comparing the effective number of predicted classes to effective number of natural classes
#============================================================================================
#restart r

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
cnn_var <- read.csv(file="./stats_cnn.csv", header=TRUE, sep=",")

# joining the two data frames.

natural_var <- natural_var %>%
  select(-c(q_H:n_eff)) %>%
  nest(q_natural = c(q_aliphatic:q_proline))

cnn_var <- cnn_var %>%
  select(-c(q_H:n_eff)) %>%
  nest(q_cnn = c(q_aliphatic:q_proline))

joined_data <- inner_join(x = natural_var, y = cnn_var, by = c('position', 'gene'))

# comparing n_eff between cnn and natural data.
n_eff_data_1 <- joined_data %>%
  select(n_eff_class.x, gene, position)
names(n_eff_data_1) <- c('n_eff_class', 'gene', 'position')
n_eff_data_1$group <- 'natural'

n_eff_data_2 <- joined_data %>%
  select(n_eff_class.y, gene, position)
names(n_eff_data_2) <- c('n_eff_class', 'gene', 'position')
n_eff_data_2$group <- 'predicted'

n_eff_data <- rbind(n_eff_data_1, n_eff_data_2)

# making df wider for plot
n_eff_data_wide_2 <- n_eff_data %>%
  pivot_wider(names_from = group, values_from = n_eff_class)

library(dplyr)
cor_2 <- n_eff_data_wide_2 %>%
  select(gene, natural, predicted) %>%
  na.omit() %>%
  group_by(gene) %>%
  summarise(cor = cor(natural, predicted)) %>%
  mutate(group = "by class")


#finding the linear regression
linear_reg_4 <- n_eff_data_wide %>%
  filter(gene == goi_1) %>%
  select(predicted, natural) %>%
  lm(formula = predicted ~ natural)
summary(linear_reg_4)

slope_4 <- linear_reg_4$coefficients[2]
intercept_4 <- linear_reg_4$coefficients[1]

# scatter plot
a2 <- n_eff_data_wide_2 %>%
  filter(gene == goi_1) %>%
  ggplot(aes(x = natural, y = predicted)) +
  geom_point(color = mustard) +
  ggtitle(label = goi_1) +
  geom_abline(slope = slope_4, intercept = intercept_4, color = "grey65" ) +
  stat_cor(label.y = 12, label.x = 8) +
  stat_regline_equation(label.y = 10.5, label.x = 8) +
  theme_cowplot()+
  scale_colour_manual(values = mustard, aesthetics = c("colour", "fill")) +
  theme(plot.title = element_text(hjust=0.5), plot.subtitle = element_text(hjust=0.5)) +
  ylab("N-eff Class Predicted") +
  xlab("N-eff Class Natural") +
  coord_cartesian(xlim = c(1, 12), ylim = c(1, 12)) +
  scale_x_continuous(breaks = seq(2, 12, by = 2)) +
  scale_y_continuous(breaks = seq(2, 12, by = 2))

#second plot: finding the linear regression
linear_reg_5 <- n_eff_data_wide %>%
  filter(gene == goi_2) %>%
  select(predicted, natural) %>%
  lm(formula = predicted ~ natural)
summary(linear_reg_5)

slope_5 <- linear_reg_5$coefficients[2]
intercept_5 <- linear_reg_5$coefficients[1]

# scatter plot
b2 <- n_eff_data_wide_2 %>%
  filter(gene == goi_2) %>%
  ggplot(aes(x = natural, y = predicted)) +
  geom_point(color = mustard) +
  ggtitle(label = goi_2) +
  geom_abline(slope = slope_5, intercept = intercept_5, color = "grey65" ) +
  stat_cor(label.y = 12, label.x = 8) +
  stat_regline_equation(label.y = 10.5, label.x = 8) +
  theme_cowplot()+
  scale_colour_manual(values = mustard, aesthetics = c("colour", "fill")) +
  theme(plot.title = element_text(hjust=0.5), plot.subtitle = element_text(hjust=0.5)) +
  ylab("N-eff Class Predicted") +
  xlab("N-eff Class Natural")+
  coord_cartesian(xlim = c(1, 12), ylim = c(1, 12)) +
  scale_x_continuous(breaks = seq(2, 12, by = 2)) +
  scale_y_continuous(breaks = seq(2, 12, by = 2))

#third plot: finding the linear regression
linear_reg_6 <- n_eff_data_wide %>%
  filter(gene == goi_3) %>%
  select(predicted, natural) %>%
  lm(formula = predicted ~ natural)
summary(linear_reg_6)

slope_6 <- linear_reg_6$coefficients[2]
intercept_6 <- linear_reg_6$coefficients[1]

# scatter plot
c2 <- n_eff_data_wide_2 %>%
  filter(gene == goi_3) %>%
  ggplot(aes(x = natural, y = predicted)) +
  geom_point(color = mustard) +
  ggtitle(label = goi_3) +
  geom_abline(slope = slope_6, intercept = intercept_6, color = "grey65" ) +
  stat_cor(label.y = 12, label.x = 8) +
  stat_regline_equation(label.y = 10.5, label.x = 8) +
  theme_cowplot()+
  scale_colour_manual(values = mustard, aesthetics = c("colour", "fill")) +
  theme(plot.title = element_text(hjust=0.5), plot.subtitle = element_text(hjust=0.5)) +
  ylab("N-eff Class Predicted") +
  xlab("N-eff Class Natural") +
  coord_cartesian(xlim = c(1, 12), ylim = c(1, 12)) +
  scale_x_continuous(breaks = seq(2, 12, by = 2)) +
  scale_y_continuous(breaks = seq(2, 12, by = 2))

#=========================================================================================================
# density plot comparing the distribution of correlation coefficients between the two groups (aa vs class)
#=========================================================================================================

correlation_coeffs <- rbind(cor_1, cor_2)

custom_colors <- c("#9875bd", "#ecb613")

d <- correlation_coeffs %>%
  ggplot(aes(x = cor, fill = group)) +
  geom_density(alpha = 0.5) + 
  scale_colour_manual(values = custom_colors, aesthetics = c("colour", "fill")) +
  #ggtitle(label = "CNN Predictions Compared to \n Alignment Consensus") +
  theme_cowplot() + 
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5), legend.position = "none") +
  ylab("Count") +
  xlab("Correlation Coefficients (R)") +
  coord_cartesian(xlim = c(0.0, 0.5), ylim = c(0, 3)) +
  scale_x_continuous(breaks = seq(0, 0.6, by = 0.1)) +
  scale_y_continuous(breaks = seq(0, 3, by = 0.5))

#=======================================================================
# making and saving final plot
#=======================================================================
p1 <- plot_grid(a, a2, b, b2, c, c2, nrow = 3, align="hv", axis = "b", labels = c('A', 'B', 'C', 'D', 'E', 'F'))
p2 <- plot_grid(d, labels = "G")

figure_2 <- plot_grid(p1, p2, nrow = 2, align = "hv", axis = "b", rel_heights = c(3, 1))

ggsave(filename = "figure_2.png", plot = figure_2, width = 10, height = 12)





