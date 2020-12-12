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

#finding mean
n_eff_data <- n_eff_data %>%
  group_by(gene, group) %>%
  summarise(mean = mean(n_eff))

# making df wider for plot
n_eff_data_wide <- n_eff_data %>%
  pivot_wider(names_from = group, values_from = n_eff)

goi = "1w7w"

#finding the linear regression
linear_reg <- n_eff_data_wide %>%
  filter(gene == goi) %>%
  select(predicted, natural) %>%
  lm(formula = predicted ~ natural)
  summary(linear_reg)

slope <- linear_reg$coefficients[2]
intercept <- linear_reg$coefficients[1]

# scatter plot
n_eff_data_wide %>%
  filter(gene == goi) %>%
  ggplot(aes(x = natural, y = predicted)) +
  geom_point() +
  ggtitle(label = goi) +
  geom_abline(slope = slope, intercept = intercept, color = "grey65" ) +
  stat_cor(label.y = 12, label.x = 8) +
  stat_regline_equation(label.y = 10.5, label.x = 8) +
  theme_cowplot() +
  theme(plot.title = element_text(hjust=0.5), plot.subtitle = element_text(hjust=0.5)) +
  ylab("N-eff Predicted") +
  xlab("N-eff Natural")


#====================================================================================
# making violin plot of correlations for all proteins ordered from highest to lowest
#====================================================================================
cor <- n_eff_data_wide %>%
  select(gene, natural, predicted) %>%
  na.omit() %>%
  group_by(gene) %>%
  summarise(cor = cor(natural, predicted))


cor %>%
  ggplot(aes(x = "", y = cor)) +
  geom_violin(fill = "grey75") + geom_sina() +
  ggtitle(label = "N-eff natural ~ N-eff predicted") +
  theme_cowplot() +
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) +
  xlab("") +
  ylab("correlation coefficients") 
