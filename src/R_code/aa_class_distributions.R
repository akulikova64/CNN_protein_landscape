# comparing amino acid class distributions

library(tidyverse)
library(yardstick)
library(ggpubr)
library(cowplot)
library(sqldf)
library(sinaplot)
library(dplyr)
library(ggforce)


cnn_data <- read.csv(file="./cnn_wt_max_freq.csv", header=TRUE, sep=",")
natural_data <- read.csv(file = "./natural_max_freq.csv", header=TRUE, sep=",")

joined_data <- rbind(x = cnn_data, y = natural_data)

#==============================================================================
#Basic stats (% match by aa class across all genes)
#===============================================================================

# check if predicted aa_class is the one found in the wt structure
matches_1 <- joined_data %>%
  select(c(gene, position, group, aa_class, class_freq)) %>%
  pivot_wider(names_from = group, values_from = c(aa_class, class_freq)) %>%
  mutate(match = aa_class_predicted == aa_class_wt) %>%
  select(gene, position, match)

# selecting the data entries where the predicted amino acid matches the 
summary_stats_1 <- matches_1 %>%
  group_by(gene, match) %>%
  summarise(count = n()) %>%
  mutate(freq = count / sum(count)) %>%
  filter(match == TRUE)

summary_stats_1 %>%
  ggplot(aes(x = "", y = freq)) +
  geom_violin(fill = "grey75") + geom_sina() +
  ggtitle(label = "CNN Accuracy by Protein") +
  theme_cowplot() + 
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) +
  xlab("") +
  ylab("Percent Accuracy \n (predicted class = wt class)")   

summary_stats_1 %>%
  summarise(mean = mean(freq))
# mean = 0.713

# check if predicted aa is the consensus
matches_2 <- joined_data %>%
  select(c(gene, position, group, aa_class, class_freq)) %>%
  pivot_wider(names_from = group, values_from = c(aa_class,class_freq)) %>%
  mutate(match = aa_class_predicted == aa_class_natural_max) %>%
  select(gene, position, match)

library(dplyr)
summary_stats_2 <- matches_2 %>%
  group_by(gene, match) %>%
  summarise(count = n()) %>%
  mutate(freq = count / sum(count)) %>%
  filter(match == TRUE)

summary_stats_2 %>%
  ggplot(aes(x = "", y = freq)) +
  geom_violin(fill = "grey75") + geom_sina() +
  ggtitle(label = "CNN Class Predictions Compared to \n Alignment Consensus Class") +
  theme_cowplot() + 
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) +
  xlab("") +
  ylab("Percent Accuracy \n (predicted class = consensus class)")   

summary_stats_2 %>%
  summarise(mean = mean(freq))
# mean = 0.331

#==============================================================================
#Comparing class n_eff predicted to class n_eff natural
#===============================================================================

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
n_eff_data_wide <- n_eff_data %>%
  pivot_wider(names_from = group, values_from = n_eff_class)

goi = "2eu8"

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
  stat_cor(label.y = 7, label.x = 3) +
  stat_regline_equation(label.y = 6, label.x = 3) +
  theme_cowplot() +
  theme(plot.title = element_text(hjust=0.5), plot.subtitle = element_text(hjust=0.5)) +
  ylab("N-eff Class Predicted") +
  xlab("N-eff Class Natural")

#====================================================================================
# making violin plot of correlations for all proteins ordered from highest to lowest
#====================================================================================
library(dplyr)
cor <- n_eff_data_wide %>%
  select(gene, natural, predicted) %>%
  na.omit() %>%
  group_by(gene) %>%
  summarise(cor = cor(natural, predicted))

cor_test <- n_eff_data_wide %>%
  filter(gene == goi) %>%
  
  
cor.test(n_eff_data_wide$natural, n_eff_data_wide$predicted, conf.level=conf.level)


cor %>%
  ggplot(aes(x = "", y = cor)) +
  geom_violin(fill = "grey75") + geom_sina() +
  ggtitle(label = "N-eff class natural ~ N-effclass predicted") +
  theme_cowplot() +
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) +
  xlab("") +
  ylab("correlation coefficients") 

# paired sample t-test
ks.test(n_eff_data_wide$predicted, n_eff_data_wide$natural)

#--------------------------------------------------------------------------------------
# 1) Probability of wt class (CNN) ~ Frequency of wt class (seq alignment) 
#--------------------------------------------------------------------------------------

# rerun the joined_data at the top
all_data <- joined_data %>%
  select(-c(aa, freq)) %>%
  pivot_wider(names_from = group, values_from = c(aa_class, class_freq))

goi = "2eu8"

#finding the linear regression
linear_reg <- all_data %>%
  filter(gene == goi) %>%
  select(class_freq_wt, class_freq_natural_wt) %>%
  lm(formula = class_freq_wt ~ class_freq_natural_wt)
summary(linear_reg)

slope <- linear_reg$coefficients[2]
intercept <- linear_reg$coefficients[1]

all_data %>%
  filter(gene == goi) %>%
  ggplot(aes(x = class_freq_wt, y = class_freq_natural_wt)) +
  geom_point() +
  ggtitle(label = goi) +
  geom_abline(slope = slope, intercept = intercept, color = "grey65" ) +
  stat_cor(label.y = 0.95, label.x = -0.1) +
  stat_regline_equation(label.y = .85, label.x = -0.1) +
  theme_cowplot() +
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) +
  ylab("wt class frequency (alignment)") +
  xlab("wt class probability (cnn)") 

library(dplyr)
cor <- all_data %>%
  select(gene, class_freq_wt, class_freq_natural_wt) %>%
  na.omit() %>%
  group_by(gene) %>%
  summarize(cor = cor(class_freq_wt, class_freq_natural_wt))

cor %>%
  ggplot(aes(x = "", y = cor)) +
  geom_violin(fill = "grey75") + geom_sina() +
  ggtitle(label = "wt class prob (CNN) ~ wt class freq (alignment)") +
  theme_cowplot() +
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) +
  xlab("") +
  ylab("correlation coefficients") 


