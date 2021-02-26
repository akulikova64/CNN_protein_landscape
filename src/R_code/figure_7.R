#this code is a continuation of figure 4
#correlation coefficients of wt prob predicted by CNN as a function of 
#neff in the alignment.

library(tidyverse)
library(cowplot)
library(sinaplot)
library(ggforce)

# reading csv files
natural_var <- read.csv(file="./stats_align_all.csv", header=TRUE, sep=",")
natural_var_20 <- read.csv(file = "./stats_align_files/stats_align_20.csv", header=TRUE, sep=",")
natural_var_40 <- read.csv(file = "./stats_align_files/stats_align_40.csv", header=TRUE, sep=",")
natural_var_60 <- read.csv(file = "./stats_align_files/stats_align_60.csv", header=TRUE, sep=",")
natural_var_80 <- read.csv(file = "./stats_align_files/stats_align_80.csv", header=TRUE, sep=",")
natural_var_100 <- read.csv(file = "./stats_align_files/stats_align_100.csv", header=TRUE, sep=",")

cnn_data <- read.csv(file = "./cnn_wt_max_freq.csv", header=TRUE, sep=",")
cnn_data <- cnn_data %>%
  filter(group == "wt") %>%
  select(position, gene, freq) %>%
  mutate(group = "predicted") %>%
  mutate(n_eff_class = NA)

#0-20
natural_var_20 <- natural_var_20 %>%
  select(position, gene, n_eff, n_eff_class) %>%
  mutate(group = "natural")


joined_20 <- rbind(natural_var_20, cnn_data) %>%
  mutate(perc_sim = "(0-20%]") 

#20-40
natural_var_40 <- natural_var_40 %>%
  select(position, gene, n_eff, n_eff_class) %>%
  mutate(group = "natural")

joined_40 <- rbind(natural_var_40, cnn_data) %>%
  mutate(perc_sim = "(20-40%]") 

#40-60
natural_var_60 <- natural_var_60 %>%
  select(position, gene, n_eff, n_eff_class) %>%
  mutate(group = "natural") 

joined_60 <- rbind(natural_var_60, cnn_data) %>%
  mutate(perc_sim = "(40-60%]")  

#60-80
natural_var_80 <- natural_var_80 %>%
  select(position, gene, n_eff, n_eff_class) %>%
  mutate(group = "natural") 

joined_80 <- rbind(natural_var_80, cnn_data) %>%
  mutate(perc_sim = "(60-80%]") 

#80-100
natural_var_100 <- natural_var_100 %>%
  select(position, gene, n_eff, n_eff_class) %>%
  mutate(group = "natural") 

joined_100 <- rbind(natural_var_100, cnn_data) %>%
  mutate(perc_sim = "(80-100%]") 

#0-100
natural_var_all <-natural_var %>%
  select(position, gene, n_eff, n_eff_class) %>%
  mutate(group = "natural")

joined_0_to_100 <- rbind(natural_var_all, cnn_var2) %>%
  mutate(perc_sim = "(0-100%]") 



all_joined <- rbind(joined_20, joined_40, joined_60, joined_80, joined_100, joined_0_to_100)

all_joined_wide <- all_joined %>%
  select(-n_eff_class) %>%
  pivot_wider(names_from = group, values_from = n_eff)

cor <- all_joined_wide %>%
  na.omit() %>%
  group_by(gene, perc_sim) %>%
  summarise(cor = cor(natural, predicted))