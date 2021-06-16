library(tidyverse)
library(cowplot)
library(ggforce)

#========================================================================================================
# here I will bin positions based on KL-divergence and find the aa distributions within each bin. 
#========================================================================================================

cnn_var <- read.csv(file = paste0("./data/PSICOV_box_20/output/stats_cnn.csv"), header=TRUE, sep=",")
natural_var <- read.csv(file="./output/output_PSICOV/stats_align_all.csv", header=TRUE, sep=",")

cnn_var2 <- cnn_var %>%
  nest(freq = q_A:q_V) %>%
  select(-c(q_aliphatic:n_eff_class)) %>%
  mutate(group = "predicted")

natural_var2 <- natural_var %>%
  nest(freq = q_A:q_V) %>%
  select(-c(q_aliphatic:n_eff_class)) %>%
  mutate(group = "natural")


