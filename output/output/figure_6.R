# comparing neff predicted vs. neff natural accors different alignments similarities. 

# reading csv files
natural_var <- read.csv(file="./stats_align_all.csv", header=TRUE, sep=",")
natural_var_20 <- read.csv(file = "./stats_align_files/stats_align_20.csv", header=TRUE, sep=",")
natural_var_40 <- read.csv(file = "./stats_align_files/stats_align_40.csv", header=TRUE, sep=",")
natural_var_60 <- read.csv(file = "./stats_align_files/stats_align_60.csv", header=TRUE, sep=",")
natural_var_80 <- read.csv(file = "./stats_align_files/stats_align_80.csv", header=TRUE, sep=",")
natural_var_100 <- read.csv(file = "./stats_align_files/stats_align_100.csv", header=TRUE, sep=",")

cnn_var <- read.csv(file="./stats_cnn.csv", header=TRUE, sep=",")

natural_var_20 <- natural_var_20 %>%
  mutate(perc_sim = "(0-20%]") 

natural_var_40 <- natural_var_40 %>%
  mutate(perc_sim = "(0-40%]") 

natural_var_60 <- natural_var_60 %>%
  mutate(perc_sim = "(0-60%]") 

natural_var_80 <- natural_var_80 %>%
  mutate(perc_sim = "(0-80%]") 

natural_var_100 <- natural_var_100 %>%
  mutate(perc_sim = "(0-100%]") 

cnn_var <- subset(cnn_var, select = -c(wt_aa))

natural_var_all <- rbind(natural_var_20, natural_var_40, natural_var_60, natural_var_80, natural_var_100)

natural_var_all <- natural_var_all %>%
  nest(q_natural = c(q_H:q_C))

cnn_var <- cnn_var %>%
  nest(q_cnn = c(q_H:q_C))

joined_data <- inner_join(x = natural_var_all, y = cnn_var, by = c('position', 'gene'))










