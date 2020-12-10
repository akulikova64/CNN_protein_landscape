library(tidyverse)
library(yardstick)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(cowplot)
library(tidyr)

# finding the Kullback-Leibler (KL) divergence
# y - CNN, x - natural
get_D_KL <- function(y, x) {
  sum <- 0
  for (i in 1:20) { 
    if (x[i] != 0 & y[i] != 0) {
      sum = sum + x[i] * log(x[i]/y[i])
    }
  }
  return(sum)
}

# a better function that Claus made:
get_D_KL2 <- function(y, x) {
  sum(ifelse(x == 0 | y == 0, 0, x*log(x/y)))
}

# reading csv files
natural_var <- read.csv(file="./natural_variability_2.csv", header=TRUE, sep=",")
natural_var_first <- read.csv(file="./natural_variability_first.csv", header=TRUE, sep=",")
natural_var_second <- read.csv(file="./natural_variability_second.csv", header=TRUE, sep=",")
cnn_var <- read.csv(file="./CNN_variability.csv", header=TRUE, sep=",")

# joining the two data frames for calculating the D_KL later on.

natural_var <- natural_var %>%
  nest(q_natural = c(q_H:q_C))

cnn_var <- cnn_var %>%
  nest(q_cnn = c(q_H:q_C))

joined_data <- inner_join(x = natural_var, y = cnn_var, by = c('position', 'gene'))

# adding a column D_KL (values calculated in function above)
joined_data <- joined_data %>%
  rowwise() %>%
  mutate(D_KL = get_D_KL2(as.numeric(q_natural[1,]), as.numeric(q_cnn[1,])))

# natural vs. natural
natural_var_first <- natural_var_first %>%
  nest(q_natural_1 = c(q_H:q_C))

natural_var_second <- natural_var_second %>%
  nest(q_natural_2 = c(q_H:q_C))

joined_natural <- inner_join(x = natural_var_first, y = natural_var_second, by = c('position', 'gene'))

joined_natural <- joined_natural %>%
  rowwise() %>%
  mutate(D_KL = get_D_KL2(as.numeric(q_natural_1[1,]), as.numeric(q_natural_2[1,])))

  
#selecting only columns needed for graphing
joined_data2 <- joined_data %>%
   select(D_KL, gene) 
joined_data2$group <- 'predicted'
 
joined_natural2 <- joined_natural %>%
   select(D_KL, gene)
joined_natural2$group <- 'natural'
 
joined_for_graph <- rbind(joined_data2, joined_natural2)

joined_for_graph <- joined_for_graph %>%
  group_by(gene, group) %>%
  summarize(D_KL = mean(D_KL))

#graphing KL divergence (Mean for each gene- 38 points for each box)
joined_for_graph %>%
  ggplot(aes(x = group, y = D_KL)) +
  geom_boxplot() +
  ggtitle(label="KL-Divergence", subtitle = "CNN data vs. natural sequences") +
  theme_cowplot() +
  theme(plot.title = element_text(hjust=0.5), plot.subtitle = element_text(hjust=0.5)) +
  ylab("Mean KL Divergence") +
  xlab("")

#ploting all 38 genes
joined_data %>%
  ggplot(aes(x = "", y = D_KL)) +
  geom_boxplot() +
  ggtitle(label = "KL-Divergence for all 38 proteins") +
  facet_wrap(~gene)+
  #theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) +
  xlab("")
  ylab("Mean KL Divergence") 
  
n_eff_data_1 <- joined_data %>%
  select(n_eff.x, gene, position)
names(n_eff_data_1) <- c('n_eff', 'gene', 'position')
n_eff_data_1$group <- 'natural'

n_eff_data_2 <- joined_data %>%
  select(n_eff.y, gene, position)
names(n_eff_data_2) <- c('n_eff', 'gene', 'position')
n_eff_data_2$group <- 'predicted'

n_eff_data <- rbind(n_eff_data_1, n_eff_data_2)

n_eff_data_3 <- n_eff_data %>%
  group_by(gene, group) %>%
  summarize(n_eff = mean(n_eff))
  
n_eff_data_3 %>%
  ggplot(aes(x = group, y = n_eff)) +
  geom_boxplot() +
  ggtitle(label="N-effective", subtitle = "CNN data vs. natural sequences") +
  theme_cowplot() +
  theme(plot.title = element_text(hjust=0.5), plot.subtitle = element_text(hjust=0.5)) +
  ylab("N-effective") +
  xlab("")

#====================================================================================================
# Replicating results from paper
#====================================================================================================

designed <- read.csv(file="./designed.csv", header=TRUE, sep=",")
evolved <- read.csv(file="./evolved.csv", header=TRUE, sep=",")

designed <- designed %>%
  nest(q_designed = c(q_H:q_C))

evolved <- evolved %>%
  nest(q_evolved = c(q_H:q_C))

joined_evolved <- inner_join(x = natural_var, y = evolved, by = c('position', 'gene'))
joined_designed <- inner_join(x = natural_var, y = designed, by = c('position', 'gene'))

#calculating the KL-Divergence
joined_evolved <- joined_evolved %>%
  rowwise() %>%
  mutate(D_KL = get_D_KL2(as.numeric(q_natural[1,]), as.numeric(q_evolved[1,])))

joined_designed <- joined_designed %>%
  rowwise() %>%
  mutate(D_KL = get_D_KL2(as.numeric(q_natural[1,]), as.numeric(q_designed[1,])))

#selecting only columns needed for graphing
joined_evolved2 <- joined_evolved %>%
  select(D_KL, gene) 
joined_evolved2$group <- 'evolved'

joined_designed2 <- joined_designed %>%
  select(D_KL, gene)
joined_designed2$group <- 'designed'

joined_natural2 <- joined_natural %>%
  select(D_KL, gene)
joined_natural2$group <- 'natural'

joined_for_graph2 <- rbind(joined_evolved2, joined_designed2, joined_natural2)

View(joined_for_graph2)

joined_for_graph2 <- joined_for_graph2 %>%
  group_by(gene, group) %>%
  summarize(D_KL = mean(D_KL))

joined_for_graph2 %>%
  ggplot(aes(x = group, y = D_KL)) +
  geom_boxplot() +
  ggtitle(label="KL-Divergence", subtitle = "designed vs. evolved vs. natural sequences") +
  theme_cowplot() +
  theme(plot.title = element_text(hjust=0.5), plot.subtitle = element_text(hjust=0.5)) +
  ylab("Mean KL Divergence") +
  xlab("")

# adding the cnn predicted results
joined_predicted2 <- joined_data %>%
  select(D_KL, gene) 
joined_predicted2$group <- 'predicted'

joined_evolved2 <- joined_evolved %>%
  select(D_KL, gene) 
joined_evolved2$group <- 'evolved'

joined_designed2 <- joined_designed %>%
  select(D_KL, gene)
joined_designed2$group <- 'designed'

joined_natural2 <- joined_natural %>%
  select(D_KL, gene)
joined_natural2$group <- 'natural'

joined_for_graph3 <- rbind(joined_natural2, joined_designed2, joined_evolved2, joined_predicted2)

joined_for_graph3 <- joined_for_graph3 %>%
  group_by(gene, group) %>%
  summarize(D_KL = mean(D_KL))
  
joined_for_graph3 %>%
  ggplot(aes(x = group, y = D_KL)) +
  geom_boxplot() +
  ggtitle(label="KL-Divergence", subtitle = "designed, evolved, predicted, natural") +
  theme_cowplot() +
  theme(plot.title = element_text(hjust=0.5), plot.subtitle = element_text(hjust=0.5)) +
  ylab("Mean KL Divergence") +
  xlab("")

#=============================================================================================================
# Graphing n_eff and frequencies (predicted vs natural). Predicted vs. consensus
#=============================================================================================================

# x- natural, y- predicted/cnn
joined_data_line <- joined_data %>%
  select(position, gene, natural = n_eff.x, predicted = n_eff.y)

joined_data_line <- joined_data_line %>%
  pivot_longer(cols = c(natural, predicted), names_to = "group", values_to = "n_eff")

joined_data_line %>%
  filter(gene == "1b4t") %>%
  ggplot(aes(x = position, y = n_eff, color = group)) +
  geom_line() +
  geom_point() +
  ggtitle(label = "1b4t") +
  theme_cowplot() +
  theme(plot.title = element_text(hjust=0.5), plot.subtitle = element_text(hjust=0.5)) +
  ylab("n-eff") +
  xlab("position") +
  coord_cartesian(xlim = c(0, 30)) +
  scale_x_continuous(breaks = seq(0, 30, by = 2)) +
  grids(axis = "x", color = "grey92", linetype = "dashed")


#=============================================================================================================
#Tests
#=============================================================================================================


#test function
func <- function(q,p){
  print(p)
  print(q)
  sum <- 0
  for (i in 1:2){
    sum = sum + (q[i] + p[i])
  }
  return(sum)
}
func <- Vectorize(func)


#test
joined_data %>%
  mutate(test = func(c(q_H.y, q_E.y), c(q_H.x, q_E.x)))

# # sampling random rows from half of the natural frequencies in joined_joined data  
# rand_rows <- c(sample.int(8187, 4094))
# 
# # creating two df with half of the natural frequency data
# first_half <- joined_data[-rand_rows,]
# second_half <- joined_data[rand_rows,]
# 
# #cleaning up data frame
# first_half <- first_half[1:24]
# second_half <- second_half[1:24]
# 
# #renaming cols (we don't need a .x because we are going to join two halfs again) 
# names(first_half) <- c('position', 'gene', 'q_H', 'q_E', 'q_D', 'q_R', 'q_K', 'q_S', 'q_T', 'q_N', 'q_Q', 'q_A', 'q_V', 'q_L', 'q_I', 'q_M', 'q_F', 'q_Y', 'q_W', 'q_P', 'q_G', 'q_C', 'entropy', 'n_eff')
# names(second_half) <- c('position', 'gene', 'q_H', 'q_E', 'q_D', 'q_R', 'q_K', 'q_S', 'q_T', 'q_N', 'q_Q', 'q_A', 'q_V', 'q_L', 'q_I', 'q_M', 'q_F', 'q_Y', 'q_W', 'q_P', 'q_G', 'q_C', 'entropy', 'n_eff')
# 
# #adding an arbitrary col with the row number for joining the data later on.
# first_half$row <- seq.int(nrow(first_half))
# second_half$row <- seq.int(nrow(second_half))
# 
# #joining data and adding the D_KL col
# new_natural <- inner_join(first_half, second_half, by = 'row')
# new_natural <- new_natural %>%
#   rowwise() %>%
#   mutate(D_KL = get_D_KL(c(q_H.y, q_E.y, q_D.y,  q_R.y, q_K.y, q_S.y, q_T.y, q_N.y, q_Q.y, q_A.y, q_V.y, q_L.y, q_I.y, q_M.y, q_F.y, q_Y.y, q_W.y, q_P.y, q_G.y, q_C.y), c(q_H.x, q_E.x, q_D.x,  q_R.x, q_K.x, q_S.x, q_T.x, q_N.x, q_Q.x, q_A.x, q_V.x, q_L.x, q_I.x, q_M.x, q_F.x, q_Y.x, q_W.x, q_P.x, q_G.x, q_C.x)))
# 
# # mutate(D_KL = get_D_KL())
# # selecting only columns needed for graphing
# joined_data2 <- joined_data %>%
#   select(D_KL) 
# joined_data2$group <- 'predicted'
# 
# new_natural2 <- new_natural %>%
#   select(D_KL)
# new_natural2$group <- 'natural'
# 
# joined_for_graph <- rbind(joined_data2, new_natural2)
# 
# 
