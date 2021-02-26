# continuation of figure 2. 
# comparing neff predicted vs. neff natural accors different alignments similarities. 
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

cnn_var <- read.csv(file="./stats_cnn.csv", header=TRUE, sep=",")
cnn_var2 <- cnn_var %>%
  select(position, gene, n_eff, n_eff_class) %>%
  mutate(group = "predicted")

#0-20
natural_var_20 <- natural_var_20 %>%
  select(position, gene, n_eff, n_eff_class) %>%
  mutate(group = "natural")

  
joined_20 <- rbind(natural_var_20, cnn_var2) %>%
  mutate(perc_sim = "(0-20%]") 

#20-40
natural_var_40 <- natural_var_40 %>%
  select(position, gene, n_eff, n_eff_class) %>%
  mutate(group = "natural")

joined_40 <- rbind(natural_var_40, cnn_var2) %>%
  mutate(perc_sim = "(20-40%]") 

#40-60
natural_var_60 <- natural_var_60 %>%
  select(position, gene, n_eff, n_eff_class) %>%
  mutate(group = "natural") 

joined_60 <- rbind(natural_var_60, cnn_var2) %>%
  mutate(perc_sim = "(40-60%]")  

#60-80
natural_var_80 <- natural_var_80 %>%
  select(position, gene, n_eff, n_eff_class) %>%
  mutate(group = "natural") 

joined_80 <- rbind(natural_var_80, cnn_var2) %>%
  mutate(perc_sim = "(60-80%]") 

#80-100
natural_var_100 <- natural_var_100 %>%
  select(position, gene, n_eff, n_eff_class) %>%
  mutate(group = "natural") 

joined_100 <- rbind(natural_var_100, cnn_var2) %>%
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

data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}

a <- cor %>%
  ggplot(aes(x = perc_sim, y = cor)) +
  geom_violin(fill = "#9875bd", alpha = 0.5) + 
  #geom_sina() +
  stat_summary(fun.data=data_summary) +
  labs(title = "Comparing Predicted Neff to Natural Neff", 
       subtitle = "Amino Acid Predictions") +
  theme_cowplot() + 
  theme(plot.title = element_text(hjust = 0.5), 
        plot.subtitle = element_text(hjust = 0.5),
        panel.grid.major.y = element_line(color = "grey92", size=0.5)
  ) +
  scale_x_discrete(
    name = "Percent Sequence Similarity of Alignment"
  ) +
  scale_y_continuous(
    name = "Correlation Coefficients",
    limits = c(-0.4, 0.6),
    breaks = seq(from = -0.4, to = 0.6, by = 0.1),
    expand = c(0, 0)) 

#Now making a plot for classes

all_joined_wide2 <- all_joined %>%
  select(-n_eff) %>%
  pivot_wider(names_from = group, values_from = n_eff_class)

cor_2 <- all_joined_wide2 %>%
  na.omit() %>%
  group_by(gene, perc_sim) %>%
  summarise(cor = cor(natural, predicted)) 

b <- cor_2 %>%
  ggplot(aes(x = perc_sim, y = cor)) +
  geom_violin(fill = "#ecb613", alpha = 0.5) + 
  #geom_sina() +
  stat_summary(fun.data=data_summary) +
  labs(title = "Comparing Predicted Neff to Natural Neff", 
       subtitle = "Within Class Predictions") +
  theme_cowplot() + 
  theme(plot.title = element_text(hjust = 0.5), 
        plot.subtitle = element_text(hjust = 0.5),
        panel.grid.major.y = element_line(color = "grey92", size=0.5)
  ) +
  scale_x_discrete(
    name = "Percent Sequence Similarity of Alignment"
  ) +
  scale_y_continuous(
    name = "Correlation Coefficients",
    limits = c(-0.4, 0.6),
    breaks = seq(from = -0.4, to = 0.6, by = 0.1),
    expand = c(0, 0)) 


figure_1 <- plot_grid(a, b, nrow = 2, align="v", labels = c('A', 'B'))

ggsave(filename = "../../analysis/figures/figure_6a.png", plot = a, width = 10, height = 4)
ggsave(filename = "../../analysis/figures/figure_6b.png", plot = b, width = 10, height = 4)






