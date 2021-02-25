library(tidyverse)
library(cowplot)
library(ggforce)

# loading data
cnn_data <- read.csv(file = "./cnn_wt_max_freq.csv", header=TRUE, sep=",")
natural_data_20 <- read.csv(file = "./natural_max_freq_files/natural_max_freq_20.csv", header=TRUE, sep=",")
natural_data_40 <- read.csv(file = "./natural_max_freq_files/natural_max_freq_40.csv", header=TRUE, sep=",")
natural_data_60 <- read.csv(file = "./natural_max_freq_files/natural_max_freq_60.csv", header=TRUE, sep=",")
natural_data_80 <- read.csv(file = "./natural_max_freq_files/natural_max_freq_80.csv", header=TRUE, sep=",")
natural_data_100 <- read.csv(file = "./natural_max_freq_files/natural_max_freq_100.csv", header=TRUE, sep=",")


joined_data_20 <- rbind(x = cnn_data, y = natural_data_20)
joined_data_40 <- rbind(x = cnn_data, y = natural_data_40)
joined_data_60 <- rbind(x = cnn_data, y = natural_data_60)
joined_data_80 <- rbind(x = cnn_data, y = natural_data_80)
joined_data_100 <- rbind(x = cnn_data, y = natural_data_100)

joined_data_20 <- joined_data_20 %>%
  mutate(perc_sim = "(0-20%]") 

joined_data_40 <- joined_data_40 %>%
  mutate(perc_sim = "(20-40%]")

joined_data_60 <- joined_data_60 %>%
  mutate(perc_sim = "(40-60%]")

joined_data_80 <- joined_data_80 %>%
  mutate(perc_sim = "(60-80%]")

joined_data_100 <- joined_data_100 %>%
  mutate(perc_sim = "(80-100%]")

all_data <- rbind(joined_data_20, joined_data_40, joined_data_60, joined_data_80, joined_data_100)

all_data_trimmed <- all_data %>%
  filter(gene != c('1dbx', '1fvg', '1k7j', '1kq6', '1kw4', '1lpy', '1ne2', '1ny1', '1pko', '1rw1', '1vhu', '1w0h', '1wkc'))

#joining just the natural data:

sim20 <- natural_data_20 %>%
  mutate(perc_sim = "(0-20%]") 

sim40 <- natural_data_20 %>%
  mutate(perc_sim = "(20-40%]") 

sim60 <- natural_data_20 %>%
  mutate(perc_sim = "(40-60%]") 

sim80 <- natural_data_20 %>%
  mutate(perc_sim = "(60-80%]") 

sim100 <- natural_data_20 %>%
  mutate(perc_sim = "(80-100%]") 

all_data2 <- rbind(sim20, sim40, sim60, sim80, sim100) #optional dataset

# checking for when predicted aa matches the consensus aa in alignment
all_data_wider <- all_data_trimmed %>%
  pivot_wider(names_from = group, values_from = c(aa, freq, aa_class, class_freq))

match_consensus <- all_data_wider %>%
  mutate(match_predict_cons = aa_predicted == aa_natural_max,
         match_wt_cons = aa_wt == aa_natural_max)

stats_for_plot <- match_consensus %>%
  group_by(gene, perc_sim) %>%
  summarise(freq_predict_cons = sum(match_predict_cons, na.rm = TRUE)/sum(!is.na(match_predict_cons)),
            freq_wt_cons = sum(match_wt_cons, na.rm = TRUE)/sum(!is.na(match_wt_cons)))

stats_for_plot2 <- stats_for_plot %>%
  pivot_longer(c(freq_predict_cons, freq_wt_cons), values_to = "freq", names_to = "type")

stats_for_plot2 %>%
  ggplot(aes(y = freq, x = perc_sim, fill = type)) +
  geom_violin(alpha = 0.5) + 
  geom_sina(size = 0.2) +
  theme_cowplot() + NULL
  facet_wrap(vars(type), ncol = 1) + NULL
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5), legend.position = "none") +
  ggtitle("CNN Predictions Compared to Alignment Consensus \n Alignment Consensus") +
  scale_x_discrete(
    name = "Percent Sequence Similarity of Alignment"
  ) +
  scale_y_continuous(
    name = "Accuracy",
    limits = c(0, 0.6),
    breaks = c(0, 0.20, 0.40, 0.60))
  
#lets see the outliers
stats_for_plot2 %>%
  filter(type == "freq_wt_cons")

ggsave(filename = "figure_5.png", plot = figure_5, width = 6, height = 4)







