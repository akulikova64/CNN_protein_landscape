library(tidyverse)
library(cowplot)
library(ggforce)

box_size = "20"
pred_wt_aa = 0.592
pred_wt_class = 0.713

# loading data
cnn_data <- read.csv(file = paste0("./data/PSICOV_box_",box_size,"/output/cnn_wt_max_freq.csv"), header=TRUE, sep=",")
natural_data_all <- read.csv(file = paste0("./data/PSICOV_box_",box_size,"/output/natural_max_freq_files/natural_max_freq_all.csv"), header=TRUE, sep=",")
natural_data_20 <- read.csv(file = paste0("./data/PSICOV_box_",box_size,"/output/natural_max_freq_files/natural_max_freq_20.csv"), header=TRUE, sep=",")
natural_data_20 <- read.csv(file = paste0("./data/PSICOV_box_",box_size,"/output/natural_max_freq_files/natural_max_freq_20.csv"), header=TRUE, sep=",")
natural_data_40 <- read.csv(file = paste0("./data/PSICOV_box_",box_size,"/output/natural_max_freq_files/natural_max_freq_40.csv"), header=TRUE, sep=",")
natural_data_60 <- read.csv(file = paste0("./data/PSICOV_box_",box_size,"/output/natural_max_freq_files/natural_max_freq_60.csv"), header=TRUE, sep=",")
natural_data_80 <- read.csv(file = paste0("./data/PSICOV_box_",box_size,"/output/natural_max_freq_files/natural_max_freq_80.csv"), header=TRUE, sep=",")
natural_data_100 <- read.csv(file = paste0("./data/PSICOV_box_",box_size,"/output/natural_max_freq_files/natural_max_freq_100.csv"), header=TRUE, sep=",")


joined_data_20 <- rbind(x = cnn_data, y = natural_data_20)
joined_data_40 <- rbind(x = cnn_data, y = natural_data_40)
joined_data_60 <- rbind(x = cnn_data, y = natural_data_60)
joined_data_80 <- rbind(x = cnn_data, y = natural_data_80)
joined_data_100 <- rbind(x = cnn_data, y = natural_data_100)
joined_data_all <- rbind(x = cnn_data, y = natural_data_all)

joined_data_all <- joined_data_all %>%
  mutate(perc_sim = "(0-100]")

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
  filter(!gene %in% c('1dbx', '1fvg', '1k7j', '1kq6', '1kw4', '1lpy', '1ne2', '1ny1', '1pko', '1rw1', '1vhu', '1w0h', '1wkc'))

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

#all_data2 <- rbind(sim20, sim40, sim60, sim80, sim100) #optional dataset

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
  pivot_longer(c(freq_predict_cons, freq_wt_cons), values_to = "freq", names_to = "condition")

custom_fills <- c("#8c7b9d", "#d2a92d")
custom_colors <- c("#655775", "#8a7228")

data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}

# plot 5
figure_5a <- stats_for_plot2 %>%
  filter(condition == "freq_predict_cons") %>%
  ggplot(aes(y = freq, x = fct_rev(perc_sim))) +
  geom_violin(fill = "#8c7b9d", color = "#655775", alpha = 0.5) + 
  geom_hline(yintercept = pred_wt_aa, linetype = "dashed", color = "#8a0f0f", alpha = 0.8, size = 0.85) +
  #geom_sina(size = 0.2) +
  stat_summary(fun.data=data_summary) +
  theme_cowplot(16) + 
  theme(plot.title = element_text(hjust = 0.5), 
        plot.subtitle = element_text(hjust = 0.5),
        axis.text = element_text(color = "black", size = 16),
        panel.grid.major.y = element_line(color = "grey92", size=0.5)
        ) +
  #labs(title = "CNN Predictions Compared to Alignment Consensus", 
       #subtitle = "Amino Acid Predictions") +
  scale_x_discrete(
    name = "Percent sequence similarity to wild type"
  ) +
  scale_y_continuous(
    name = "Accuracy",
    limits = c(0, 1.0),
    breaks = seq(from = 0, to = 1.0, by = 0.1),
    expand = c(0, 0))
  
figure_5a

ggsave(filename = paste0("./analysis/figures/figure_5a_box_",box_size,".png"), plot = figure_5a, width = 8, height = 4)

#within aa class predictions plot.

class_match <- all_data_wider %>%
  mutate(match_predict_cons = aa_class_predicted == aa_class_natural_max,
         match_wt_cons = aa_class_wt == aa_class_natural_max)

stats_for_class_plot <- class_match %>%
  group_by(gene, perc_sim) %>%
  summarise(freq_predict_cons = sum(match_predict_cons, na.rm = TRUE)/sum(!is.na(match_predict_cons)),
            freq_wt_cons = sum(match_wt_cons, na.rm = TRUE)/sum(!is.na(match_wt_cons)))

stats_for_class_plot2 <- stats_for_class_plot %>%
  pivot_longer(c(freq_predict_cons, freq_wt_cons), values_to = "freq", names_to = "condition")

# plot 5b
figure_5b <- stats_for_class_plot2 %>%
  filter(condition == "freq_predict_cons") %>%
  ggplot(aes(y = freq, x = fct_rev(perc_sim))) +
  geom_violin(fill = "#d2a92d", color = "#8a7228", alpha = 0.5) + 
  geom_hline(yintercept = pred_wt_class, linetype = "dashed", color = "#8a0f0f", alpha = 0.8, size = 0.85) +
  #geom_sina(size = 0.2) +
  stat_summary(fun.data=data_summary) +
  theme_cowplot(16) + 
  theme(plot.title = element_text(hjust = 0.5), 
        plot.subtitle = element_text(hjust = 0.5),
        axis.text = element_text(color = "black", size = 16),
        panel.grid.major.y = element_line(color = "grey92", size=0.5)) +
  #labs(title = "CNN Predictions Compared to Alignment Consensus", 
       #subtitle = "Within Class Predictions") +
  scale_x_discrete(
    name = "Percent sequence similarity to wild type") +
  scale_y_continuous(
    name = "Accuracy",
    limits = c(0, 1.0),
    breaks = seq(from = 0, to = 1.0, by = 0.1),
    expand = c(0, 0))

figure_5b

ggsave(filename = paste0("./analysis/figures/figure_5b_",box_size,".png"), plot = figure_5b, width = 8, height = 4)

figure_5 <- plot_grid(figure_5a, figure_5b, nrow = 2, align="h", labels = c('a', 'b'))
ggsave(filename = paste0("./analysis/figures/figure_5_",box_size,".png"), plot = figure_5, width = 8, height = 8)

# #testing...
# 
# a <- tibble(x = 1:4, y = c(TRUE, FALSE, TRUE, TRUE))
# 
# a %>%
#   group_by(y) %>%
#   summarise(freq = sum(y)/sum(!is.na(y)))
# 
# #it works ... cool



