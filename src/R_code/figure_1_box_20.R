library(tidyverse)
library(cowplot)
library(ggforce)

#======================
#***BOX SIZE 20*****
#======================

# code for figure 1. 
#==============================================================================
# a) Basic stats (% match across all genes)
#===============================================================================

# set working directory to: "Desktop/Natural_var_project/"
# loading data
cnn_data <- read.csv(file = "./data/PSICOV_box_20/output/cnn_wt_max_freq.csv", header=TRUE, sep=",")
natural_data <- read.csv(file = "./data/PSICOV_box_20/output/natural_max_freq_files/natural_max_freq_all.csv", header=TRUE, sep=",")

joined_data <- rbind(x = cnn_data, y = natural_data)

joined_data_trimmed <- joined_data %>%
  filter(!gene %in% c('1dbx', '1eaz', '1fvg', '1k7j', '1kq6', '1kw4', '1lpy', '1ne2', '1ny1', '1pko', '1rw1', '1vhu', '1w0h', '1wkc', '2tps'))

# check if predicted aa is the one found in the wt structure
joined_data_wider <- joined_data_trimmed %>%
  pivot_wider(names_from = group, values_from = c(aa, freq, aa_class, class_freq))

joined_data_wider <- na.omit(joined_data_wider)


match_wt <- joined_data_wider %>%
  mutate(match_predict_wt = aa_predicted == aa_wt) 

#data entries where the predicted amino acid matches the wt
stats_1 <- match_wt %>%
  group_by(gene) %>%
  summarise(freq_predict_wt = sum(match_predict_wt, na.rm = TRUE)/sum(!is.na(match_predict_wt))) %>%
  mutate(group = "predicted aa = wt aa") %>%
  mutate(x_label = "aa \n predictions")

# now, we need to add accuracy within each class (predicted class == wt class)

match_wt_class <- joined_data_wider %>%
  mutate(match_predict_wt_class = aa_class_predicted == aa_class_wt)
  
  

# selecting the data entries where the predicted amino acid class matches the wt aa class.
stats_2 <- match_wt_class %>%
  group_by(gene) %>%
  summarise(freq_predict_wt = sum(match_predict_wt_class, na.rm = TRUE)/sum(!is.na(match_predict_wt_class))) %>%
  mutate(group = "predicted class = wt class") %>%
  mutate(x_label = "class \n predictions")



joined_single_and_class <- rbind(stats_1, stats_2)

custom_fills <- c("#8c7b9d", "#d2a92d")
custom_colors <- c("#655775", "#8a7228")

data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}

plot_20_a <- joined_single_and_class %>%
  ggplot(aes(y = freq_predict_wt, x = x_label, fill = group, color = group)) +
  geom_violin(alpha = 0.6, size = 0.7) +
  stat_summary(fun.data=data_summary, color = "black", alpha = 0.7) +
  #ggtitle(label = "20A Box: Predicting Wild Type") +
  theme_cowplot(16) + 
  theme(plot.title = element_text(hjust = 0, size=16), 
        plot.subtitle = element_text(hjust = 0.5),
        panel.grid.major.y = element_line(color = "grey92", size=0.5),
        legend.position = "none") +
  scale_fill_manual(values = custom_fills) +
  scale_color_manual(values = custom_colors) +
  scale_y_continuous(
    name = "Accuracy",
    limits = c(0.0, 1.0),
    breaks = seq(0, 1.0, by = 0.1),
    expand = c(0, 0)) +
  scale_x_discrete(
    name = "",
    labels = c("Amino Acid \n Prediction", "Class \n Prediction"))

  
plot_20_a

stats_1 %>%
  summarise(mean = mean(freq_predict_wt))
# mean for box_size 20 is = 0.592

stats_2 %>%
  summarise(mean = mean(freq_predict_wt))
# mean for box_size 20 = 0.713

#===================================================================================
# b) check if predicted aa is the most frequent aa in the alignment (consensus)
#===================================================================================
match_cons <- joined_data_wider %>%
  mutate(match_predict_cons = aa_predicted == aa_natural_max)

stats_3 <- match_cons %>%
  group_by(gene) %>%
  summarise(freq_predict_cons = sum(match_predict_cons, na.rm = TRUE)/sum(!is.na(match_predict_cons))) %>%
  mutate(group = "predicted aa = consensus") %>%
  mutate(x_label = "aa \n predictions")

  
match_cons_class <- joined_data_wider %>%
  mutate(match_predict_cons_class = aa_class_predicted == aa_class_natural_max) 
  
stats_4 <- match_cons_class %>%
  group_by(gene) %>%
  summarise(freq_predict_cons = sum(match_predict_cons_class, na.rm = TRUE)/sum(!is.na(match_predict_cons_class))) %>%
  mutate(group = "predicted class = consensus") %>%
  mutate(x_label = "class \n predictions")


joined_consensus <- rbind(stats_3, stats_4)

plot_20_b <- joined_consensus %>%
  ggplot(aes(y = freq_predict_cons, x = x_label, fill = group, color = group)) +
  geom_violin(alpha = 0.6, size = 0.7) + 
  stat_summary(fun.data=data_summary, color = "black", alpha = 0.7) +
  #ggtitle(label = "20A Box: Predicting Consensus") +
  theme_cowplot(16) + 
  theme(plot.title = element_text(hjust = 0, size = 16), 
        plot.subtitle = element_text(hjust = 0.5),
        panel.grid.major.y = element_line(color = "grey92", size=0.5),
        legend.position = "none") +
  scale_fill_manual(values = custom_fills) +
  scale_color_manual(values = custom_colors) +
  scale_y_continuous(
    name = "Accuracy",
    limits = c(0.0, 1.0),
    breaks = seq(0, 1.0, by = 0.1),
    expand = c(0, 0)) +
  scale_x_discrete(
    name = "",
    labels = c("Amino Acid \n Prediction", "Class \n Prediction"))
  
plot_20_b

stats_3 %>%
  summarise(mean = mean(freq_predict_cons))
# mean = 0.383

stats_4 %>%
  summarise(mean = mean(freq_predict_cons))
# mean = 0.552

#==================================================================
# mini bar chart (when pred = cons, what % does not match with wt)
#==================================================================
match_cons <- joined_data_wider %>%
  mutate(match_predict_cons = aa_predicted == aa_natural_max) %>%
  mutate(cons_is_wt = aa_natural_max == aa_wt)

bar_1 <- match_cons %>%
  filter(match_predict_cons == TRUE) %>%
  summarise(freq_cons_not_wt = sum(!cons_is_wt, na.rm = TRUE)/sum(!is.na(cons_is_wt))) %>%
  mutate(group = "consensus aa is not the wt") %>%
  mutate(x_label = "aa \n predictions")

match_cons_class <- joined_data_wider %>%
  mutate(match_predict_cons_class = aa_class_predicted == aa_class_natural_max) %>%
  mutate(cons_class_is_wt = aa_class_natural_max == aa_class_wt)

bar_2 <- match_cons_class %>%
  filter(match_predict_cons_class == TRUE) %>%
  summarise(freq_cons_not_wt = sum(!cons_class_is_wt, na.rm = TRUE)/sum(!is.na(cons_class_is_wt))) %>%
  mutate(group = "consensus class is not the wt") %>%
  mutate(x_label = "class \n predictions")


joined_bars <- rbind(bar_1, bar_2) %>%
  mutate(freq = 100*(1-freq_cons_not_wt)) 

plot_20_c <- joined_bars %>% 
  ggplot(aes(x = group, y = freq, fill = group, color = group)) +
  geom_col(alpha = 0.6) +
  scale_fill_manual(values = c("#8c7b9d", "#d2a92d")) +
  scale_color_manual(values = c("#655775", "#8a7228")) +
  scale_x_discrete(
    name = "",
    labels = c("aa", "Class"))+
  scale_y_continuous(
    name = "Frequency at which \n predicted consensus is wild type",
    limits = c(0, 100),
    breaks = seq(0, 100, by = 10),
    labels = paste0(seq(0, 100, by = 10), "%"),
    expand = c(0,0)) +
  theme_cowplot(16) + 
  theme(
    panel.grid.major.y = element_line(color = "grey92", size=0.5),
    legend.position = "none")

plot_20_c

#==============================================================
# making and saving the plots
#==============================================================

figure_1 <- plot_grid(plot_20_a, plot_20_b, plot_20_c, rel_widths = c(2, 2, 1), nrow = 1, align="h", labels = c('a', 'b', 'c'))

ggsave(filename = "./analysis/figures/figure_1_20_new.png", plot = figure_1, width = 12, height = 5)

