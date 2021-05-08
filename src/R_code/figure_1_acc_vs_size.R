library(tidyverse)
library(cowplot)
library(ggforce)

# I will make a figure of accuracy as a function of box size.

cnn_data_12 <- read.csv(file = "./data/PSICOV_box_12/output/cnn_wt_max_freq.csv", header=TRUE, sep=",")
cnn_data_20 <- read.csv(file = "./data/PSICOV_box_20/output/cnn_wt_max_freq.csv", header=TRUE, sep=",")
cnn_data_30 <- read.csv(file = "./data/PSICOV_box_30/output/cnn_wt_max_freq.csv", header=TRUE, sep=",")
cnn_data_40 <- read.csv(file = "./data/PSICOV_box_40/output/cnn_wt_max_freq.csv", header=TRUE, sep=",")

natural_data_12 <- read.csv(file = "./data/PSICOV_box_12/output/natural_max_freq_files/natural_max_freq_all.csv", header=TRUE, sep=",")
natural_data_20 <- read.csv(file = "./data/PSICOV_box_20/output/natural_max_freq_files/natural_max_freq_all.csv", header=TRUE, sep=",")
natural_data_30 <- read.csv(file = "./data/PSICOV_box_30/output/natural_max_freq_files/natural_max_freq_all.csv", header=TRUE, sep=",")
natural_data_40 <- read.csv(file = "./data/PSICOV_box_40/output/natural_max_freq_files/natural_max_freq_all.csv", header=TRUE, sep=",")

# adding a box_size column to cnn data
cnn_data_12 <- cnn_data_12 %>%
  mutate(box_size = "12")

cnn_data_20 <- cnn_data_20 %>%
  mutate(box_size = "20")

cnn_data_30 <- cnn_data_30 %>%
  mutate(box_size = "30")

cnn_data_40 <- cnn_data_40 %>%
  mutate(box_size = "40")

# adding a box_size column to natural data
natural_data_12 <- natural_data_12 %>%
  mutate(box_size = "12")

natural_data_20 <- natural_data_20 %>%
  mutate(box_size = "20")

natural_data_30 <- natural_data_30 %>%
  mutate(box_size = "30")

natural_data_40 <- natural_data_40 %>%
  mutate(box_size = "40")


joined_data <- rbind(x = cnn_data_12, cnn_data_20, cnn_data_30, cnn_data_40, y = natural_data_12, natural_data_20, natural_data_30, natural_data_40)

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
  group_by(gene, box_size) %>%
  summarise(freq_predict_wt = sum(match_predict_wt, na.rm = TRUE)/sum(!is.na(match_predict_wt))) %>%
  mutate(group = "predicted aa = wt aa") %>%
  mutate(x_label = "aa \n predictions")

# now, we need to add accuracy within each class (predicted class == wt class)

match_wt_class <- joined_data_wider %>%
  mutate(match_predict_wt_class = aa_class_predicted == aa_class_wt)



# selecting the data entries where the predicted amino acid class matches the wt aa class.
stats_2 <- match_wt_class %>%
  group_by(gene, box_size) %>%
  summarise(freq_predict_wt = sum(match_predict_wt_class, na.rm = TRUE)/sum(!is.na(match_predict_wt_class))) %>%
  mutate(group = "predicted class = wt class") %>%
  mutate(x_label = "class \n predictions")


#for stats_1: amino acids
custom_fill_1 <- "#8c7b9d"
custom_color_1 <- "#655775"

#for stats_2: amino acid classes
custom_fill_2 <- "#d2a92d"
custom_color_2 <- "#8a7228"

data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}

plot_a <- stats_1 %>%
  ggplot(aes(y = freq_predict_wt, x = box_size)) +
  geom_violin(alpha = 0.5, size = 0.7, fill = custom_fill_1, color = custom_color_1 ) +
  stat_summary(fun.data=data_summary, color = "black", alpha = 0.7) +
  theme_cowplot(12) + 
  theme(plot.title = element_text(hjust = 0, size=12), 
        plot.subtitle = element_text(hjust = 0.5),
        panel.grid.major.y = element_line(color = "grey92", size=0.5),
        legend.position = "none") +
  scale_y_continuous(
    name = "Accuracy",
    limits = c(0.0, 1.0),
    breaks = seq(0, 1.0, by = 0.1),
    expand = c(0, 0)) +
  scale_x_discrete(
    name = "Box Size (Å)")

plot_a

plot_b <- stats_2 %>%
  ggplot(aes(y = freq_predict_wt, x = box_size)) +
  geom_violin(alpha = 0.5, size = 0.7, fill = custom_fill_2, color = custom_color_2 ) +
  stat_summary(fun.data=data_summary, color = "black", alpha = 0.7) +
  theme_cowplot(12) + 
  theme(plot.title = element_text(hjust = 0, size=12), 
        plot.subtitle = element_text(hjust = 0.5),
        panel.grid.major.y = element_line(color = "grey92", size=0.5),
        legend.position = "none") +
  scale_y_continuous(
    name = "Accuracy",
    limits = c(0.0, 1.0),
    breaks = seq(0, 1.0, by = 0.1),
    expand = c(0, 0)) +
  scale_x_discrete(
    name = "Box Size (Å)")

plot_b

figure_1_all_boxes <- plot_grid(plot_a, plot_b, nrow = 2, align="h", labels = c('a', 'b'))
ggsave(filename = paste0("./analysis/figures/figure_1_all_boxes.png"), plot = figure_1_all_boxes, width = 8, height = 8)

#trying a combined plot:

joined_single_and_class <- rbind(stats_1, stats_2)

custom_fills <- c("#8c7b9d", "#d2a92d")
custom_colors <- c("#655775", "#8a7228")

plot_c <- joined_single_and_class %>%
  ggplot(aes(y = freq_predict_wt, x = box_size, fill = group, color = group)) +
  geom_violin(
    alpha = 0.5, 
    size = 0.7,
    position=position_dodge(width = 0.6)) +
  stat_summary(
    fun.data=data_summary, 
    color = "black", 
    alpha = 0.7,
    position=position_dodge(width = 0.6)) +
  theme_cowplot(12) + 
  theme(plot.title = element_text(hjust = 0, size=12), 
        plot.subtitle = element_text(hjust = 0.5),
        panel.grid.major.y = element_line(color = "grey92", size=0.5),
        legend.title = element_blank(),
        legend.position = "top",
        legend.direction = "vertical",
        legend.text = element_text(size = 12)) +
  scale_fill_manual(values = custom_fills) +
  scale_color_manual(values = custom_colors) +
  scale_y_continuous(
    name = "Accuracy",
    limits = c(0.0, 1.0),
    breaks = seq(0, 1.0, by = 0.1),
    expand = c(0, 0)) +
  scale_x_discrete(
    name = "Box Size (Å)")
plot_c


ggsave(filename = paste0("./analysis/figures/figure_1_plot_c.png"), plot = plot_c, width = 8, height = 4)


#===================================================================================
# b) check if predicted aa is the most frequent aa in the alignment (consensus)
#===================================================================================
match_cons <- joined_data_wider %>%
  mutate(match_predict_cons = aa_predicted == aa_natural_max)

stats_3 <- match_cons %>%
  group_by(gene, box_size) %>%
  summarise(freq_predict_cons = sum(match_predict_cons, na.rm = TRUE)/sum(!is.na(match_predict_cons))) %>%
  mutate(group = "predicted aa = consensus") %>%
  mutate(x_label = "aa \n predictions")


match_cons_class <- joined_data_wider %>%
  mutate(match_predict_cons_class = aa_class_predicted == aa_class_natural_max)

stats_4 <- match_cons_class %>%
  group_by(gene, box_size) %>%
  summarise(freq_predict_cons = sum(match_predict_cons_class, na.rm = TRUE)/sum(!is.na(match_predict_cons_class))) %>%
  mutate(group = "predicted class = consensus") %>%
  mutate(x_label = "class \n predictions")


joined_consensus <- rbind(stats_3, stats_4)

plot_d <- joined_consensus %>%
  ggplot(aes(y = freq_predict_cons, x = box_size, fill = group, color = group)) +
  geom_violin(
    alpha = 0.5, 
    size = 0.7,
    position=position_dodge(width = 0.6)) +
  stat_summary(
    fun.data=data_summary, 
    color = "black", 
    alpha = 0.7,
    position=position_dodge(width = 0.6)) +
  theme_cowplot(12) + 
  theme(plot.title = element_text(hjust = 0, size = 12), 
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
    name = "Box Size (Å)")
plot_d

figure_1_all <- plot_grid(plot_c, plot_d, nrow = 1, axis = "t", align="h", labels = c('a', 'b'))
ggsave(filename = paste0("./analysis/figures/figure_1_boxsize.png"), plot = figure_1_all, width = 10, height = 4)
