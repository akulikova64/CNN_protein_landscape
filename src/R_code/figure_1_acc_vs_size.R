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

stat_data_1 <- stats_1 %>%
  select(-c(gene, x_label)) %>%
  group_by(box_size, group) %>%
  summarise(estimate = mean(freq_predict_wt),
            std_error = sd(freq_predict_wt)/sqrt(length(freq_predict_wt)))

plot_a <- stats_1 %>%
  ggplot(aes(y = freq_predict_wt, x = box_size)) +
  geom_violin(alpha = 0.5, size = 0.7, fill = custom_fill_1, color = custom_color_1 ) +
  #stat_summary(fun.data=data_summary, color = "black", alpha = 0.7) +
  geom_pointrange(data = stat_data_1, aes(x = box_size,
                                          y = estimate,
                                          ymin = estimate - std_error,
                                          ymax = estimate + std_error),
                  color = "black", alpha = 0.7, size = 0.3) +
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

stat_data_2 <- joined_single_and_class %>%
  select(-c(gene, x_label)) %>%
  group_by(box_size, group) %>%
  summarise(estimate = mean(freq_predict_wt),
            std_error = sd(freq_predict_wt)/sqrt(length(freq_predict_wt)))

plot_c <- joined_single_and_class %>%
  ggplot(aes(y = freq_predict_wt, x = box_size, fill = group, color = group)) +
  geom_violin(
    alpha = 0.5, 
    size = 0.7,
    position=position_dodge(width = 0.6)) +
  # stat_summary(
  #   fun.data=data_summary, 
  #   color = "black", 
  #   alpha = 0.7,
  #   position=position_dodge(width = 0.6)) +
  geom_pointrange(data = stat_data_2, aes(x = box_size,
                                          y = estimate,
                                          ymin = estimate - std_error,
                                          ymax = estimate + std_error),
                  color = "black", alpha = 0.7, size = 0.3,
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

stat_data_3 <- joined_consensus %>%
  select(-c(gene, x_label)) %>%
  group_by(box_size, group) %>%
  summarise(estimate = mean(freq_predict_cons),
            std_error = sd(freq_predict_cons)/sqrt(length(freq_predict_cons)))

plot_d <- joined_consensus %>%
  ggplot(aes(y = freq_predict_cons, x = box_size, fill = group, color = group)) +
  geom_violin(
    alpha = 0.5, 
    size = 0.7,
    position=position_dodge(width = 0.6)) +
  # stat_summary(
  #   fun.data=data_summary, 
  #   color = "black", 
  #   alpha = 0.7,
  #   position=position_dodge(width = 0.6)) +
  geom_pointrange(data = stat_data_3, aes(x = box_size,
                                          y = estimate,
                                          ymin = estimate - std_error,
                                          ymax = estimate + std_error),
                  color = "black", alpha = 0.7, size = 0.3,
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


#=======================================================================================
# Performing a paired t-test between box sizes (predicting the wild type)
#=======================================================================================


# for single amino acids:
for_ttests <- stats_1 %>%
  select(c(gene, box_size, freq_predict_wt))

averages <- for_ttests %>%
  group_by(box_size) %>%
  summarise(mean = mean(freq_predict_wt))
averages

for_paired_test <- for_ttests %>%
  pivot_wider(names_from = box_size, values_from = freq_predict_wt) %>%
  rename("box_12" = `12`,
         "box_20" = `20`,
         "box_30" = `30`,
         "box_40" = `40`)

results_1 <- t.test(for_paired_test$box_12, for_paired_test$box_20, paired = TRUE, alternative = "two.sided")
results_1
#p-value = 1.108e-10 (there is a significant difference)

results_2 <- t.test(for_paired_test$box_12, for_paired_test$box_30, paired = TRUE, alternative = "two.sided")
results_2
#p-value = 2.2e-16 (there is a significant difference)

results_3 <- t.test(for_paired_test$box_12, for_paired_test$box_40, paired = TRUE, alternative = "two.sided")
results_3
#p-value = 7.738e-13 (there is a significant difference)

results_4 <- t.test(for_paired_test$box_20, for_paired_test$box_30, paired = TRUE, alternative = "two.sided")
results_4
#p-value = 9.121e-05 (there is a significant difference)

results_5 <- t.test(for_paired_test$box_20, for_paired_test$box_40, paired = TRUE, alternative = "two.sided")
results_5
#p-value = 0.137 (there is no difference)

results_6 <- t.test(for_paired_test$box_30, for_paired_test$box_40, paired = TRUE, alternative = "two.sided")
results_6
#p-value = 0.0073 (there is a significant difference)

#==================================================================================================
# for classes:

for_ttests_2 <- stats_2 %>%
  select(c(gene, box_size, freq_predict_wt))

averages_2 <- for_ttests_2 %>%
  group_by(box_size) %>%
  summarise(mean = mean(freq_predict_wt))
averages_2
#

for_paired_test <- for_ttests_2 %>%
  pivot_wider(names_from = box_size, values_from = cor) %>%
  rename("box_12" = `12`,
         "box_20" = `20`,
         "box_30" = `30`,
         "box_40" = `40`)


results_1b <- t.test(for_paired_test$box_12, for_paired_test$box_20, paired = TRUE, alternative = "two.sided")
results_1b
#p-value = 0.4402 (no difference)

results_2b <- t.test(for_paired_test$box_12, for_paired_test$box_30, paired = TRUE, alternative = "two.sided")
results_2b
#p-value = 0.0093 (there is a significant difference)

results_3b <- t.test(for_paired_test$box_12, for_paired_test$box_40, paired = TRUE, alternative = "two.sided")
results_3b
#p-value = 3.5730e-06 (there is a significant difference)

results_4b <- t.test(for_paired_test$box_20, for_paired_test$box_30, paired = TRUE, alternative = "two.sided")
results_4b
#p-value = 0.0276 (there is a significant difference)

results_5b <- t.test(for_paired_test$box_20, for_paired_test$box_40, paired = TRUE, alternative = "two.sided")
results_5b
#p-value = 3.1390e-05 (there is a significant difference)

results_6b <- t.test(for_paired_test$box_30, for_paired_test$box_40, paired = TRUE, alternative = "two.sided")
results_6b
#p-value = 0.0067 (there is a significant difference)

#=======================================================================================
# Performing a paired t-test between box sizes (predicting the consensus)
#=======================================================================================


# for single amino acids:
for_ttests <- stats_3 %>%
  select(c(gene, box_size, freq_predict_cons))

averages <- for_ttests %>%
  group_by(box_size) %>%
  summarise(mean = mean(freq_predict_cons))
averages

for_paired_test <- for_ttests %>%
  pivot_wider(names_from = box_size, values_from = freq_predict_cons) %>%
  rename("box_12" = `12`,
         "box_20" = `20`,
         "box_30" = `30`,
         "box_40" = `40`)

results_1 <- t.test(for_paired_test$box_12, for_paired_test$box_20, paired = TRUE, alternative = "two.sided")
results_1
#p-value = 0.0026 (there is a significant difference)

results_2 <- t.test(for_paired_test$box_12, for_paired_test$box_30, paired = TRUE, alternative = "two.sided")
results_2
#p-value = 6.568e-06 (there is a significant difference)

results_3 <- t.test(for_paired_test$box_12, for_paired_test$box_40, paired = TRUE, alternative = "two.sided")
results_3
#p-value = 0.000149 (there is a significant difference)

results_4 <- t.test(for_paired_test$box_20, for_paired_test$box_30, paired = TRUE, alternative = "two.sided")
results_4
#p-value = 0.1082 (no difference)

results_5 <- t.test(for_paired_test$box_20, for_paired_test$box_40, paired = TRUE, alternative = "two.sided")
results_5
#p-value = 0.1678 (no difference)

results_6 <- t.test(for_paired_test$box_30, for_paired_test$box_40, paired = TRUE, alternative = "two.sided")
results_6
#p-value = 0.8333 (no difference)

#==================================================================================================
# for classes:

for_ttests_2 <- stats_4 %>%
  select(c(gene, box_size, freq_predict_cons))

averages_2 <- for_ttests_2 %>%
  group_by(box_size) %>%
  summarise(mean = mean(freq_predict_cons))
averages_2


for_paired_test <- for_ttests_2 %>%
  pivot_wider(names_from = box_size, values_from = freq_predict_cons) %>%
  rename("box_12" = `12`,
         "box_20" = `20`,
         "box_30" = `30`,
         "box_40" = `40`)


results_1b <- t.test(for_paired_test$box_12, for_paired_test$box_20, paired = TRUE, alternative = "two.sided")
results_1b
#p-value = 0.2573 (no difference)

results_2b <- t.test(for_paired_test$box_12, for_paired_test$box_30, paired = TRUE, alternative = "two.sided")
results_2b
#p-value = 0.000309 (there is a significant difference)

results_3b <- t.test(for_paired_test$box_12, for_paired_test$box_40, paired = TRUE, alternative = "two.sided")
results_3b
#p-value = 0.0075516 (there is a significant difference)

results_4b <- t.test(for_paired_test$box_20, for_paired_test$box_30, paired = TRUE, alternative = "two.sided")
results_4b
#p-value = 0.002141 (there is a significant difference)

results_5b <- t.test(for_paired_test$box_20, for_paired_test$box_40, paired = TRUE, alternative = "two.sided")
results_5b
#p-value = 0.06939 (no difference)

results_6b <- t.test(for_paired_test$box_30, for_paired_test$box_40, paired = TRUE, alternative = "two.sided")
results_6b
#p-value = 0.2147 (no difference)



