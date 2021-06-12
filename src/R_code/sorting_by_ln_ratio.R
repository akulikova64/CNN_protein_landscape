library(tidyverse)
library(cowplot)
library(ggforce)

#===============================================================================
# Analysis for protein engineering applications.
#===============================================================================
# Looking at sites where WT has low prob and max_prob is high
#===============================================================================
# ***BOX SIZE 20*****
#===============================================================================
# a) Basic stats (% match across all genes)
#===============================================================================

# set working directory to: "Desktop/Natural_var_project/"
# loading data
cnn_data <- read.csv(file = "./data/PSICOV_box_20/output/cnn_wt_max_freq.csv", header=TRUE, sep=",")
natural_data <- read.csv(file = "./data/PSICOV_box_20/output/natural_max_freq_files/natural_max_freq_all.csv", header=TRUE, sep=",")

joined_data <- rbind(x = cnn_data, y = natural_data)

joined_data_trimmed <- joined_data %>%
  filter(!gene %in% c('1dbx', '1eaz', '1fvg', '1k7j', '1kq6', '1kw4', '1lpy', '1ne2', '1ny1', '1pko', '1rw1', '1vhu', '1w0h', '1wkc', '2tps'))

wide <- joined_data_trimmed %>%
  select(-c(aa_class, class_freq)) %>%
  pivot_wider(names_from = group, values_from = c(freq, aa)) 

ratio <- wide %>%
  mutate(ln_ratio = log(freq_predicted/freq_wt)) %>%
  na.omit() %>%
  arrange(desc(ln_ratio))

get_ratio_bin <- function(x) {
  
  if (x > 1) {
    return(">1")
  }
  if (x < 1 & x != 0) {
    return("<1")
  }
  if (x == 0) {
    return("0")
  }
  
}

counts <- ratio %>%
  mutate(ratio_bin = map(ln_ratio, get_ratio_bin)) %>%
  group_by(gene) %>%
  count(ratio_bin) %>%
  pivot_wider(names_from = ratio_bin, values_from = n)
# looking at these counts, this separation is good enough


ratio_new <- ratio %>%
  mutate(ratio_bin = map(ln_ratio, get_ratio_bin)) %>%
  filter(ratio_bin != "0")

match_cons <- ratio_new %>%
  mutate(match_predict_cons = aa_predicted == aa_natural_max)

for_plot <- match_cons %>%
  select(c(gene, position, ratio_bin, match_predict_cons )) %>%
  group_by(gene, ratio_bin) %>%
  summarise(freq_predict_cons = sum(match_predict_cons, na.rm = TRUE)/sum(!is.na(match_predict_cons))) %>%
  mutate(group = ifelse(ratio_bin == ">1", "high", "low"))

custom_fills <- c("#8c7b9d", "#d2a92d")
custom_colors <- c("#655775", "#8a7228")

data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}

plot_a <- for_plot %>%
  ggplot(aes(y = freq_predict_cons, x = group, fill = group, color = group)) +
  geom_violin(alpha = 0.6, size = 0.7) + 
  stat_summary(fun.data=data_summary, color = "black", alpha = 0.7) +
  theme_cowplot(16) + 
  theme(plot.title = element_text(hjust = 0, size = 16), 
        plot.subtitle = element_text(hjust = 0.5),
        panel.grid.major.y = element_line(color = "grey92", size=0.5),
        legend.position = "none") +
  scale_fill_manual(values = custom_fills) +
  scale_color_manual(values = custom_colors) +
  scale_y_continuous(
    name = "Accuracy",
    limits = c(0.0, 0.55),
    breaks = seq(0, 0.55, by = 0.1),
    expand = c(0, 0)) +
  scale_x_discrete(
    name = "",
    labels = c("ln(max/wt)>1", "ln(max/wt)<1 \n (not 0)"))

plot_a

ggsave(filename = "./analysis/figures/sorted_by_ln.png", plot = plot_a, width = 7, height = 5)

#=================================================================
# Looking at top 10 highest ln() and top 10 lowest ln()
#=================================================================

threshold <- 10

ratio_clean <- ratio %>%
  select(c(gene, position, aa_predicted, aa_natural_max, aa_wt, ln_ratio)) %>%
  mutate(match_wt_cons = aa_wt == aa_natural_max) %>%
  filter(match_wt_cons != TRUE) %>%
  select(-c(match_wt_cons))

ratio <- ratio %>%
  mutate(match_wt_cons = aa_wt == aa_natural_max) %>%
  filter(match_wt_cons != TRUE) %>%
  select(-c(match_wt_cons))

zeros_top <- ratio %>%
  select(c(gene, position, aa_predicted, aa_natural_max, aa_wt, ln_ratio, freq_predicted)) %>%
  filter(ln_ratio == 0) %>%
  group_by(gene) %>%
  slice_max(n=threshold, order_by = freq_predicted, with_ties = FALSE) %>%
  mutate(bin = "zeros_top") %>%
  select(-freq_predicted)

zeros_bottom <- ratio %>%
  select(c(gene, position, aa_predicted, aa_natural_max, aa_wt, ln_ratio, freq_predicted)) %>%
  filter(ln_ratio == 0) %>%
  group_by(gene) %>%
  slice_min(n=threshold, order_by = freq_predicted, with_ties = FALSE) %>%
  mutate(bin = "zeros_bottom") %>%
  select(-freq_predicted)

ratio_filt <- ratio_clean %>%
  filter(ln_ratio != 0)

top_10 <- ratio_filt %>%
  group_by(gene) %>%
  slice_max(n=threshold, order_by = ln_ratio, with_ties = FALSE) %>%
  mutate(bin = "top_10")

low_10 <- ratio_filt %>%
  group_by(gene) %>%
  slice_min(n=threshold, order_by = ln_ratio, with_ties = FALSE) %>%
  mutate(bin = "low_10")

joined_bins <- rbind(top_10, low_10, zeros_top, zeros_bottom)

matches <- joined_bins %>%
  mutate(match_predict_cons = aa_predicted == aa_natural_max)

for_plot_b <- matches %>%
  select(c(gene, position, bin, match_predict_cons )) %>%
  group_by(gene, bin) %>%
  summarise(freq_predict_cons = sum(match_predict_cons, na.rm = TRUE)/sum(!is.na(match_predict_cons)))

custom_fills <- c("#8c7b9d", "#d2a92d", "#80b380", "#80b380")
custom_colors <- c("#655775", "#8a7228", "#396039", "#396039")


plot_b <- for_plot_b %>%
  ggplot(aes(y = freq_predict_cons, x = fct_relevel(bin, "top_10", "low_10", "zeros_bottom", "zeros_top"), fill = bin, color = bin)) +
  geom_violin(alpha = 0.6, size = 0.7, bw = 0.05) + 
  stat_summary(fun.data=data_summary, color = "black", alpha = 0.7) +
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
    labels = c(paste("top",threshold,"\n ln(max/wt) > 0"), paste("bottom",threshold,"\n ln(max/wt) ≠ 0 "), paste("bottom",threshold,"max \n ln(max/wt) = 0"), paste("top",threshold,"max \n ln(max/wt) = 0")))

plot_b

ggsave(filename = paste("./analysis/figures/sorted_top_",threshold,".png"), plot = plot_b, width = 8, height = 4)

