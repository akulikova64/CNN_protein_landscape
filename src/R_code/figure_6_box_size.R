# comparing neff predicted vs. neff natural across different box sizes for similarity group (60-80%)
library(tidyverse)
library(cowplot)
library(broom)
library(ggforce)

# useful function for getting mean and standard of deviation (for violin plots):
data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}

#set working directory to:
#C:\Users\avch\Desktop\Natural_var_project

# reading csv files
#only using the 40-60% similarity group
natural_var <- read.csv(file = "./output/output_PSICOV/stats_align_files/stats_align_60.csv", header=TRUE, sep=",")

cnn_var_12 <- read.csv(file= paste0("./data/PSICOV_box_12/output/stats_cnn.csv"), header=TRUE, sep=",")
cnn_var_20 <- read.csv(file= paste0("./data/PSICOV_box_20/output/stats_cnn.csv"), header=TRUE, sep=",")
cnn_var_30 <- read.csv(file= paste0("./data/PSICOV_box_30/output/stats_cnn.csv"), header=TRUE, sep=",")
cnn_var_40 <- read.csv(file= paste0("./data/PSICOV_box_40/output/stats_cnn.csv"), header=TRUE, sep=",")


natural_var <- natural_var %>%
  select(position, gene, n_eff, n_eff_class) %>%
  mutate(group = "natural")

#box_size 12
cnn_var_12 <- cnn_var_12 %>%
  select(position, gene, n_eff, n_eff_class) %>%
  mutate(group = "predicted")

joined_12 <- rbind(natural_var, cnn_var_12) %>%
  mutate(box_size = "12") 

#box_size 20
cnn_var_20 <- cnn_var_20 %>%
  select(position, gene, n_eff, n_eff_class) %>%
  mutate(group = "predicted")

joined_20 <- rbind(natural_var, cnn_var_20) %>%
  mutate(box_size = "20") 

#box_size 30
cnn_var_30 <- cnn_var_30 %>%
  select(position, gene, n_eff, n_eff_class) %>%
  mutate(group = "predicted")

joined_30 <- rbind(natural_var, cnn_var_30) %>%
  mutate(box_size = "30") 

#box_size 40
cnn_var_40 <- cnn_var_40 %>%
  select(position, gene, n_eff, n_eff_class) %>%
  mutate(group = "predicted")

joined_40 <- rbind(natural_var, cnn_var_40) %>%
  mutate(box_size = "40") 

#joining all the box_sizes
all_joined <- rbind(joined_12, joined_20, joined_30, joined_40)

all_joined_wide <- all_joined %>%
  select(-n_eff_class) %>%
  pivot_wider(names_from = group, values_from = n_eff)

# fitting a linear model to data (getting R^2 and p-values)
lm_summary <- all_joined_wide %>%
  na.omit() %>%
  nest(data = -c(gene, box_size)) %>%
  mutate(
    fit = map(data, ~lm(natural ~ predicted, data = .x)),
    glance_out = map(fit, glance)
  ) %>%
  select(gene, box_size, glance_out) %>%
  unnest(cols = glance_out)
lm_summary

# getting rid of the genes that are not present in all 4 groups:
lm_summary <- lm_summary %>%
  select(gene, box_size, r.squared, p.value) 

# making two dataframes for correlation coefficients and p-values:
cor_coeffs <- lm_summary %>%
  select(gene, box_size, r.squared)

p_values <- lm_summary %>%
  select(gene, box_size, p.value)

new_p_values <- p_values %>%
  nest(data = -c(box_size)) %>%
  mutate(
    new_p = map(data, ~p.adjust(.x$p.value, method = "fdr", n = length(.x$p.value)))
  ) %>%
  unnest(cols = c(data, new_p))
new_p_values

# removing genes that are not found across all 5 similarity groups:
cor_coeffs_wide <- cor_coeffs %>%
  pivot_wider(names_from = box_size, values_from = r.squared)

#here is where the error is. 
p_values_wide <- new_p_values %>%
  select(-p.value) %>%
  pivot_wider(names_from = box_size, values_from = new_p)
p_values_wide

cor_coeffs_reduced <- na.omit(cor_coeffs_wide)
p_values_reduced <- na.omit(p_values_wide)

# making both dataframes longer again for plotting:
cor_coeffs <- cor_coeffs_reduced %>%
  pivot_longer(
    cols = -gene, 
    names_to = "box_size", 
    values_to = c("r_squared"))

p_values <- p_values_wide %>%
  pivot_longer(
    cols = -gene, 
    names_to = "box_size", 
    values_to = c("p_value"))

# labeling the significant p-values:
p_values <- p_values %>%
  mutate(signif = ifelse(p_value <= 0.05, TRUE, FALSE))

#FINDING CORRELATION COEFFS:
# an alternative method for finding the correlations:
cor <- all_joined_wide %>%
  na.omit() %>%
  group_by(gene, box_size) %>%
  summarise(cor = cor(natural, predicted)) 


# filtering for the correlations that are significant:
joined_cors_2 <- inner_join(p_values, cor)

#comparing all_data_1 (from figure_6.R) and all_data_2 from this script (rerun figure_6.R to get the first joined_cors_1)

reduced <- right_join(joined_cors_1, joined_cors_2)
reduced <- na.omit(reduced) 
reduced <- reduced %>%
  select(gene) %>%
  distinct() %>%
  mutate(flag = TRUE)

for_all <- left_join(cor, reduced, by = "gene")
for_all <- na.omit(for_all) 
for_all <- for_all %>%
  select(-flag)

joined_final <- left_join(joined_cors_2, reduced, by = "gene")
joined_final <- na.omit(joined_final) 
joined_final <- joined_final %>%
  select(-flag)

joined_final <- joined_final%>%
  group_by(gene) %>%
  mutate(
    # pick y value corresponding to y3
    color_y = sum(cor * (box_size == "12")),
    dx = rnorm(n(), mean = 0, sd = .05),
    dy = rnorm(n(), mean = 0, sd = .05),
    x_value = as.numeric(factor(box_size)))

# filtering for the correlations that are significant:
sig_cor <- joined_final %>%
  filter(signif == TRUE) %>%
  select(cor, box_size, gene, dx, dy, x_value, color_y)

# filtering for the correlations that are **NOT** significant:
not_signif <- joined_final %>%
  filter(signif == FALSE) %>%
  select(cor, box_size, gene, dx, dy, x_value, color_y)

#===========================================================================================
# plot_a (correlations bw neff predicted and neff natural across different box sizes)
#=============================================================================================
# sig_cor <- sig_cor %>%
#   group_by(gene) %>%
#   mutate(
#     # pick y value corresponding to y3
#     color_y = sum(cor * (box_size == "12"))
#   ) 
# 
# all_data <- for_all %>%
#   group_by(gene) %>%
#   mutate(
#     # pick y value corresponding to y3
#     color_y = sum(cor * (box_size == "12"))
#   ) 


plot_a <- ggplot() +
  geom_path(
    data = not_signif,
    aes(x = as.numeric(factor(box_size))+dx, y = cor+dy, group = gene, color = color_y),
    size = 0.25) +
  geom_path(
    data = sig_cor,
    aes(x = as.numeric(factor(box_size))+dx, y = cor+dy, group = gene, color = color_y),
    size = 0.25) +
  geom_point(
    data = not_signif,
    aes(x = as.numeric(factor(box_size))+dx, y = cor+dy, group = gene),
    shape = 21, 
    color = "black",
    fill = "white",
    size = 2) +
  geom_point(
    data = sig_cor,
    aes(x = as.numeric(factor(box_size))+dx, y = cor+dy, group = gene, fill = color_y),
    shape = 21, 
    color = "black",
    size = 2) +
  scale_x_continuous(
    name = "Box Size (Å)",
    limits = c(0.5, 4.5),
    labels = c("12", "20", "30", "40"),
    breaks = (seq(from = 1, to = 4, by = 1)),
    expand = c(0, 0)) +
  scale_y_continuous(
    name = "Correlation Coefficients",
    limits = c(-0.2, 0.7),
    breaks = seq(from = -0.2, to = 0.6, by = 0.1),
    expand = c(0, 0)) +
  scale_color_gradient(
    aesthetics = c("color", "fill"), 
    high = "#ffd966", 
    low = "#080845") +
  theme_bw(16) +
  theme(
    legend.position = "none",
    axis.text = element_text(color = "black", size = 16),
    panel.grid.minor = element_blank())

plot_a

ggsave(filename = paste0("./analysis/figures/figure_6_box_size.png"), plot = plot_a, width = 8, height = 4)

#=======================================================================================
# Performing a paired t-test between the 12A box and 20A box. 
#=======================================================================================

for_ttests <- joined_final %>%
  select(c(gene, box_size, cor))

averages <- for_ttests %>%
  group_by(box_size) %>%
  summarise(mean = mean(cor))


for_paired_test <- for_ttests %>%
  pivot_wider(names_from = box_size, values_from = cor) %>%
  rename("box_12" = `12`,
         "box_20" = `20`,
         "box_30" = `30`,
         "box_40" = `40`)


results_1 <- t.test(for_paired_test$box_12, for_paired_test$box_20, paired = TRUE, alternative = "two.sided")
results_1
#p-value = 0.4402 (can't reject null, no difference)

results_2 <- t.test(for_paired_test$box_12, for_paired_test$box_30, paired = TRUE, alternative = "two.sided")
results_2
#p-value = 0.0093 (there is a significant difference)

results_3 <- t.test(for_paired_test$box_12, for_paired_test$box_40, paired = TRUE, alternative = "two.sided")
results_3
#p-value = 3.5730e-06 (there is a significant difference)

results_4 <- t.test(for_paired_test$box_20, for_paired_test$box_30, paired = TRUE, alternative = "two.sided")
results_4
#p-value = 0.0276 (there is a significant difference)

results_5 <- t.test(for_paired_test$box_20, for_paired_test$box_40, paired = TRUE, alternative = "two.sided")
results_5
#p-value = 3.1390e-05 (there is a significant difference)

results_6 <- t.test(for_paired_test$box_30, for_paired_test$box_40, paired = TRUE, alternative = "two.sided")
results_6
#p-value = 0.0067 (there is a significant difference)


#===========================================================================================
# correlating n_eff natural vs. predicted for all positions (not by protein)
#=============================================================================================

natural_var <- read.csv(file = "./output/output_PSICOV/stats_align_files/stats_align_60.csv", header=TRUE, sep=",")

cnn_var_12 <- read.csv(file= paste0("./data/PSICOV_box_12/output/stats_cnn.csv"), header=TRUE, sep=",")
cnn_var_20 <- read.csv(file= paste0("./data/PSICOV_box_20/output/stats_cnn.csv"), header=TRUE, sep=",")
cnn_var_30 <- read.csv(file= paste0("./data/PSICOV_box_30/output/stats_cnn.csv"), header=TRUE, sep=",")
cnn_var_40 <- read.csv(file= paste0("./data/PSICOV_box_40/output/stats_cnn.csv"), header=TRUE, sep=",")


natural_var <- natural_var %>%
  select(position, gene, n_eff, n_eff_class) %>%
  mutate(group = "natural")

#box_size 12
cnn_var_12 <- cnn_var_12 %>%
  select(position, gene, n_eff, n_eff_class) %>%
  mutate(group = "predicted")

joined_12 <- rbind(natural_var, cnn_var_12) %>%
  mutate(box_size = "12") 

#box_size 20
cnn_var_20 <- cnn_var_20 %>%
  select(position, gene, n_eff, n_eff_class) %>%
  mutate(group = "predicted")

joined_20 <- rbind(natural_var, cnn_var_20) %>%
  mutate(box_size = "20") 

#box_size 30
cnn_var_30 <- cnn_var_30 %>%
  select(position, gene, n_eff, n_eff_class) %>%
  mutate(group = "predicted")

joined_30 <- rbind(natural_var, cnn_var_30) %>%
  mutate(box_size = "30") 

#box_size 40
cnn_var_40 <- cnn_var_40 %>%
  select(position, gene, n_eff, n_eff_class) %>%
  mutate(group = "predicted")

joined_40 <- rbind(natural_var, cnn_var_40) %>%
  mutate(box_size = "40") 

#joining all the box_sizes
all_joined <- rbind(joined_12, joined_20, joined_30, joined_40)

all_joined_wide <- all_joined %>%
  select(-n_eff_class) %>%
  pivot_wider(names_from = group, values_from = n_eff)

# fitting a linear model to data (getting R^2 and p-values)
lm_summary <- all_joined_wide %>%
  na.omit() %>%
  nest(data = -c(box_size)) %>%
  mutate(
    fit = map(data, ~lm(natural ~ predicted, data = .x)),
    glance_out = map(fit, glance)
  ) %>%
  select(box_size, glance_out) %>%
  unnest(cols = glance_out)
lm_summary

get_n_eff_bin <- function(x) {
  
  if (x > 0.0 & x <= 4.0) {
    return("(0-4]")
  }
  else if (x > 4.0 & x <= 8.0) {
    return("(4-8]")
  }
  else if (x > 8.0 & x <= 12.0) {
    return("(8-12]")
  }
  else if (x > 12.0 & x <= 16.0) {
    return("(12-16]")
  }
  else if (x > 16.0 & x <= 20.0) {
    return("(16-20]")
  }
}

with_bin <- all_joined_wide %>%
  na.omit() %>%
  mutate(pred_n_eff_bin = map_chr(predicted, get_n_eff_bin))

stat_data <- with_bin %>%
  select(-c(position, gene)) %>%
  group_by(pred_n_eff_bin, box_size) %>%
  summarise(estimate = mean(natural),
            std_error = sd(natural)/sqrt(length(natural)))
  

plot_n_eff <- with_bin %>%
  ggplot(aes(y = natural, x = fct_relevel(pred_n_eff_bin, "(0-4]", "(4-8]", "(8-12]", "(12-16]", "(16-20]" ))) +
  geom_sina(alpha = 0.6, size = 0.2, bw = 0.3, fill = "#99a88a", color = "#4d5841") +
  geom_pointrange(data = stat_data, aes(x = pred_n_eff_bin,
                                        y = estimate,
                                        ymin = estimate - 1.96*std_error,
                                        ymax = estimate + 1.96*std_error),
                  color = "black", alpha = 0.7, size = 0.4) +
  #stat_summary(fun.data=data_summary, color = "black", alpha = 0.7) +
  facet_wrap(vars(box_size)) +
  theme_cowplot(16) + 
  theme(plot.title = element_text(hjust = 0, size = 16), 
        plot.subtitle = element_text(hjust = 0.5),
        panel.grid.major.y = element_line(color = "grey92", size=0.5),
        legend.position = "none") +
  scale_y_continuous(
    name = "Natural variation (n-eff)",
    limits = c(0, 20),
    breaks = seq(0, 20, by = 2),
    expand = c(0, 0)) +
  scale_x_discrete(
    name = "Predicted variation (n-eff)")

plot_n_eff


ggsave(filename = "./analysis/figures/n_eff_vs_n_eff.png", plot = plot_n_eff, width = 8, height = 5)

