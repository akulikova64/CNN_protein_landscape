library(tidyverse)
library(yardstick)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(cowplot)
library(tidyr)
library(sqldf)
library(sinaplot)
library(ggforce)

cnn_data <- read.csv(file="./cnn_wt_max_freq.csv", header=TRUE, sep=",")
natural_data <- read.csv(file = "./natural_max_freq.csv", header=TRUE, sep=",")

joined_data <- rbind(x = cnn_data, y = natural_data)

joined_data %>%
  filter(gene == "1b4t") %>%
  #filter(group == c("predicted", "wt", "natural_max")) %>%
  filter(group == "predicted" | group == "wt" | group == "natural_max") %>%
  ggplot(aes(x = position, y = freq, color = group)) +
  geom_line() +
  geom_point() + geom_text(aes(label = aa), hjust = -1, vjust = 0) +
  ggtitle(label = "1b4t") +
  theme_cowplot() +
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) +
  ylab("frequency") +
  xlab("position") +
  coord_cartesian(xlim = c(0, 30)) +
  scale_x_continuous(breaks = seq(0, 30, by = 2)) +
  grids(axis = "x", color = "grey92", linetype = "dashed")

index_data <- joined_data %>%
  pivot_wider(names_from = group, values_from = c(aa,freq)) %>%
  mutate(match = aa_predicted != aa_natural_wt) %>%
  select(gene, position, match)

# selecting the data entries where the predicted amino acid matches the 
filtered_data <- joined_data %>%
  left_join(index_data) %>%
  filter(match == TRUE) %>%
  select(-match)


# pred <- filtered_data %>%
#   select(gene, position, aa = aa_predicted, freq = freq_predicted) %>%
#   mutate(group = "predicted")
# 
# wt <- filtered_data %>%
#   select(gene, position, aa = aa_wt, freq = freq_wt) %>%
#   mutate(group = "wt")
# 
# nat <- filtered_data %>%
#   select(gene, position, aa = aa_natural, freq = freq_natural) %>%
#   mutate(group = "natural")
# 
# filtered <- rbind(pred, wt, nat)
  
# pivot_longer(c(freq_predicted, freq_wt, freq_natural), names_to = "group", values_to = c("freq")) %>%
# pivot_longer(c(aa_predicted, aa_wt, aa_natural), names_to = "group", values_to = c("aa"), names_repair= "minimal")

filtered_data %>%
  filter(gene == "1b4t") %>%
  filter(group == "predicted" | group == "wt" | group == "natural_max") %>%
  ggplot(aes(x = position, y = freq, color = group)) +
  geom_line() +
  geom_point() + geom_text(aes(label = aa), hjust = -1, vjust = 0) +
  ggtitle(label = "1b4t") +
  theme_cowplot() +
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) +
  ylab("frequency") +
  xlab("position") +
  coord_cartesian(xlim = c(0, 30)) +
  scale_x_continuous(breaks = seq(0, 30, by = 2)) +
  grids(axis = "x", color = "grey92", linetype = "dashed")

#rerun above code to get original joined_data    
joined_data %>%
  filter(gene == "1b4t") %>%
  ggplot(aes(x = freq_predicted, y = freq_wt)) +
  geom_point() +
  ggtitle(label = "1b4t") +
  theme_cowplot() +
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) +
  ylab("frequency wt") +
  xlab("frequency predicted") 

index_data2 <- joined_data %>%
  pivot_wider(names_from = group, values_from = c(aa,freq))

#total freq of wt in natural alignment vs cnn wt prob. 
index_data2 %>%
  filter(gene == "1b4t") %>%
  ggplot(aes(x = freq_wt, y = freq_natural_wt)) +
  geom_point() +
  ggtitle(label = "1b4t") +
  theme_cowplot() +
  geom_abline(slope = 1, color = "grey65") +
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) +
  ylab("frequency natural of wt") +
  xlab("frequency wt") 

index_data3 <- filtered_data %>%
  pivot_wider(names_from = group, values_from = c(aa,freq))

all_data <- joined_data %>%
  pivot_wider(names_from = group, values_from = c(aa,freq))

#----------------------------------------------------------------------------------------------------
# 1) Probability of wt residue (CNN) ~ Frequency of wt residue (seq alignment)  
goi = "1b4t"

#finding the linear regression
linear_reg <- all_data %>%
  filter(gene == goi) %>%
  select(freq_wt, freq_natural_wt) %>%
  lm(formula = freq_wt ~ freq_natural_wt)
summary(linear_reg)

slope <- linear_reg$coefficients[2]
intercept <- linear_reg$coefficients[1]

all_data %>%
  filter(gene == goi) %>%
  ggplot(aes(x = freq_wt, y = freq_natural_wt)) +
  geom_point() +
  ggtitle(label = goi) +
  geom_abline(slope = slope, intercept = intercept, color = "grey65" ) +
  stat_cor(label.y = 0.95, label.x = -0.1) +
  stat_regline_equation(label.y = .85, label.x = -0.1) +
  theme_cowplot() +
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) +
  ylab("wt frequency (alignment)") +
  xlab("wt probability (cnn)") 

library(dplyr)
cor <- all_data %>%
  select(gene, freq_wt, freq_natural_wt) %>%
  na.omit() %>%
  group_by(gene) %>%
  summarize(cor = cor(freq_wt, freq_natural_wt))

cor %>%
  ggplot(aes(x = "", y = cor)) +
  geom_violin(fill = "grey75") + geom_sina() +
  ggtitle(label = "wt prob (CNN) ~ wt freq (alignment)") +
  theme_cowplot() +
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) +
  xlab("") +
  ylab("correlation coefficients") 

#---------------------------------------------------------------------------------------------
# 2) Probability of wt residue (CNN) ~ Neff or entropy (seq alignment)

natural_var <- read.csv(file="./stats_align_all.csv", header=TRUE, sep=",")
gene_name <- "1w7w"

index_data4 <- all_data %>%
  select(position, gene, freq_wt) %>%
  left_join(natural_var) %>%
  select(position, gene, freq_wt, n_eff)

linear_reg <- index_data4 %>%
  filter(gene == gene_name) %>%
  select(freq_wt, n_eff) %>%
  lm(formula = n_eff ~ freq_wt)
summary(linear_reg)

slope <- linear_reg$coefficients[2]
intercept <- linear_reg$coefficients[1]

index_data4 %>%
  filter(gene == gene_name) %>%
  ggplot(aes(x = freq_wt, y = n_eff)) +
  geom_point() +
  ggtitle(label = gene_name) +
  geom_abline(slope = slope, intercept = intercept, color = "grey65" ) +
  stat_cor(label.x = 0, label.y = 11) +
  stat_regline_equation(label.x = 0, label.y = 10) +
  theme_cowplot() +
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) +
  ylab("n-eff (alignment)") +
  xlab("wt probability (cnn)") 

library(dplyr)
cor2 <- index_data4 %>%
  select(gene, freq_wt, n_eff) %>%
  na.omit() %>%
  group_by(gene) %>%
  summarize(cor = cor(freq_wt, n_eff))

cor2 %>%
  ggplot(aes(x = "", y = cor)) +
  geom_violin(fill = "grey75") + geom_sina() + 
  ggtitle(label = "wt prob (CNN) ~ n-eff (alignment)") +
  theme_cowplot() +
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) +
  xlab("") +
  ylab("correlation coefficients") 

  # paired sample t-test
  t_test_1 <- t.test(before, after, paired = FALSE)

#----------------------------------------------------------------------------------------------
# 3) Probability of wt class (CNN) ~ Frequency of wt class (seq alignment)

natural_var <- read.csv(file="./stats_align_all.csv", header=TRUE, sep=",")
cnn_var <- read.csv(file="./stats_cnn.csv", header=TRUE, sep=",")

natural_var <- natural_var %>%
  nest(q_natural = c(q_H:q_C))

cnn_var <- cnn_var %>%
  nest(q_cnn = c(q_H:q_C))

# looking at classes
joined_data3 <- inner_join(x = natural_var, y = cnn_var, by = c('position', 'gene')) 

#natural_var and cnn_var from variability_38_genes
classes <- joined_data3 %>%
  mutate(proline = map(q_natural, function(x) x$q_P),
         aromatic = map(q_natural, function(x) x$q_F + x$q_Y + x$q_W),
         negative = map(q_natural, function(x) x$q_D + x$q_E),
         positive = map(q_natural, function(x) x$q_K + x$q_R + x$q_H),
         aliphatic = map(q_natural, function(x) x$q_G + x$q_A + x$q_V + x$q_L + x$q_M + x$q_I),
         polar = map(q_natural, function(x) x$q_S + x$q_T + x$q_C + x$q_N + x$q_Q)
  ) %>%
  nest(q_nat_class = c(proline:polar)) %>%
  
  mutate(proline = map(q_cnn, function(x) x$q_P),
         aromatic = map(q_cnn, function(x) x$q_F + x$q_Y + x$q_W),
         negative = map(q_cnn, function(x) x$q_D + x$q_E),
         positive = map(q_cnn, function(x) x$q_K + x$q_R + x$q_H),
         aliphatic = map(q_cnn, function(x) x$q_G + x$q_A + x$q_V + x$q_L + x$q_M + x$q_I),
         polar = map(q_cnn, function(x) x$q_S + x$q_T + x$q_C + x$q_N + x$q_Q)
  ) %>%
  nest(q_cnn_class = c(proline:polar)) %>%
  
  mutate(wt_class = map(wt_aa, function(q){
  if (q %in% c("GLY", "ALA", "VAL", "LEU", "MET", "ILE")) return ("aliphatic")
  if (q %in% c("SER", "THR", "CYS", "ASN", "GLN")) return ("polar")
  if (q %in% c("LYS", "ARG", "HIS")) return ("positive")
  if (q %in% c("ASP", "GLU")) return ("negative")
  if (q %in% c("PHE", "TYR", "TRP")) return ("aromatic")
  if (q %in% c("PRO")) return ("proline")
}))

linear_reg <- wt_nat_wt_cnn %>%
  filter(gene == gene_name) %>%
  select(wt_class_nat, wt_class_cnn) %>%
  lm(formula = wt_class_nat ~ wt_class_cnn)
summary(linear_reg)

slope <- linear_reg$coefficients[2]
intercept <- linear_reg$coefficients[1]

wt_nat_wt_cnn <- classes %>%
  mutate(nat = map2(q_nat_class, wt_class, function(x, y) x[,y,drop=TRUE][[1]]),
         cnn = map2(q_cnn_class, wt_class, function(x, y) x[,y,drop=TRUE][[1]])
         ) %>%
  #select(gene, wt_class_nat, wt_class_cnn)
  mutate(wt_class_nat = sapply(nat, function(x) x[[1]])) %>%
  mutate(wt_class_cnn = sapply(cnn, function(x) x[[1]]))
  #mutate(jblp = mapply(function(x,y) x[[1]] + y[[1]], wt_class_nat, wt_class_cnn))

wt_nat_wt_cnn %>%  
  filter(gene == "1b4t") %>%
  ggplot(aes(x = wt_class_cnn, y = wt_class_nat)) +
  geom_point() +
  ggtitle(label = "1b4t") +
  stat_cor(label.x = 0, label.y = 1) +
  stat_regline_equation(label.x = 0, label.y = .9) +
  geom_abline(slope = slope, intercept = intercept, color = "grey65" ) +
  theme_cowplot() +
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) +
  ylab("wt class frequency (alignment)") +
  xlab("wt class probability (cnn)") 

# looking at any class changes  
max_cnn_wt_nat <- classes %>%
  mutate(max_cnn = map(q_cnn_class, function(x) max(as.numeric(as.character(unlist(x))))),
         wt_cnn = map2(q_cnn_class, wt_class, function(x, y) x[,y,drop=TRUE][[1]])
  ) %>%
  mutate(max_cnn = sapply(max_cnn, function(x) x[[1]])) %>%
  mutate(wt_cnn = sapply(wt_cnn, function(x) x[[1]])) %>%
  filter(max_cnn != wt_cnn) %>%
  pivot_longer(c(max_cnn, wt_cnn), names_to = "group", values_to = "freq")

max_cnn_wt_nat %>%
  filter(gene == "1b4t") %>%
  #select(position, max_cnn, wt_cnn) %>%
  ggplot(aes(x = position, y = freq, color = group)) +
  geom_line() +
  geom_point() + #geom_text(aes(label = ), hjust = -1, vjust = 0) +
  ggtitle(label = "1b4t") +
  theme_cowplot() +
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) +
  ylab("frequency") +
  xlab("position") +
  coord_cartesian(xlim = c(0, 30)) +
  scale_x_continuous(breaks = seq(0, 30, by = 2)) +
  grids(axis = "x", color = "grey92", linetype = "dashed")

#now lets use joined data from above
classes2 <- joined_data %>%
  mutate(class = map(aa, function(q){
    if (q %in% c("G", "A", "V", "L", "M", "I")) return ("aliphatic")
    if (q %in% c("S", "T", "C", "N", "Q")) return ("polar")
    if (q %in% c("K", "R", "H")) return ("positive")
    if (q %in% c("D", "E")) return ("negative")
    if (q %in% c("F", "Y", "W")) return ("aromatic")
    if (q %in% c("P")) return ("proline")
  }))

#counting up the number of class mismatches (cnn max vs cnn wt)
count_mismatches <- classes %>%
  mutate(max_cnn = map(q_cnn_class, function(x) colnames(x, max(as.numeric(as.character(unlist(x))))))
  #mutate(max_cnn = map(q_cnn_class, function(x) unlist(x))
         ) 
  mutate(max_cnn = sapply(max_cnn, function(x) x[[1]]))
  filter(max_cnn != wt_cnn)

count_mismatches <- classes %>%
  select(position, gene, q_cnn_class, wt_class) %>%
  unnest(q_cnn_class) %>%
  pivot_longer(c(proline:polar), names_to = "class_names", values_to = "class_freqs") %>%
  mutate(freqs = sapply(class_freqs, function(x) x[[1]])) 

count_mismatches <- count_mismatches[order(count_mismatches$gene,count_mismatches$position, -count_mismatches$freqs),]
count_mismatches <- count_mismatches %>%  
  group_by(gene, position) %>%
  summarise(wt_class = first(wt_class), cnn_class = first(class_names)) %>%
  mutate(match = (wt_class == cnn_class))
  
count_mismatches2 <- count_mismatches %>%
  filter(match == FALSE) %>%
  group_by(gene) %>%
  count() %>%
  left_join(
    (count_mismatches %>%
      group_by(gene) %>%
      count()
  ),by="gene")
  
count_mismatches2 <- count_mismatches2 %>%
  mutate(proportion = n.x / n.y)

count_mismatches2 %>%
  ggplot(aes(x = "", y = proportion)) +
  geom_violin(fill = "grey75") + geom_sina() +
  ggtitle(label = "Class mismatches (cnn max freq class vs cnn wt class") +
  theme_cowplot() +
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) +
  xlab("") +
  ylab("proportion mismatch") 

#now find mismatches for amino acids (cnn max vs cnn wt)
aa_mismatches <- classes %>%
  select(gene, position, wt_aa, q_cnn)
aa_mismatches <- aa_mismatches %>%
  unnest(q_cnn) %>%
  pivot_longer(c(q_H:q_C), names_to = "aa", values_to = "freqs") %>%
  mutate(freqs = sapply(freqs, function(x) x[[1]])) 

aa_mismatches <- aa_mismatches[order(aa_mismatches$gene, aa_mismatches$position, -aa_mismatches$freqs),]
aa_mismatches <- aa_mismatches %>%  
  group_by(gene, position) %>%
  summarise(wt_aa = first(wt_aa), aa = first(aa)) #%>%
  #mutate(match = (wt_class == cnn_class))
aa_mismatches <- aa_mismatches %>%  
  mutate(wt_aa_fixed = map(wt_aa, function(q){
    if (q == 'HIS') return ('q_H') 
    if (q == 'GLU') return ('q_E') 
    if (q == 'ASP') return ('q_D') 
    if (q == 'ARG') return ('q_R') 
    if (q == 'LYS') return ('q_K') 
    if (q == 'SER') return ('q_S') 
    if (q == 'THR') return ('q_T') 
    if (q == 'ASN') return ('q_N') 
    if (q == 'GLN') return ('q_Q') 
    if (q == 'ALA') return ('q_A') 
    if (q == 'VAL') return ('q_V') 
    if (q == 'LEU') return ('q_L') 
    if (q == 'ILE') return ('q_I') 
    if (q == 'MET') return ('q_M') 
    if (q == 'PHE') return ('q_F') 
    if (q == 'TYR') return ('q_Y') 
    if (q == 'TRP') return ('q_W') 
    if (q == 'PRO') return ('q_P') 
    if (q == 'GLY') return ('q_G') 
    if (q == 'CYS') return ('q_C')
  })) %>%
  mutate(match = (wt_aa_fixed == aa))

aa_mismatches2 <- aa_mismatches %>%
  filter(match == FALSE) %>%
  group_by(gene) %>%
  count() %>%
  left_join(
    (aa_mismatches %>%
       group_by(gene) %>%
       count()
    ),by="gene")

aa_mismatches2 <- aa_mismatches2 %>%
  mutate(proportion = n.x / n.y)

count_mismatches2 %>%
  ggplot(aes(x = "", y = proportion)) +
  geom_violin(fill = "grey75") + geom_sina() +
  ggtitle(label = "aa mismatches (cnn max freq aa vs cnn wt aa") +
  theme_cowplot() +
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) +
  xlab("") +
  ylab("proportion mismatch") 


#-----------------------------------------------------------------------------------------------------
# 4) average of differences between neff natural and neff CNN over a 5 residue sliding window

joined_data5 <- joined_data3 %>%
  mutate(diff = abs(n_eff.x - n_eff.y), by_five = floor((position-1)/5)) %>%
  group_by(gene, by_five) %>%
  summarize(mean = mean(diff), position = min(position))

joined_data5 %>%
  #filter(gene == "1b4t") %>%
  ggplot(aes(x = position, y = mean)) +
  geom_line() + #geom_point() +
  ggtitle(label = "1b4t (bin by 5)") +
  #geom_abline(slope = 1, color = "grey65" ) +
  theme_cowplot() +
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) +
  facet_wrap(~gene) +
  xlab("position") +
  ylab("n-eff mean diff (cnn vs natural)") 

#--------------------------------------------------------------------------------------
# comparing to consensus
#---------------------------------------------------------------------------------------
# probability of consensus in cnn data ~ pobability of cnn winner
r1 <- joined_data3 %>%
  select(gene, position, q_cnn, q_natural) %>%
  mutate(max_cnn = map(q_cnn, function(x) max(as.numeric(as.character(unlist(x)))))
  ) %>%
  mutate(max_cnn = sapply(max_cnn, function(x) x[[1]])) %>%
  unnest(q_natural) 

r1 <- r1 %>%
  pivot_longer(c(q_H:q_C), names_to = "aa_nat", values_to = "freq_nat") %>%
  mutate(freq_nat = sapply(freq_nat, function(x) x[[1]])) 

r1 <- r1[order(r1$gene, r1$position, -r1$freq_nat),]

r1 <- r1 %>%
  mutate(rn = row_number())
r1 <- subset(r1, (rn-1) %% 20 == 0 ) 
r1 <- r1 %>%
  mutate(nat_cons_cnn = map2(q_cnn, aa_nat, function(x, y) x[,y,drop=TRUE][[1]])
         ) %>%
  mutate(nat_cons_cnn = sapply(nat_cons_cnn, function(x) x[[1]]))

r1 %>% 
  filter(gene == "2aiu") %>%
  ggplot(aes(x = max_cnn, y = nat_cons_cnn)) +
  geom_point() +
  ggtitle(label = "2aiu") +
  geom_abline(slope = 1, intercept = 0, color = "grey65" ) +
  stat_cor(label.y = 0.95, label.x = 0) +
  stat_regline_equation(label.y = .85, label.x = 0) +
  theme_cowplot() +
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) +
  ylab("cnn prob of consensus from alignment") +
  xlab("cnn max probability") 

cor3 <- r1 %>%
  select(gene, max_cnn, nat_cons_cnn) %>%
  na.omit() %>%
  group_by(gene) %>%
  summarize(cor = cor(max_cnn, nat_cons_cnn))


cor3 %>%
  ggplot(aes(x = "", y = cor)) +
  geom_violin(fill = "grey75") + geom_sina() + 
  ggtitle(label = "prob of natural consensus (CNN) ~ max prob (CNN)") +
  theme_cowplot() +
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) +
  xlab("") +
  ylab("correlation coefficients") 


# freq of consensus in alignment ~ freq of cnn winner in alignment
r2 <- joined_data3 %>%
  select(gene, position, q_cnn, q_natural) %>%
  mutate(max_nat = map(q_natural, function(x) max(as.numeric(as.character(unlist(x)))))
  ) %>%
  mutate(max_nat = sapply(max_nat, function(x) x[[1]])) %>%
  unnest(q_cnn) 

r2 <- r2 %>%
  pivot_longer(c(q_H:q_C), names_to = "aa_cnn", values_to = "freq_cnn") %>%
  mutate(freq_cnn = sapply(freq_cnn, function(x) x[[1]])) 

r2 <- r2[order(r2$gene, r2$position, -r2$freq_cnn),]

r2 <- r2 %>%
  mutate(rn = row_number())
r2 <- subset(r2, (rn-1) %% 20 == 0 ) 
r2 <- r2 %>%
  mutate(nat_freq_cnn_max = map2(q_natural, aa_cnn, function(x, y) x[,y,drop=TRUE][[1]])
  ) %>%
  mutate(nat_freq_cnn_max = sapply(nat_freq_cnn_max, function(x) x[[1]]))

r2 %>% 
  filter(gene == "1b4t") %>%
  ggplot(aes(x = nat_freq_cnn_max, y =  max_nat)) +
  geom_point() +
  ggtitle(label = "1b4t") +
  geom_abline(slope = 1, intercept = 0, color = "grey65" ) +
  stat_cor(label.y = 0.95, label.x = 0.25) +
  stat_regline_equation(label.y = .85, label.x = 0.25) +
  theme_cowplot() +
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) +
  ylab("max freq in natural alignment") +
  xlab("freq of cnn max in natural alignment") 

cor4 <- r2 %>%
  select(gene, max_cnn, nat_cons_cnn) %>%
  na.omit() %>%
  group_by(gene) %>%
  summarize(cor = cor(max_cnn, nat_cons_cnn))


cor3 %>%
  ggplot(aes(x = "", y = cor)) +
  geom_violin(fill = "grey75") + geom_sina() + 
  ggtitle(label = "prob of natural consensus (CNN) ~ max prob (CNN)") +
  theme_cowplot() +
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) +
  xlab("") +
  ylab("correlation coefficients") 


 
  

#==============================================================================
#Basic stats (% match across all genes)
#===============================================================================

# check if predicted aa is the one found in the wt structure
matches_1 <- joined_data %>%
  pivot_wider(names_from = group, values_from = c(aa,freq)) %>%
  mutate(match = aa_predicted == aa_wt) %>%
  select(gene, position, match)

# selecting the data entries where the predicted amino acid matches the 
summary_stats_1 <- matches_1 %>%
  group_by(gene, match) %>%
  summarise(count = n()) %>%
  mutate(freq = count / sum(count)) %>%
  filter(match == TRUE)

summary_stats_1 %>%
  ggplot(aes(x = "", y = freq)) +
  geom_violin(fill = "grey75") + geom_sina() +
  ggtitle(label = "CNN Accuracy by Protein") +
  theme_cowplot() + 
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) +
  xlab("") +
  ylab("Percent Accuracy \n (predicted = wt)")   

summary_stats_1 %>%
  summarise(mean = mean(freq))
# mean = 0.713

# check if predicted aa is the one found in the wt structure
matches_2 <- joined_data %>%
  pivot_wider(names_from = group, values_from = c(aa,freq)) %>%
  mutate(match = aa_predicted == aa_natural_max) %>%
  select(gene, position, match)

summary_stats_2 <- matches_2 %>%
  group_by(gene, match) %>%
  summarise(count = n()) %>%
  mutate(freq = count / sum(count)) %>%
  filter(match == TRUE)

summary_stats_2 %>%
  ggplot(aes(x = "", y = freq)) +
  geom_violin(fill = "grey75") + geom_sina() +
  ggtitle(label = "CNN Predictions Compared to \n Alignment Consensus") +
  theme_cowplot() + 
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) +
  xlab("") +
  ylab("Percent Accuracy \n (predicted = consensus)")   

summary_stats_2 %>%
  summarise(mean = mean(freq))
# mean = 0.331





