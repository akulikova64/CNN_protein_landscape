# check the number of genes in files to track if any genes were lost


natural_data <- read.csv(file = "./data/PSICOV_box_20/output/natural_max_freq_files/natural_max_freq_all.csv", header=TRUE, sep=",")

stats_align_all <- read.csv(file = "./output/output_PSICOV/stats_align_all.csv", header=TRUE, sep=",")
cnn_wt_max_freq <- read.csv(file = "./data/PSICOV_box_20/output/cnn_wt_max_freq.csv", header=TRUE, sep=",")


count_stats <- stats_align_all %>%
  select(gene) %>%
  group_by(gene) %>%
  count() #149 genes

count_cnn <- cnn_wt_max_freq %>%
  select(gene) %>%
  group_by(gene) %>%
  count() %>%
  select(-n) #143 genes

count_cnn
