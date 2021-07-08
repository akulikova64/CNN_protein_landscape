library(tidyverse)

#counting the consensus ties

tie_data <- read.csv(file = "./data/PSICOV_box_20/output/consensus_ties.csv", header=TRUE, sep=",")

counts <- tie_data %>%
  mutate(ties = ifelse(tie_count > 0, TRUE, FALSE)) %>%
  group_by(ties) %>%
  count()

ties_only <- tie_data %>%
  filter(tie_count > 0)
