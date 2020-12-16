# comparing amino acid class distributions

library(tidyverse)
library(yardstick)
library(ggpubr)
library(cowplot)
library(sqldf)
library(sinaplot)
library(dplyr)
library(ggforce)

# reading csv files
natural_var <- read.csv(file="./stats_align_all.csv", header=TRUE, sep=",")
cnn_var <- read.csv(file="./stats_cnn.csv", header=TRUE, sep=",")