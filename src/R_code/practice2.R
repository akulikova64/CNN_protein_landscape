library(tidyverse)
library(cowplot)

# paractice with map()

calc <- function(x) {
  cat = c("a", "b")
  dog = c("c", "d")
  
  if (x %in% cat) {
    return("cat")
  } 
  if (x %in% dog) {
    return("dog")
  }
}

table <- tibble(col1 = c("a", "b", "c", "d"))

new_table <- table %>%
  mutate(col2 = map_chr(col1, calc))
