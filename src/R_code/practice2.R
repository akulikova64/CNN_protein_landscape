# paractice with map()

calc <- function(x) {
  return("done")
}

table <- tibble(a = c(1, 2, 3), b = c(4, 5, 6))


new_table <- table %>%
  mutate(new = map(a, calc))
