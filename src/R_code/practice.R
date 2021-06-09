table <- tibble(x = c(TRUE, TRUE, TRUE, TRUE, FALSE), y = c(FALSE, TRUE, TRUE, TRUE, TRUE))

#when x is TRUE, what is the false rate for y? (should be equal to 3)

data <- table %>%
  filter(x == TRUE) %>%
  summarise(freq = sum(!y, na.rm = TRUE)/sum(!is.na(y)))


#new practice problem:

library(stats)

table <- tibble(x = c("A", "A", "A", "B", "B"), p_value = c(0.5, 0.2, 0.8, 0.5, 0.4))

new_p_values <- table %>%
  nest(data = -c(x)) %>%
  mutate(
    new_p = map(data, ~p.adjust(.x$p_value, method = "fdr", n = length(.x$p_value)))
  ) %>%
  unnest(cols = c(data, new_p))

