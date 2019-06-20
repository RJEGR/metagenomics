library(devtools)
install_github("juliasilge/tidytext")

# AND APPLY inverse document frequency (idf) ANALYSIS

library(dplyr)

library(janeaustenr)
original_books <- austen_books() %>%
  group_by(book) %>%
  mutate(line = row_number()) %>%
  ungroup()

library(tidytext)
tidy_books <- original_books %>%
  unnest_tokens(word, text)

tidy_books <- tidy_books %>%
  anti_join(get_stopwords())

tidy_books %>%
  count(word, sort = TRUE)

library(tidyr)
get_sentiments("bing")

janeaustensentiment <- tidy_books %>%
  inner_join(get_sentiments("bing"), by = "word") %>% 
  count(book, index = line %/% 80, sentiment) %>% 
  spread(sentiment, n, fill = 0) %>% 
  mutate(sentiment = positive - negative)

library(ggplot2)
ggplot(janeaustensentiment, aes(index, sentiment, fill = book)) +
  geom_bar(stat = "identity", show.legend = FALSE) +
  facet_wrap(~book, ncol = 2, scales = "free_x")
