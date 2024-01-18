# Load library 
library(tidytext) 
library(janeaustenr)
library(scales)

library(stringr)
library(wordcloud)
library(RColorBrewer)
library(wordcloud2) 
library(tm)
library(tidyr)
library(tidyverse)
library(factoextra)
library(FactoMineR)
library(ggfortify)  # Plot PCA 
library(wrangle)
library(fmsb) # Radar chart 
library(packcircles)
library("htmlwidgets") 
library(webshot)
webshot::install_phantomjs(force = TRUE)


# Set the working directory 
setwd("~/Project_Multiomics-talk/")  


# Multi-omics paper, total 
full_article <- read.csv("./all_articles_multiomics.csv") 

full_title_year <- full_article %>% select(Title, Publication.Year, PMID) 

# Tokenize 
tidy_title_year <- full_title_year %>%
    unnest_tokens(word, Title)

# tmic_title <- merged_df_unique %>% 
#     filter(database == "human") %>% unnest_tokens(word, Title) 

data("stop_words") # Exclude stop words  

tidy_title_year <- tidy_title_year %>% 
    anti_join(stop_words) 

# tidy_title_year %>% count(word, sort = TRUE) 

# tmic_title_year <- tmic_title %>% 
#     anti_join(stop_words)

# Calculate word frequency
frequency <- tidy_title_year %>% 
    group_by(Publication.Year) %>% 
    count(word, sort = TRUE) %>% 
    mutate(prop = n / sum(n)) %>%
    select(-n)


freq_tmic <- tmic_title_year %>% 
    group_by(Publication.Year) %>% 
    count(word, sort = TRUE) %>% 
    mutate(prop = n / sum(n)) %>%
    select(-n)

head(freq_tmic)


frequency_2023 <- frequency %>%
    filter(Publication.Year == 2023) %>% 
    rename(yr_2023 = prop)

head(frequency_2023)

frequency_2022 <- frequency %>% 
    filter(Publication.Year == 2022) %>% 
    rename(yr_2022 = prop)

freq_2023_2022 <- frequency_2023 %>% 
    full_join(frequency_2022, by = join_by(word)) 



freq_2023_vs_rest <- frequency_2023 %>% 
    full_join(frequency, by = join_by(word)) %>% 
    select(-1) %>% 
    rename(prop_2023 = yr_2023, publication.year = Publication.Year.y) %>% 
    filter(publication.year != 2023)


freq_2023_tmic <- frequency_2023 %>% 
    full_join(freq_tmic, by = join_by(word)) %>% 
    select(-1) %>% 
    rename(publication.year = Publication.Year.y, tmic_prop = prop) 

head(freq_2023_tmic)



# Scatter plot 
ggplot(freq_2023_tmic, 
       aes(x = tmic_prop, y = yr_2023),
       color = abs(yr_2023 - tmic_prop)) +
    geom_abline(color = "gray40", lty = 2) +
    geom_jitter(alpha = 0.1, size = 2.5, width = 0.3, height = 0.3) + 
    geom_text(aes(label = word), check_overlap = TRUE, vjust = 1.5) +
    scale_x_log10(labels = percent_format()) +
    scale_y_log10(labels = percent_format()) +
    scale_color_gradient(limits = c(0, 1), low = "darkslategray4", high = "gray50") +
    facet_wrap(~publication.year) +
    theme(legend.position = "none")
    
# Pearson correlation 
cor.test(data = freq_2023_tmic[freq_2023_tmic$publication.year == 2023, ],
         ~ tmic_prop + yr_2023) 



#------------------Analyze tf-idf-----------------------

title_words <- full_title_year %>%
    unnest_tokens(word, Title) %>% 
    count(Publication.Year, word, sort = TRUE) 

total_words <- title_words %>% 
    group_by(Publication.Year) %>% 
    summarize(total = sum(n)) 

title_words <- left_join(title_words, total_words)

# head(title_words)

title_tf_idf <- title_words %>% 
    bind_tf_idf(word, Publication.Year, n) %>% 
    select(-total) %>% 
    arrange(desc(tf_idf)) 


title_tf_idf %>% 
    group_by(Publication.Year) %>% 
    slice_max(tf_idf, n = 15) %>% 
    ungroup() %>%
    ggplot(aes(tf_idf, fct_reorder(word, tf_idf), fill = Publication.Year)) +
    geom_col(show.legend = FALSE) + 
    facet_wrap(~Publication.Year, ncol = 3, scales = "free") +
    labs(x = "tf-idf", y = NULL) 




#-------------------Tokenizing by n-gram--------------------------

title_bigrams <- full_title_year %>%
    unnest_tokens(bigram, Title, token = "ngrams", n = 2) %>% 
    filter(!is.na(bigram)) 

head(title_bigrams)

full_title_2023 <- full_title_year %>% filter(Publication.Year == 2023)



# Visualize a network of bigrams with ggraph 

library(igraph)
library(ggraph)

# Custom funcitons 
count_bigrams <- function(dataset) {
    dataset %>%
        unnest_tokens(bigram, Title, token = "ngrams", n = 2) %>%
        separate(bigram, c("word1", "word2"), sep = " ") %>%
        filter(!word1 %in% stop_words$word,
               !word2 %in% stop_words$word) %>%
        count(word1, word2, sort = TRUE)
}

visualize_bigrams <- function(bigrams) {
    set.seed(2016)
    a <- grid::arrow(type = "closed", length = unit(.15, "inches"))
    
    bigrams %>%
        graph_from_data_frame() %>%
        ggraph(layout = "fr") +
        geom_edge_link(aes(edge_alpha = n), show.legend = FALSE, arrow = a) +
        geom_node_point(color = "lightblue", size = 5) +
        geom_node_text(aes(label = name), vjust = 1, hjust = 1) +
        theme_void()
}

count_bigrams_title <- count_bigrams(full_title_2023)

count_bigrams_title %>% 
    filter(n > 10,
           !str_detect(word1, "\\d"),
           !str_detect(word2, "\\d")) %>% 
    visualize_bigrams()




#----------------------------------MetaboloAnalyst-----------------------------------------

table(merged_df_unique$database, merged_df_unique$Publication.Year)

head(citation_multiomics) 








#----------------Topic Modeling Entire Multi-omics Papers-----------------------

# Ref: https://www.tidytextmining.com/topicmodeling#library-heist 

# Get title-word counts 
pmid_word_counts <- tidy_title_year %>% 
    count(PMID, word, sort = TRUE)

# Convert into DocumentTermMatrix object 
title_dtm <- pmid_word_counts %>% 
    cast_dtm(PMID, word, n)

library(topicmodels)

title_lda <- LDA(title_dtm, k = 4, control = list(seed = 1234))  

title_topics <- tidy(title_lda, matrix = "beta")
title_topics

# Find top terms within each topic 
top_terms <- title_topics %>% 
    group_by(topic) %>%
    slice_max(beta, n = 15) %>% 
    ungroup() %>% 
    arrange(topic, -beta)


# Plot top terms associated with each topic 
top_terms %>% 
    mutate(term = reorder_within(term, beta, topic)) %>% 
    ggplot(aes(beta, term, fill = factor(topic))) + 
    geom_col(show.legend = FALSE) + 
    facet_wrap(~ topic, scales = "free") +
    scale_y_reordered()


# Classify titles into 3 categories 
title_gamma <- tidy(title_lda, matrix = "gamma")
title_gamma 











