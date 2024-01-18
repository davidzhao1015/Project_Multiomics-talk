library(tidyverse)

# Input data: Multi-omics citation of HMDB, Drugbank and MetaboloAnalyst 
head(citation_multiomics) 

# Input data: entire citation of HMDB, Drugbank and MetaboloAnalyst 
head(merged_df_unique)

unique(merged_df_unique$database)

# Input data: entire multiomics papers 
head(total_article) 



# Stacked bar plot x represents 3 tools, y total and multiomics- citations 

multiomics_pmcid <- total_article$PMID 

merged_df_unique_tag <- merged_df_unique %>% 
    mutate(is_multiomics = case_when(PMID %in% multiomics_pmcid ~ "yes",
                                     !(PMID %in% multiomics_pmcid) ~ "no"))

citations_by_db <- merged_df_unique_tag %>% 
    filter(database %in% c("metaboanalyst", "omicsnet", "paintsomics", "xcms", "workflow", "omicsanalyst", "metabolights")) %>% 
    filter(Publication.Year >= 0) %>% 
    group_by(database) %>% 
    count(is_multiomics) %>% 
    ungroup()

head(citations_by_db)

# Count bar 
citations_by_db %>% 
    ggplot(aes(x = database, y = n, fill = (is_multiomics == "no"))) +
    geom_bar(position = "stack", stat = "identity") +
    coord_flip()
    

# Proportion bar 
ggplot(citations_by_db, aes(x = database, y = n, fill = (is_multiomics == "no"))) +
    geom_bar(position = "fill", stat = "identity") +
    coord_flip()


# Bubble chart 

head(merged_df_unique_tag)

merged_df_percent_omics <- merged_df_unique_tag %>% 
    filter(database %in% c("metaboanalyst", "omicsnet", "paintsomics", "xcms", "workflow", "omicsanalyst", "metabolights")) %>% 
    group_by(database, Publication.Year) %>% 
    count(is_multiomics) %>% 
    ungroup() %>% 
    pivot_wider(names_from = is_multiomics, values_from = n, values_fill = 0) %>% 
    mutate(total = yes + no) %>% 
    mutate(percent_omics = yes/(yes + no)) 

ggplot(merged_df_percent_omics, 
       aes(x = Publication.Year, 
           y = percent_omics, 
           size = total, color = database)) +
           geom_point(alpha = .8)



#--------------------Word Frequency--------------------------


sub_by_db <- merged_df_unique_tag %>% 
    filter(database %in% c("metaboanalyst", "omicsnet", "paintsomics", "xcms", "workflow", "omicsanalyst", "metabolights"))


# Tokenize 
tidy_sub_by_db <- sub_by_db %>%
    unnest_tokens(word, Title)

data("stop_words") # Exclude stop words  

tidy_sub_by_db <- tidy_sub_by_db %>% 
    anti_join(stop_words) 

head(tidy_sub_by_db)

# Calculate word frequency
frequency_sub_db <- tidy_sub_by_db %>% 
    group_by(database) %>% 
    count(word, sort = TRUE) %>% 
    mutate(prop = n / sum(n)) %>%
    select(-n)

head(frequency_sub_db)

frequency_sub_db2 <- frequency_sub_db %>% 
    pivot_wider(names_from = database, values_from = prop, values_fill = 0) 

freq_matrix <- frequency_sub_db2 %>% 
    column_to_rownames("word") %>% 
    as.matrix()

res <- cor(freq_matrix)

res2 <- rcorr(freq_matrix)

corrplot(res, type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 45) 

heatmap(x = res, symm = TRUE)

# corrplot(res2$r, p.mat = res2$P)



# Line plot - x is years (and dates) while y is citation numbers 
merged_df_unique_tag %>% 
    filter(database %in% c("human", "drugbank", "metaboanalyst")) %>% 
    filter(Publication.Year >= 0 & Publication.Year < 2023) %>% 
    group_by(database, Publication.Year, is_multiomics) %>%
    count(is_multiomics) %>% 
    ggplot(aes(x = Publication.Year, y = n, color =  database))+
    geom_line(aes(linetype = is_multiomics)) 

# Line plot - x is year while y is multiomics citation numbers 
merged_df_unique_tag %>% 
    filter(database %in% c("human")) %>% 
    filter(Publication.Year >= 0 & Publication.Year < 2023) %>% 
    filter(is_multiomics == "yes") %>% 
    group_by(database, Publication.Year) %>%
    count(is_multiomics) %>% 
    ggplot(aes(x = Publication.Year, y = n, fill =  database))+
    geom_line(aes(color = database)) 


# Line plot - y is cumulative sum while x is creation.date 
merged_df_unique_tag %>% 
    filter(database %in% c("human", "drugbank")) %>% 
    filter(Publication.Year >= 0 & Publication.Year < 2023) %>% 
    filter(is_multiomics == "yes") %>% 
    group_by(database, Publication.Year) %>% 
    count(is_multiomics) %>% 
    ungroup() %>% 
    group_by(database) %>% 
    mutate(cumulative_sum = cumsum(n)) %>%  
    ggplot(aes(x = Publication.Year, y = cumulative_sum, group =  database))+
    geom_line(aes(color = database)) 


#------------------Analyze tf-idf----------------------- 

citation_title_year <- merged_df_unique_tag %>% 
    filter(database %in% c("human")) %>% 
    filter(Publication.Year >= 0 & Publication.Year < 2023) %>% 
    filter(is_multiomics == "yes") %>% 
    mutate(after_2018 = case_when(Publication.Year >= 2018 ~ "yes",
                                  Publication.Year < 2018 ~ "no")) %>% 
    select(Title, after_2018) 

citation_title_words <- citation_title_year %>%
    unnest_tokens(word, Title) %>% 
    count(after_2018, word, sort = TRUE) 

citation_total_words <- citation_title_words %>% 
    group_by(after_2018) %>% 
    summarize(total = sum(n)) 

citation_title_words <- left_join(citation_title_words, citation_total_words)

citation_title_tf_idf <- citation_title_words %>% 
    bind_tf_idf(word, after_2018, n) %>% 
    select(-total) %>% 
    arrange(desc(tf_idf)) 

head(citation_title_tf_idf)

citation_title_tf_idf %>%
    group_by(after_2018) %>%
    slice_max(tf_idf, n = 10) %>%
    ungroup() %>%
    ggplot(aes(tf_idf, fct_reorder(word, tf_idf), fill = after_2018)) +
    geom_col(show.legend = FALSE) +
    facet_wrap(~after_2018, ncol = 2, scales = "free") +
    labs(x = "tf-idf", y = NULL)

list_diff_2018 <- citation_title_tf_idf %>%
    group_by(after_2018) %>%
    slice_max(tf_idf, n = 10) %>%
    pull(word)



#-------------------Tokenizing by n-gram--------------------------

citation_title_bigrams <- citation_title_year %>% 
    unnest_tokens(bigram, Title, token = "ngrams", n = 2) %>% 
    filter(!is.na(bigram)) 

head(citation_title_bigrams)



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

count_bigrams <- count_bigrams(citation_title_year) 

count_bigrams %>% 
    filter(n > 10,
           !str_detect(word1, "\\d"),
           !str_detect(word2, "\\d")) %>% 
    visualize_bigrams()


#----------------------MetaboloAnalyst------------------------------------------

citation_title_year_ma <- merged_df_unique_tag %>% 
    filter(database %in% c("metaboanalyst", "omicsnet", "paintsomics")) %>% 
    filter(is_multiomics == "yes") %>% 
    select(Title, database) 

# head(citation_title_year_ma)

citation_title_words_ma <- citation_title_year_ma %>%
    unnest_tokens(word, Title) %>% 
    count(database, word, sort = TRUE) 

citation_total_words_ma <- citation_title_words_ma %>% 
    group_by(database) %>% 
    summarize(total = sum(n)) 

citation_title_words_ma <- left_join(citation_total_words_ma, citation_title_words_ma)

citation_title_tf_idf_ma <- citation_title_words_ma %>% 
    bind_tf_idf(word, database, n) %>% 
    select(-total) %>% 
    arrange(desc(tf_idf)) 

citation_title_tf_idf_ma %>%
    group_by(database) %>%
    slice_max(tf_idf, n = 10) %>%
    ungroup() %>%
    ggplot(aes(tf_idf, fct_reorder(word, tf_idf), fill = database)) +
    geom_col(show.legend = FALSE) +
    facet_wrap(~database, ncol = 3, scales = "free") +
    labs(x = "tf-idf", y = NULL)



# n-gram for metaboloAnalyst 
citation_title_bigrams_ma <- citation_title_year_ma %>% 
    unnest_tokens(bigram, Title, token = "ngrams", n = 2) %>% 
    filter(!is.na(bigram)) 

head(citation_title_bigrams_ma)

ma_separated <- citation_title_bigrams_ma %>% 
    separate(bigram, c("word1", "word2"), sep = " ") 

ma_filtered <- ma_separated %>% 
    filter(!word1 %in% stop_words$word) %>% 
    filter(!word2 %in% stop_words$word) 

bigram_counts_ma <- ma_filtered %>% 
    count(word1, word2, sort = TRUE) 

ma_filtered %>% 
    filter(word2 == "microbiota") %>% 
    count(database, word1, sort = T)

ma_united <- ma_filtered %>% 
    unite(bigram, word1, word2, sep = " ")

bigram_tf_idf_ma <- ma_united %>% 
    count(database, bigram) %>% 
    bind_tf_idf(bigram, database, n) %>% 
    arrange(desc(tf_idf)) 

head(bigram_tf_idf_ma)

bigram_tf_idf_ma %>% 
    group_by(database) %>% 
    slice_max(tf_idf, n = 7) %>% 
    ungroup() %>%
    ggplot(aes(tf_idf, bigram, fill = database)) +
    geom_col(show.legend = FALSE) + 
    facet_wrap(~database, ncol = 3, scales = "free") +
    labs(x = "tf-idf", y = NULL) 







