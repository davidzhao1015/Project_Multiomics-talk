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


# Input data: Multi-omics citation of HMDB, Drugbank and MetaboloAnalyst 
head(citation_multiomics) 

# Input data: entire citation of HMDB, Drugbank and MetaboloAnalyst 
head(merged_df_unique)

unique(merged_df_unique$database)


# Input data: entire multiomics papers 
head(total_article) 



#----------------Word frequency in entire citation of HMDB and drugbank------------

multiomics_pmcid <- total_article$PMID 

merged_df_unique_tag <- merged_df_unique %>% 
    mutate(is_multiomics = case_when(PMID %in% multiomics_pmcid ~ "yes",
                                     !(PMID %in% multiomics_pmcid) ~ "no")) %>%
    mutate(after_2018 = case_when(Publication.Year >= 2018 ~ "yes",
                                  Publication.Year < 2018 ~ "no"))

citations_by_db <- merged_df_unique_tag %>% 
    filter(database %in% c("human"))  

head(citations_by_db)

full_title_year2 <- citations_by_db %>% 
    select(Title, database, is_multiomics, Publication.Year) %>% 
    filter(is_multiomics == "yes") 

# Tokenize 
tidy_title_year2 <- full_title_year2 %>%
    unnest_tokens(word, Title)

# Exclude stop words 
data("stop_words") 

tidy_title_year2 <- tidy_title_year2 %>% 
    anti_join(stop_words) 

tidy_title_year2 %>% count(word, sort = TRUE) 

# Calculate word frequency
frequency2 <- tidy_title_year2 %>%
    group_by(Publication.Year) %>%
    count(word, sort = TRUE)
    mutate(prop = n / sum()) %>%
    select(-n) %>%
    ungroup()
    
frequency2 %>% 
    filter(word %in% list_diff_2018) %>% 
    ggplot(aes(x = Publication.Year, y = n)) +
    geom_col() +
    geom_vline(xintercept = 2018, color = "red") + 
    facet_wrap(~word)
    
list_before2018 <- frequency2 %>%
    slice_max(prop, n = 15) %>%
    filter(after_2018 == "no") %>%
    pull(word)

list_after2018 <- frequency2 %>%
    slice_max(prop, n = 15) %>%
    filter(after_2018 == "yes") %>%
    pull(word)

list_after2018[!(list_after2018 %in% list_before2018)]


frequency3 <- frequency2 %>% 
    pivot_wider(id_cols = word, names_from = after_2018, values_from = prop, values_fill = 0)

ggplot(frequency3, 
       aes(x = yes, y = no), color = abs(yes - no)) +
    geom_abline(color = "gray40", lty = 2) +
    geom_jitter(alpha = 0.1, size = 2.5, width = 0.3, height = 0.3) + 
    geom_text(aes(label = word), check_overlap = TRUE, vjust = 1.5) +
    scale_x_log10(labels = percent_format()) +
    scale_y_log10(labels = percent_format()) +
    scale_color_gradient(limits = c(0, 0.0000001), high = "red", low = "gray50") +
    theme(legend.position = "none") +
    theme_classic()

cor.test(data = frequency3, ~ yes + no) 



# Cancer - HMDB and drugbank users both favoriate 
top_drugbank_wordfreq <- frequency3 %>% arrange(desc(drugbank)) %>% pull(word)
list_drug <- top_drugbank_wordfreq[1:15]

top_human_wordfreq <- frequency3 %>% arrange(desc(human)) %>% pull(word)
list_human <- top_human_wordfreq[1:15]

list_human[list_human %in% list_drug] 


# ti-idf in citations based on multiomics: HMDB vs drugbank 
human_title_year2 <- citations_by_db %>% 
    select(Title, database, is_multiomics) %>% 
    filter(is_multiomics == "yes")

head(human_title_year2)

 human_title_words <- human_title_year2 %>%
    unnest_tokens(word, Title) %>% 
    count(database, word, sort = TRUE) 

human_total_words <- human_title_words %>% 
    group_by(database) %>% 
    summarize(total = sum(n)) 

human_title_words <- left_join(human_title_words, human_total_words)

human_title_tf_idf <- human_title_words %>% 
    bind_tf_idf(word, database, n) %>% 
    select(-total) %>% 
    arrange(desc(tf_idf)) 

human_title_tf_idf %>% 
    group_by(database) %>% 
    slice_max(tf_idf, n = 10) %>% 
    ungroup() %>%
    ggplot(aes(tf_idf, fct_reorder(word, tf_idf), fill = database)) +
    geom_col(show.legend = FALSE) + 
    facet_wrap(~database, ncol = 2, scales = "free") +
    labs(x = "tf-idf", y = NULL) 



# ti-idf in citations of drugbank: non-multiomics vs multiomics 

drug_title_year2 <- citations_by_db %>% 
    select(Title, database, is_multiomics) %>% 
    filter(database == "drugbank")

head(drug_title_year2)

drug_title_words <- drug_title_year2 %>%
    unnest_tokens(word, Title) %>% 
    count(is_multiomics, word, sort = TRUE) 

drug_total_words <- drug_title_words %>% 
    group_by(is_multiomics) %>% 
    summarize(total = sum(n)) 

drug_title_words <- left_join(drug_title_words, drug_total_words)

drug_title_tf_idf <- drug_title_words %>% 
    bind_tf_idf(word, is_multiomics, n) %>% 
    select(-total) %>% 
    arrange(desc(tf_idf)) 

drug_title_tf_idf %>% 
    group_by(is_multiomics) %>% 
    slice_max(tf_idf, n = 7) %>% 
    ungroup() %>%
    ggplot(aes(tf_idf, fct_reorder(word, tf_idf), fill = is_multiomics)) +
    geom_col(show.legend = FALSE) + 
    facet_wrap(~is_multiomics, ncol = 2, scales = "free") +
    labs(x = "tf-idf", y = NULL) 


# Chi-sqaure test 
table(citations_by_db$database, citations_by_db$is_multiomics)

chisq <- chisq.test(citations_by_db$database, citations_by_db$is_multiomics)

print(chisq) # Significant associated 

round(chisq$residuals, 3) # HMDB is positively associated with is_multiomics 

library(corrplot)

corrplot(chisq$residuals, is.cor = FALSE)







