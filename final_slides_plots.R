# This R codes is for visualization of my presentation slides about multi-omics and TMIC's bioinformatics tools 

# References:
# - R graphics cookbook
# - ggplot2: Elegant graphics for data analysis (3e) 
# - Text mining with R: A tidy approach 

# Load library 
library(tidytext) 
library(janeaustenr)
library(corrplot)
library("Hmisc")
library(ggrepel)
library(stringr)
library(wordcloud)
library(RColorBrewer)
library(wordcloud2) 
library(gganimate)
library(gifski)
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



#===========================================================================

# Bar plot illustrates multi-omics publications over time 

#===========================================================================

summary_total_article <- total_article %>% 
    group_by(Publication.Year) %>% 
    count(count = n()) %>% 
    ungroup() %>% 
    mutate(cum_count = cumsum(n))

# TODO 
# Axis label 
# 4000 - 4K 
# Bar color 
# Bar spacing and width 

summary_total_article %>%
    filter(Publication.Year < 2024) %>%
    ggplot(aes(x = Publication.Year, y = cum_count)) +
    geom_col(aes(fill = (2020 <= Publication.Year & Publication.Year < 2023)), 
             width = 0.7) +
    geom_smooth(method = "loess", se = FALSE, colour = "grey50", linewidth = 0.5) + 
    scale_x_continuous(breaks = seq(2009, 2023, 1)) +
    scale_fill_manual(values = c('grey', 'purple')) + 
    labs(x = NULL, y = 'Papers on PubMed, cumulative',
         caption = 'The plot depicts the rising trend of multi-omics publications in the PubMed database from 2019 to 2023.') +
    theme(legend.position = "none", 
          panel.background = element_rect(fill = 'white'),
          axis.line = element_line(colour = "grey50"),
          plot.caption.position = "panel",
          plot.caption = element_text(colour = "grey50"))




#===========================================================================

# Word Cloud based on entire publications in PubMed during 2009 and 2024 

#===========================================================================

get_word_frequnecy <- function(input, years) {
    
    total_article0 <- input
    
    total_article_2022 <- total_article0 %>% filter(Publication.Year %in% years) 
    title_human_db_all <- total_article_2022$Title
    # Create a corpus
    title_human_db_corpus_all <- tm::Corpus(VectorSource(title_human_db_all))
    # Clean text
    title_human_db_corpus_all <- title_human_db_corpus_all %>%
        tm_map(stripWhitespace)
    title_human_db_corpus_all <- tm_map(title_human_db_corpus_all,
                                        content_transformer(tolower))
    # Create a document-term-matrix
    doc_term_matrix_all <- TermDocumentMatrix(title_human_db_corpus_all)
    matrix_all <- as.matrix(doc_term_matrix_all)
    words_all <- sort(rowSums(matrix_all), decreasing = TRUE)
    df_all <- data.frame(word = names(words_all), freq = words_all)
    
    return(df_all)
}

years <- seq(2009, 2024, 1) # Customize the year range 

total_w_freq <- get_word_frequnecy(total_article, years) 

data("stop_words")

total_w_freq <- total_w_freq %>% 
    anti_join(stop_words) %>% 
    filter(word != c("the","for","and","with"))

set.seed(1234)
wordcloud(words = total_w_freq$word,
          freq = total_w_freq$freq,
          fixed.asp = TRUE,
          min.freq = 10,
          max.words = 100,
          random.order = TRUE,
          random.color = FALSE,
          col = hcl.colors(n = length(total_w_freq$word), 
                           palette = "viridis",
                           alpha = 0.9,
                           rev = TRUE))

# my_graph <- wordcloud2(total_w_freq, size= 1.6) 
# saveWidget(my_graph,"tmp.html",selfcontained = F)
# webshot("tmp.html","fig_1.pdf", delay =5, vwidth = 480, vheight=480)  


#===========================================================================

# Bar plot: Ranking of citations of TMIC's databases  

#===========================================================================

merged_df_unique2 <- merged_df_unique[ ,c("PMID", "database")] %>% 
    filter(database %in% c("csf", 
                           "drugbank", 
                           "ecoli", 
                           "exposure", 
                           "fecal", 
                           "human", 
                           "serum", 
                           "urine"))

tmic_db_pmcid <- merged_df_unique2$PMID
multiomics_pmcid <- total_article$PMID 

total_article2 <- total_article %>% 
    left_join(merged_df_unique2, by = "PMID") %>% 
    mutate(citated_tmic = case_when(PMID %in% tmic_db_pmcid ~ "yes",
                                    !(PMID %in% tmic_db_pmcid) ~ "no")) 

sum(total_article2$citated_tmic == "yes")/ length(total_article2$citated_tmi)

merged_df_unique_tag <- merged_df_unique %>% 
    mutate(is_multiomics = case_when(PMID %in% multiomics_pmcid ~ "yes",
                                     !(PMID %in% multiomics_pmcid) ~ "no"))

citations_by_db_all <- merged_df_unique_tag %>% 
    filter(database %in% c("csf", 
                           "drugbank", 
                           "ecoli", 
                           "exposure", 
                           "fecal", 
                           "human", 
                           "serum", 
                           "urine")) %>% 
    filter(Publication.Year >= 0 & Publication.Year < 2024) %>% 
    group_by(database) %>% 
    count(is_multiomics) %>% 
    ungroup()

head(citations_by_db_all)

# TODO
# Bar width and spacing 
# Bar color 
# Font size 
# Axis line, top 
# Reduce space between axis label and ticks 

cols <- c(rep('grey', 6), '#FF00F7', '#1E98FD')  

citations_by_db_all %>% 
    filter(is_multiomics == "yes") %>% 
    mutate(database = fct_reorder(database, n)) %>% 
    ggplot(aes(x = database, y = n, fill = database)) +
    geom_bar(stat = "identity", 
             alpha = .6, 
             width = .7) +
    geom_text(aes(label = n), vjust = 0.3, hjust = -0.3, colour = "grey50") + 
    scale_x_discrete(labels = c("human" = "HMDB",
                                "drugbank" = "DrugBank",
                                "serum" = "Serum Metabolome DB",
                                "exposure" = "Toxic Exposome DB",
                                "urine" = "Urine Metabolome DB",
                                "fecal" = "Fecal Metabolome DB",
                                "ecoli" = "E.coli Metabolome DB",
                                "csf" = "CSF Metabolome DB")) +
    theme(panel.background = element_rect(fill = 'white'),
          legend.position = 'none',
          axis.ticks.y  = element_blank(), 
          axis.line.x  = element_line(color = 'grey50'),
          plot.margin = margin(t = 10, r = 25, b = 10, l = 25)) + 
    scale_fill_manual(values = cols) + 
    xlab(NULL) + 
    ylab('Citations by Multi-Omics Papers on PubMed (2015-2023)') + 
    coord_flip() 
    


#===========================================================================

# Scatter plot: Word frequency between HMDB and drugbank   

#===========================================================================

multiomics_pmcid <- total_article$PMID 

merged_df_unique_tag <- merged_df_unique %>% 
    mutate(is_multiomics = case_when(PMID %in% multiomics_pmcid ~ "yes",
                                     !(PMID %in% multiomics_pmcid) ~ "no")) %>%
    mutate(after_2018 = case_when(Publication.Year >= 2018 ~ "yes",
                                  Publication.Year < 2018 ~ "no"))

hmdb_drug <- merged_df_unique_tag %>% 
    filter(database %in% c("human", "drugbank"))  

hmdb_drug_omics <- hmdb_drug %>% 
    select(Title, database, is_multiomics, Publication.Year) %>% 
    filter(is_multiomics == "yes") 

# table(hmdb_drug_omics$database) 

hmdb_drug_omics_token <- hmdb_drug_omics %>%
    unnest_tokens(word, Title)

hmdb_drug_omics_token <- hmdb_drug_omics_token %>% 
    anti_join(stop_words) 

hmdb_drug_omics_token2 <- hmdb_drug_omics_token %>%
    group_by(database) %>%
    count(word, sort = TRUE) %>% 
    mutate(prop = n /sum(n)) %>%
    select(-n) %>%
    ungroup()

hmdb_drug_omics_token2 %>%
    group_by(database) %>%
    slice_max(prop, n = 15) %>% 
    ggplot(aes(x = fct_reorder(word, prop), y = prop)) +
    geom_col() +
    coord_flip() +
    facet_wrap(~database, scales = "free")

hmdb_drug_token_wide <- hmdb_drug_omics_token2 %>% 
    pivot_wider(names_from = database, values_from = prop, values_fill = 0) 

cor.test(data = hmdb_drug_token_wide, ~ human + drugbank) 


com_words <- hmdb_drug_token_wide %>% 
    filter(!(word %in% c('omics', 'multi', '2', '19'))) %>% 
    filter(human != 0) %>% 
    filter(drugbank != 0) %>% 
    mutate(abs_diff = abs(human - drugbank)) %>% 
    slice_min(n = 9, order_by = abs_diff) %>% 
    pull(word)

library(scales)
set.seed(123) 

hmdb_drug_token_wide %>% 
    filter(!(word %in% c('omics', 'multi', '2', '19'))) %>% 
    filter(human != 0) %>% 
    filter(drugbank != 0) %>% 
    mutate(diff_db = human - drugbank) %>% 
    mutate(common_word = case_when(word %in% com_words ~ word,
                                   !(word %in% com_words) ~ NA)) %>% 
    ggplot(aes(x = human, y = drugbank)) +
    geom_abline(slope =  1, intercept = 0, color = "gray40", lty = 2) +
    geom_jitter(aes(color = diff_db),
                alpha = 0.5, 
                size = 2.5,
                width = 0.2, 
                height = 0.2) + 
    scale_x_continuous(name = 'HMDB', 
                       limits = c(0.001, 0.04), 
                       trans = "log10", 
                       labels = percent) +
    scale_y_continuous(name = "DrugBank", 
                       limits = c(0.001, 0.04), 
                       trans = "log10", 
                       labels = percent) +
    geom_text_repel(aes(label = word), 
              na.rm = TRUE) +
    scale_colour_gradient2(low = "#FF00F7", 
                           mid = 'green',
                           high = "#1E98FD",
                           midpoint = 0) +
    theme_classic() + 
    theme(legend.position = "none") 



#===========================================================================

# Bar plot: tf-idf HMDB and Drugbank   

#===========================================================================

human_title_year2 <- hmdb_drug %>% 
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

db.labs <- c("HMDB", "DrugBank") 
names(db.labs) <- c("human", "drugbank") 

human_title_tf_idf %>% 
    group_by(database) %>% 
    slice_max(tf_idf, n = 10) %>% 
    ungroup() %>%
    ggplot(aes(tf_idf, fct_reorder(word, tf_idf), fill = database)) +
    geom_col(show.legend = FALSE) + 
    facet_wrap(~database, 
               ncol = 2, 
               scales = "free", 
               labeller = labeller(database = db.labs)) +
    labs(x = "tf-idf", y = NULL, 
         caption = 'The words in titles specific to multi-omics citations of DrugBank and HMDB') +
    theme_minimal() +
    theme(plot.background = element_blank(),
          panel.grid = element_blank(),
          plot.caption.position = "plot") + 
    scale_fill_manual(values = c("#FF00F7", "#1E98FD"))


#===========================================================================

# Line plot: citation over time 

#===========================================================================

merged_df_unique_tag %>% 
    filter(database %in% c("human")) %>% 
    filter(Publication.Year >= 0 & Publication.Year < 2023) %>% 
    filter(is_multiomics == "yes") %>% 
    group_by(database, Publication.Year) %>% 
    count(is_multiomics) %>% 
    ungroup() %>% 
    group_by(database) %>% 
    mutate(cumulative_sum = cumsum(n)) %>%  
    mutate(versions = c(NA, NA, 'HMDB 4.0', NA, NA, NA, 'HMDB 5.0')) %>% 
    ggplot(aes(x = Publication.Year, y = cumulative_sum, group =  database))+
    geom_line(colour = "#1E98FD") +
    geom_text(aes(label = versions), color = "grey50", nudge_x = -0.6, nudge_y = 0.8) + 
    scale_x_continuous(breaks = seq(2015, 2022, 1)) +
    labs(x = NULL, y = 'Multi-omics Citations, Cumulative') + 
    theme_classic() +
    theme(panel.grid = element_blank(),
          legend.position = "none")



#===========================================================================

# Bar plot: HMDB word frequency over time  

#===========================================================================

# tf-idf analysis 
citation_title_year <- merged_df_unique_tag %>% 
    filter(database == "human") %>% 
    filter(Publication.Year >= 0 & Publication.Year < 2023) %>% 
    filter(is_multiomics == "yes") %>% 
    mutate(after_2018 = case_when(Publication.Year >= 2018 ~ "yes",
                                  Publication.Year < 2018 ~ "no")) %>% 
    select(Title, database, Publication.Year, after_2018)  

citation_title_words <- citation_title_year %>%
    unnest_tokens(word, Title) %>% 
    count(after_2018, word, sort = TRUE) 

data("stop_words") 

citation_title_words <- citation_title_words %>% anti_join(stop_words) 

citation_total_words <- citation_title_words %>% 
    group_by(after_2018) %>% 
    count(word) %>% 
    summarize(total = sum(n)) 

citation_title_words <- left_join(citation_title_words, citation_total_words)

citation_title_tf_idf <- citation_title_words %>% 
    bind_tf_idf(word, after_2018, n) %>% 
    select(-total) %>% 
    arrange(desc(tf_idf)) 

head(citation_title_tf_idf)

list_diff_2018 <- citation_title_tf_idf %>%
    filter(after_2018 == "yes") %>% 
    slice_max(tf_idf, n = 7) %>%
    pull(word)


citations_hmdb <- merged_df_unique_tag %>%
    filter(database == "human") %>%
    filter(is_multiomics == "yes") %>% 
    select(Title, Publication.Year) 
    
citations_hmdb_token <- citations_hmdb %>% unnest_tokens(word, Title)

data("stop_words")

citations_hmdb_token1 <- citations_hmdb_token %>% anti_join(stop_words)

hmdb_w_freq <- citations_hmdb_token1 %>%
    group_by(Publication.Year) %>%
    count(word, sort = TRUE) %>% 
    ungroup()

hmdb_w_freq1 <- hmdb_w_freq %>% 
    pivot_wider(names_from = Publication.Year, 
                values_from = n, 
                values_fill = 0) %>% 
    mutate('2017' = rep(0, nrow(.))) 

hmdb_w_freq2 <- hmdb_w_freq1 %>% 
    pivot_longer(cols = '2023':'2017', 
                 names_to = 'Publication.Year', 
                 values_to = 'n') %>%  
    group_by(word) %>% 
    arrange(Publication.Year) %>%
    mutate(cumulative_n = cumsum(n)) %>% 
    ungroup()

drop_list_w <- c('approach', 'research', 'study', 'reveals')

hmdb_w_freq2 %>% 
    filter(word %in% list_diff_2018 & !(word %in% drop_list_w)) %>% 
    ggplot(aes(x = Publication.Year, y = cumulative_n, group = word)) +
    geom_point(colour = "navy") + 
    geom_line(colour = "grey50") + 
    theme_bw() + 
    labs(x = NULL, y = 'Citation Number, Cumulative', 
         caption = 'Title Trends: Growing Popularity of Words in Multi-Omics Papers Citing HMDB Since 2018.') + 
    facet_wrap(~word, scales = 'fixed', nrow = 2) +
    theme(panel.background = element_blank(),
          panel.grid = element_blank(),
          strip.text = element_text(face = 'bold', size = rel(1.5)),
          strip.background = element_rect(fill = 'lightblue', colour = 'black', size = 1))  



#===========================================================================

# Bubble plot: MetaboloAnalyst and counterparts 

#===========================================================================

# multiomics_pmcid <- total_article$PMID 

# merged_df_unique_tag <- merged_df_unique %>%
#     mutate(is_multiomics = case_when(PMID %in% multiomics_pmcid ~ "yes",
#                                      !(PMID %in% multiomics_pmcid) ~ "no"))
# 
# ma_citation <- merged_df_unique_tag %>% 
#     filter(database %in% c("metaboanalyst", "omicsnet", "paintsomics", 
#                            "xcms", "workflow", "omicsanalyst", "metabolights")) %>% 
#     filter(Publication.Year >= 0) %>% 
#     group_by(database) %>% 
#     count(is_multiomics) %>% 
#     ungroup()
# 
# head(ma_citation)

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


ma_database <- merged_df_unique[ ,c('PMID', 'database', 'Publication.Year')] 

ma_df_bubbleplt <- total_article %>% 
    right_join(ma_database, by = 'PMID') %>% 
    select(-Publication.Year.x) %>% 
    rename(Publication.Year = Publication.Year.y) %>% 
    filter(database %in% c("metaboanalyst", "omicsnet", "paintsomics", "xcms", "workflow", "omicsanalyst", "metabolights")) %>% 
    mutate(is_multiomics = case_when(PMID %in% multiomics_pmcid ~ "yes",
                                     !(PMID %in% multiomics_pmcid) ~ "no")) %>% 
    group_by(database, is_multiomics) %>% 
    count(Publication.Year) %>% 
    ungroup() %>% 
    pivot_wider(names_from = is_multiomics, values_from = n, values_fill = 0) %>% 
    rename(no_multiomics = no, yes_multiomics = yes) %>% 
    mutate(total = no_multiomics + yes_multiomics) %>% 
    group_by(database) %>% 
    arrange(Publication.Year) %>% 
    mutate(cumulative_yes_multiomics = cumsum(yes_multiomics)) %>% 
    mutate(cumulative_total = cumsum(total)) %>% 
    ungroup()


# TODO 
# X axis range 2012 to 2023 
# Legend name
# Database order in legend 
# Colors - two clusters 
# Plot background color 
# Plot grind lines 
# Hightlight MetaboloAnalyst 

ma_df_bubbleplt %>% 
    ggplot(aes(x = Publication.Year, 
               y = cumulative_yes_multiomics, 
               colour = database)) +
    geom_point(alpha = .8, 
               aes(size = cumulative_total),
               fill = 'black') +
    scale_x_continuous(breaks = seq(2012, 2023, 1)) +
    scale_color_brewer(palette = 'Set1',
                       limits = c('metaboanalyst',
                                  'xcms',
                                  'workflow',
                                  'metabolights',
                                  'omicsanalyst',
                                  'paintsomics',
                                  'omicsnet'),
                       name = 'Web apps') +
    scale_size_continuous(name = 'Total citations, cumulative') +
    guides(colour = guide_legend(order = 1),
           size = guide_legend(order = 2)) +
    theme_classic() 


#===========================================================================

# Point plot animated: MetaboloAnalyst and counterparts 

#===========================================================================
    
# TODO 
# Big point size
# Tags 
# Theme 
# Axis labels 

webapp_anim <- ggplot(ma_df_bubbleplt, 
                      aes(x = cumulative_total, 
                          y = cumulative_yes_multiomics,
                          colour = database)) +
    geom_point(size = 10) +
    scale_color_brewer(palette = 'Set1',
                       limits = c('metaboanalyst',
                                  'xcms',
                                  'workflow',
                                  'metabolights',
                                  'omicsanalyst',
                                  'paintsomics',
                                  'omicsnet'),
                       name = 'Web apps') +
    transition_time(Publication.Year) +
    labs(title = 'Year: {frame_time}', 
         x = 'Total Citations (cumulative)', 
         y = 'Multiomics Citations (cumulative)')  + 
    geom_label_repel(aes(label = database), 
                     show.legend = FALSE, 
                     size = 5, 
                     nudge_x = 0.5,
                     nudge_y = 0.5) +
    theme_classic() + 
    theme(legend.position = "none") +
    theme(plot.title = element_text(color = 'black', size = 25, face = 'bold'))

animate(webapp_anim,
        duration = 15,
        fps = 40, width = 500, height = 500,
        end_pause = 80,
        renderer = gifski_renderer())

anim_save("webapp_anim2.gif")

    

#===========================================================================

# Heat map: MetaboloAnalyst and counterparts 

#===========================================================================

sub_by_db <- merged_df_unique_tag %>% 
    filter(database %in% c("metaboanalyst", "omicsnet", "paintsomics", 
                           "xcms", "workflow", "omicsanalyst", "metabolights"))

tidy_sub_by_db <- sub_by_db %>%
    unnest_tokens(word, Title)

data("stop_words")  

tidy_sub_by_db <- tidy_sub_by_db %>% 
    anti_join(stop_words) 

head(tidy_sub_by_db)

# Calculate word frequency
frequency_sub_db <- tidy_sub_by_db %>% 
    group_by(database) %>% 
    count(word, sort = TRUE) %>% 
    mutate(prop = n / sum(n)) %>%
    select(-n)

frequency_sub_db2 <- frequency_sub_db %>% 
    pivot_wider(names_from = database, values_from = prop, values_fill = 0) 

freq_matrix <- frequency_sub_db2 %>% 
    column_to_rownames("word") %>% 
    as.matrix()

res <- cor(freq_matrix)

heatmap(x = res, 
        # col = cm.colors(49, alpha = .8), 
        col = hcl.colors(49, palette = 'Light Grays', alpha = 1), 
        symm = TRUE,
        cexRow = 1.6,
        cexCol = 1.6,
        margins = c(9.5, 9.5),
        RowSideColors = c('blue', 'red', 'red', 'red', 'yellow', 'yellow', 'yellow'),
        Colv = NA) 


#===========================================================================

# Bar plot: MetaboloAnalyst and counterparts 

#===========================================================================

webapp_cluster1 <- c('omicsanalyst', 'paintsomics', 'omicsnet')
webapp_cluster2 <- c('xcms', 'workflow', 'metabolights') 

ma_clusters <- merged_df_unique_tag %>% 
    filter(database %in% webapp_cluster1 | database %in% webapp_cluster2 | database == 'metaboanalyst') %>% 
    mutate(cluster = case_when(database %in% webapp_cluster1 ~ "c1",
                               database %in% webapp_cluster2 ~ "c2",
                               database == 'metaboanalyst' ~ 'metaboanalyst')) %>% 
    select(Title, database, cluster) 

ma_clusters_token <- ma_clusters %>%
    unnest_tokens(word, Title) %>% 
    count(cluster, word, sort = TRUE) 

ma_clusters_token_total <- ma_clusters_token %>% 
    group_by(cluster) %>% 
    summarize(total = sum(n)) 

citation_title_words_ma <- left_join(ma_clusters_token_total, ma_clusters_token)

citation_title_tf_idf_ma <- citation_title_words_ma %>% 
    bind_tf_idf(word, cluster, n) %>% 
    select(-total) %>% 
    arrange(desc(tf_idf)) 

citation_title_tf_idf_ma %>%
    group_by(cluster) %>%
    slice_max(tf_idf, n = 10) %>%
    ungroup() %>%
    ggplot(aes(tf_idf, fct_reorder(word, tf_idf), fill = cluster)) +
    geom_col(show.legend = FALSE) +
    facet_wrap(~cluster, ncol = 3, scales = "free") +
    labs(x = "tf-idf", y = NULL) 

# n-gram for metaboloAnalyst 
citation_title_bigrams_ma <- ma_clusters %>% 
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
    count(cluster, bigram) %>% 
    bind_tf_idf(bigram, cluster, n) %>% 
    arrange(desc(tf_idf)) 

bigram_tf_idf_ma %>% 
    group_by(cluster) %>% 
    slice_max(tf_idf, n = 7) %>% 
    ungroup() %>%
    ggplot(aes(tf_idf, fct_reorder(bigram, tf_idf), fill = cluster)) +
    geom_col(show.legend = FALSE, width = 0.7) + 
    facet_wrap(~factor(cluster, levels = c('c2', 'metaboanalyst', 'c1'),
                       labels = c("Cluster 1", "MetaboAnalyst", "Cluster 2")), 
               ncol = 3, 
               scales = "free") +
    scale_fill_manual(values = c('yellow', 'red', 'blue')) + 
    theme_minimal() + 
    theme(strip.text = element_text(face = 'bold', size = rel(1.5)),
          strip.background = element_rect(fill = 'lightblue', colour = 'black', size = 1)) + 
    theme(panel.grid  = element_blank(),
          axis.line.x.bottom = element_line(colour = "grey50"),
          axis.ticks.x.bottom = element_line(colour = "grey50")) + 
    theme(panel.spacing = unit(0.9, 'cm', data = NULL)) +  
    theme(axis.text.y = element_text(size = 11)) + 
    theme(axis.text.y = element_text(hjust = 1)) + 
    labs(x = "tf-idf", y = NULL, 
         caption = 'Words specific to different web apps in titles of multi-omics ciations (2012-2023)') +
    theme(plot.caption = element_text(size = 10, colour = "grey50")) 






