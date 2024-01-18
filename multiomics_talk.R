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

#===============================================

# Entire publication up to 2013 

#=============================================== 

total_article <- read.csv("./all_articles_multiomics.csv") 

# Create a vector containing only the text
table(total_article$Publication.Year)   


#---------------------Line plot entire papers---------------------------------

head(total_article)

summary_total_article <- total_article %>% 
    group_by(Create.Date) %>% 
    count(count = n()) %>% 
    ungroup() %>% 
    mutate(cum_count = cumsum(n))

# summary_total_article %>% 
#     filter(Create.Date < 2024) %>% 
#     ggplot(aes(x = Create.Date, y = cum_count)) +
#     geom_col() +
#     geom_vline(aes(xintercept = "2020/02/11"), colour = "red") 
#     # geom_smooth(method = "lm", formula = y ~ poly(x, 3), se = FALSE) 
    

# Publication in 2022 
get_word_frequnecy <- function(input, year) {
    
    total_article0 <- input
    
    total_article_2022 <- total_article0 %>% filter(Publication.Year == year) 
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



# Publication 2022 
df_all_2022 <- get_word_frequnecy(input = total_article, 2022) 
    
df_all_2022 <- df_all_2022 %>% 
    mutate(total = sum(freq)) %>%
    mutate(percentage = freq/total) %>%
    arrange(desc(percentage))


# Publication in 2023  
df_all_2023 <- get_word_frequnecy(input = total_article, 2023)

df_all_2023 <- df_all_2023 %>% 
    mutate(total = sum(freq)) %>%
    mutate(percentage = freq/total) %>%
    arrange(desc(percentage))


#------------2022 vs 2023------------------------

# Word rates 
df_all_2022_trim <- df_all_2022 %>% 
    filter(!(word %in% unwanted_word_list)) 
hot_word_2022_all <- df_all_2022_trim$word[1:20]  


df_all_2023_trim <- df_all_2023 %>% 
    filter(!(word %in% unwanted_word_list)) 
hot_word_2023_all <- df_all_2023_trim$word[1:20]  

hot_word_2023_all[!(hot_word_2023_all %in% hot_word_2022_all)]   





# Compare 2022 versus 2023 
compare_2022_2023 <- df_all_2022 %>% inner_join(df_all_2023, by = "word")

compare_2022_2023_percentage <- compare_2022_2023 %>% 
    select(1, 4, 7) %>% 
    gather(key ="year", value = "word_percent", 2:3)

# ggplot(compare_2022_2023_percentage, aes(x = word, y = word_percent, fill = year)) +
#     geom_bar(stat = "identity", position = "dodge") +
#     coord_flip() +
#     facet_wrap(~year)

compare_2022_2023_diff <- compare_2022_2023 %>% mutate(diff = (percentage.y - percentage.x)/percentage.x) 

compare_2022_2023_diff_top <- compare_2022_2023_diff %>% 
    arrange(desc(diff)) 

compare_2022_2023_diff_top$word[which(compare_2022_2023_diff_top$word %in% hot_word_2023)] 
compare_2022_2023_diff_top$diff[which(compare_2022_2023_diff_top$word %in% hot_word_2023)] 


# ggplot(compare_2022_2023_diff_top, aes(x = word, y = diff))+
#     geom_bar(stat = "identity") +
#     coord_flip()






#----------------------------Word cloud---------------------------------


# Create an unwanted word list - including verb, adjective, preposition. 
unwanted_word_list <- c("and", "the", "for", "data", "analysis", "approach", "from", "reveals", 
                        "using", "approaches", "study", "with", "based", "analyses", 
                        "altered", "identifies", "novel", "identify", "profiling", "mining",
                        "modeling", "revealing", "reveal", "against", "via", "learning", "targeting",
                        "through", "into", "exploration", "repurposing", "expanding", "exploring",
                        "defining", "integrative", "clinical", "potential", "revealed", "associated", "new") 

index_all <- which(df_all$word %in% unwanted_word_list)
df2_all <- df_all[-index_all, ]

set.seed(1234)
wordcloud(words = df2_all$word,
          freq = df2_all$freq,
          scale = c(3, 1.5),
          min.freq = 5,
          max.words = 50,
          random.order = F,
          random.color = F,
          col=terrain.colors(length(df2_all$word) , alpha=0.9) , rot.per=0.3 )

my_graph <- wordcloud2(df2_all, size= 1.6) 
saveWidget(my_graph,"tmp.html",selfcontained = F)
webshot("tmp.html","fig_1.pdf", delay =5, vwidth = 480, vheight=480)  


#==============================================================

# Task 2: How many TMIC paper in total multiomics papers?

#==============================================================
tmic_ciation <- read.csv("./citation_pubmed_files/tmic-paper.csv")  
# sum(tmic_ciation$PMID %in% total_article$PMID == 1) 




#==============================================================

# Task 3: Database ranking 

#==============================================================
file_list <- list.files("./citation_pubmed_files") 

# Read in individual csv files 
read_file <- function (file_path, file_list_element) {
    single_df <- read.csv(file_path)
    repeat_n <- nrow(single_df)
    db_name <- str_extract(file_list_element, "^[a-z]+") 
    single_df$database <- rep(db_name, repeat_n) 
    return(single_df) 
}

master_list <- list()
for (j in 1:length(file_list)) {
    file_path <- paste0("./citation_pubmed_files/", as.character(file_list[j]))
    master_list[[j]] <- read_file(file_path, file_list[j]) 
}

# Merge all data frames 
i <- 2
merged_df <- master_list[[1]]
m <- length(master_list) + 1 

while (i < m) {
    
    print(i)

    merged_df <- rbind(merged_df, master_list[[i]])

    i <- i + 1
}

# Remove duplicate rows 
merged_df_unique <- unique(merged_df)


table(merged_df_unique$database, merged_df_unique$Publication.Year)


# Citation of multiomics 
citation_multiomics <- merged_df_unique[(merged_df_unique$PMID %in% total_article$PMID), ] 

table(citation_multiomics$database)

citation_multiomics <- citation_multiomics %>%
    dplyr::filter(database %in% c("drugbank", "human", "metaboanalyst"))


head(citation_multiomics)




#--------------------Human metabolism database---------------------------------

citation_multiomics_human <- citation_multiomics %>% filter(database == "human")

human_word_2022 <- get_word_frequnecy(citation_multiomics_human, 2022) 
human_word_2023 <- get_word_frequnecy(citation_multiomics_human, 2023) 


# Word rates 
human_word_2022_top <- human_word_2022 %>% 
    mutate(total = sum(freq)) %>%
    mutate(percentage = freq/total) %>%
    arrange(desc(percentage)) %>% 
    filter(!(word %in% unwanted_word_list)) 

hot_word_2022 <- human_word_2022_top$word[1:20] 

human_word_2023_top <- human_word_2023 %>% 
    mutate(total = sum(freq)) %>%
    mutate(percentage = freq/total) %>%
    arrange(desc(percentage)) %>% 
    filter(!(word %in% unwanted_word_list)) 

hot_word_2023 <- human_word_2023_top$word[1:20]  

hot_word_2023[hot_word_2023 %in% hot_word_2023_all] 
hot_word_2022[hot_word_2022 %in% hot_word_2022_all] 

citation_multiomics_human$Title[citation_multiomics_human$Publication.Year == 2023]


#=========================================================

# Motivation - Identify interesting points 

#=========================================================

# Create a subset for human database alone 
full_citation_df <- table(merged_df_unique$database, merged_df_unique$Publication.Year) %>% as.data.frame()
names(full_citation_df) <- c("database", "year", "citations_full") 
full_citation_human <- full_citation_df %>% filter(database == "human") 

# Create a subset for human database alone 
omics_citation_df <- table(citation_multiomics$database, citation_multiomics$Publication.Year) %>% as.data.frame()
names(omics_citation_df) <- c("database", "year", "citations_omics") 
omics_citation_human <- omics_citation_df %>% filter(database == "human")  

# Merge full citation and omic citations by human and year 
merged_human_citation <- omics_citation_human %>% 
    inner_join(full_citation_human, by  = c("year", "database")) %>% 
    gather(key = "citation_type", value = "paper_number", 3:4)

# Rate by year 
omics_citation_human %>% 
    inner_join(full_citation_human, by  = c("year", "database")) %>% 
    mutate(rate = citations_omics/citations_full) %>% 
    ggplot(aes(x = year, y = rate)) +
    geom_point()

# Citation by types 
ggplot(merged_human_citation, 
       aes(x = year, y = paper_number, group = citation_type)) +
    geom_line() +
    geom_point() +
    facet_wrap(~citation_type, scales = "free_y")


df_barplot <- table(citation_multiomics$database) %>% as.data.frame()

names(df_barplot) <- c("database", "citations")

# Bar plot 
ggplot(df_barplot, aes(x = database, y = citations)) +
    geom_bar(stat = "identity")


# Paper number by year 
line_plot_year <- table(citation_multiomics$Publication.Year) %>% as.data.frame()

names(line_plot_year) <- c("year", "citations") 

ggplot(line_plot_year, aes(x = year, y = citations)) +
    geom_bar(stat = "identity")


# Citation by year and by database 
citation_year_db <- table(citation_multiomics$Publication.Year, citation_multiomics$database) %>% as.data.frame()

names(citation_year_db) <- c('year', "database", "citations")

ggplot(citation_year_db, aes(x = year, y = citations)) +
    geom_bar(stat = "identity") +
    coord_flip()



# # Create a vector containing only the text
# title_human_db <- citation_multiomics$Title
# 
# # Create a corpus
# title_human_db_corpus <- tm::Corpus(VectorSource(title_human_db))
# 
# # Clean text
# title_human_db_corpus <- title_human_db_corpus %>%
#     tm_map(stripWhitespace)
# 
# title_human_db_corpus <- tm_map(title_human_db_corpus, content_transformer(tolower))
# 
# 
# # Create a document-term-matrix
# doc_term_matrix <- TermDocumentMatrix(title_human_db_corpus)
# matrix <- as.matrix(doc_term_matrix)
# words <- sort(rowSums(matrix), decreasing = TRUE)
# 
# df <- data.frame(word = names(words), freq = words)
# 
# index <- which(df$word %in% unwanted_word_list)
# df2 <- df[-index, ]
# 
# head(df2, 15)
# 
# set.seed(1234)
# wordcloud(words = df2$word,
#           freq = df2$freq,
#           scale = c(3, 1),
#           min.freq = 5,
#           max.words = 70,
#           random.order = F,
#           rot.per = 0.45,
#           colors = brewer.pal(3, "Accent"),
#           colorblindFriendly = T)



## Text mining - Word frequency 
major_tools_index <- which(table(citation_multiomics$database) > 10)

# tmic_tools <- c("drugbank", "human", "metaboanalyst") 

major_tools <- unique(citation_multiomics$database)[major_tools_index]

word_freq_list <- list()

for (k in 1: length(major_tools)) {
    
    tool <- major_tools[k]
    
    sub_citation <- citation_multiomics[(citation_multiomics$database == tool), ] 
    
    # Create a vector containing only the text 
    title_human_db <- sub_citation$Title
    
    # Create a corpus
    title_human_db_corpus <- tm::Corpus(VectorSource(title_human_db))
    
    # Clean text 
    title_human_db_corpus <- title_human_db_corpus %>% 
        tm_map(stripWhitespace) 
    title_human_db_corpus <- tm_map(title_human_db_corpus, 
                                    content_transformer(tolower)) 
    
    # Create a document-term-matrix 
    doc_term_matrix <- TermDocumentMatrix(title_human_db_corpus) 
    matrix <- as.matrix(doc_term_matrix)
    words <- sort(rowSums(matrix), decreasing = TRUE) 
    
    df <- data.frame(word = names(words), freq = words) 
    
    index <- which(df$word %in% unwanted_word_list) 
    
    df2 <- df[-index, ]
    
    rep_number <- nrow(df2)
    
    df2$database <- rep(tool, rep_number)
    
    word_freq_list[[k]] <- df2
}


# Merge individual word freq tables 
f <- 2
merged_wordfreq_df <- word_freq_list[[1]]
limit <- length(word_freq_list) + 1 

while (f < limit) {
    
    print(f)
    
    merged_wordfreq_df <- rbind(merged_wordfreq_df, word_freq_list[[f]])
    
    f <- f + 1
}


# Reshape from long to wide 
merged_wordfreq_wider <- merged_wordfreq_df %>%
    pivot_wider(
        names_from = database, 
        values_from = freq,
        values_fill = 0)

# head(merged_wordfreq_wider)



#-------------------------------------PCA---------------------------------------

# Prepare data table 
merged_wordfreq_wider_t <- t(merged_wordfreq_wider) %>% as.data.frame()
colnames(merged_wordfreq_wider_t) <- merged_wordfreq_wider_t[1, ]
merged_wordfreq_wider_t <- merged_wordfreq_wider_t[-1, ]
merged_wordfreq_wider_t2 <- apply(merged_wordfreq_wider_t, 2, as.numeric) 
rownames(merged_wordfreq_wider_t2) <- rownames(merged_wordfreq_wider_t) 
merged_wordfreq_wider_t2 <- as.data.frame(merged_wordfreq_wider_t2)

# head(merged_wordfreq_wider_t2)

# Remove common words 
drop_words <- c("multiomics", "multi-omics", "network")

merged_wordfreq_wider_t3 <- merged_wordfreq_wider_t2[ , -which(colnames(merged_wordfreq_wider_t2) %in% drop_words)]

colnames(merged_wordfreq_wider_t3)

# Remove constant variables 
constant_col <- names(wrangle::constant(merged_wordfreq_wider_t3))
merged_wordfreq_wider_t4 <- merged_wordfreq_wider_t3 %>% select(-constant_col) 


# Normalization 
merged_wordfreq_normalized <- scale(merged_wordfreq_wider_t4)

# Correlation between words 
# corr_matrix <- cor(merged_wordfreq_normalized)
# pca_res2 <- princomp(corr_matrix)


# Compuate PCA 
pca_res <- prcomp(merged_wordfreq_wider_t4, scale. = F) 

## Variable 
factoextra::fviz_pca_var(pca_res)

## Contribution of variables to PCs
fviz_cos2(pca_res, choice = "var", axes = 1:2, top = 11, fill = "lightgray")

# Individuals 
fviz_pca_ind(pca_res)

# Biplot 
fviz_pca_biplot(pca_res, 
                label = "var",
                ggtheme = theme_minimal(), 
                col.var = "cos2")





#-----------------------------------Radar chart--------------------------------- 

top_word <- c("drug", "metabolites", "metabolomics", "metabolism", "target", 
              "networks", "research", "gene", "gut", "systems", "disease")

radar_df <- merged_wordfreq_normalized %>% 
    as.data.frame() %>% 
    select(top_word) 

min = rep(-1.5, 11)
max = rep(1.5, 11)

radar_df2 <- rbind(max, min, radar_df) %>% as.data.frame() 

# Color vector
colors_border=c( rgb(0.2,0.5,0.5,0.9), rgb(0.8,0.2,0.5,0.9) , rgb(0.7,0.5,0.1,0.9) )
colors_in=c( rgb(0.2,0.5,0.5,0.4), rgb(0.8,0.2,0.5,0.4) , rgb(0.7,0.5,0.1,0.4) )

fmsb::radarchart(radar_df2, axistype = 1,
                 #custom polygon
                 pcol=colors_border , pfcol=colors_in , plwd=4 , plty=1,
                 #custom the grid
                 cglcol="grey", cglty=1, axislabcol="grey", caxislabels=seq(0,20,5), cglwd=0.8,
                 #custom labels
                 vlcex=0.8)

# Add a legend
legend(x=0.7, y=1, 
       legend = rownames(radar_df2[-c(1,2),]), 
       bty = "n", pch=20 , 
       col=colors_in , 
       text.col = "grey", 
       cex=1.2, pt.cex=3)






#------------------------Publication by year------------------------------------

all_by_year <- table(total_article$Publication.Year) %>% as.data.frame()
colnames(all_by_year) <- c("year", "all")


tmic_by_year <- table(citation_multiomics$Publication.Year) %>% as.data.frame()
colnames(tmic_by_year) <- c("year", "tmic_citation")

combine_year <- tmic_by_year %>% 
    left_join(all_by_year, by = "year") 

combine_year <- combine_year %>% 
    mutate(rate = tmic_citation / all)

ggplot(combine_year, aes(x = year, y = rate)) +
    geom_bar(stat = "identity") +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    xlab("") + 
    ylab("") + 
    theme_minimal() 


#--------------------------Circle packing---------------------------------------

# Refer to https://r-graph-gallery.com/circle-packing.html 

df_circle <- table(citation_multiomics$database) %>% as.data.frame()

names(df_circle) <- c("group", "paper")

df_circle$group <- as.character(df_circle$group)

df_circle[4, ] <- c("non-tmic", 3784) 

packing <- circleProgressiveLayout(df_circle$paper, sizetype = 'area') 
df_circle <- cbind(df_circle, packing)

dat.gg <- circleLayoutVertices(packing, npoints = 50)

ggplot() + 
    geom_polygon(data = dat.gg, 
                 aes(x, y, group = id, fill=as.factor(id)), 
                 colour = "black", alpha = 0.6) +
    # scale_fill_manual(values = magma(nrow(df_circle))) +
    # geom_text(data = df_circle, aes(x, y, size=value, label = group)) +
    scale_size_continuous(range = c(1,4)) +
    theme_void() + 
    theme(legend.position="none") +
    coord_equal()


#------------------------Popular words frequency-------------------------------

# Calculate word rate = word / total words 

# df2_all # word frequency in all papers

# merged_wordfreq_df # word freq in tmic-citatoins 

# All paper 
all_paper_number <- nrow(total_article)

df2_all_rate <- df2_all %>% 
    mutate(sum_paper = all_paper_number) %>% 
    mutate(rate_all = freq / sum_paper) %>% 
    arrange(desc(rate_all))

head(df2_all_rate, 10)


# TMIC citation 
tmic_citation_number <- nrow(citation_multiomics)

grouped_tmic <- merged_wordfreq_df %>% 
    group_by(word) %>% 
    summarise(freq_3 = sum(freq)) 

head(grouped_tmic)
    
df_tmic_rate <- grouped_tmic %>% 
    mutate(sum_paper = tmic_citation_number) %>% 
    mutate(rate_tmic = freq_3/ sum_paper) %>%
    arrange(desc(rate_tmic))

head(df_tmic_rate, 10)


# Merge two data tables 

top20_df2_all_rate <-  df2_all_rate[1:20, ]

merged_word_freq <- top20_df2_all_rate %>% 
    inner_join(df_tmic_rate, by = "word") %>% 
    rename(sum_all = sum_paper.x, sum_tmic = sum_paper.y, freq_all = freq, freq_tmic = freq_3)

head(merged_word_freq)

gather_word_freq <- merged_word_freq %>% 
    select(word, rate_all, rate_tmic) %>% 
    gather(type, rate, -word) 


ggplot(gather_word_freq, aes(x = word, y = rate, fill = type)) +
    geom_bar(stat = "identity", position = "dodge") +
    coord_flip()



#---------------------------------EXPERIMENT------------------------------------

# Paper by top words of TMIC's citation 
major_tools <- unique(citation_multiomics$database)[major_tools_index]

word_matrix_list <- list()

for (k in 1: length(major_tools)) {
    
    tool <- major_tools[k]
    
    sub_citation <- citation_multiomics[(citation_multiomics$database == tool), ] 
    
    # Create a vector containing only the text 
    title_human_db <- sub_citation$Title
    
    # Create a corpus
    title_human_db_corpus <- tm::Corpus(VectorSource(title_human_db))
    
    # Clean text 
    title_human_db_corpus <- title_human_db_corpus %>% 
        tm_map(stripWhitespace) 
    title_human_db_corpus <- tm_map(title_human_db_corpus, 
                                    content_transformer(tolower)) 
    
    # Create a document-term-matrix 
    doc_term_matrix <- TermDocumentMatrix(title_human_db_corpus) 
    matrix <- as.matrix(doc_term_matrix)
    
    matrix_t <- t(matrix) %>% as.data.frame() 
    
    matrix_t$PMID <- sub_citation$PMID
    
    df_matrix <- as.data.frame(matrix_t) 
    
    df_matrix <- df_matrix %>%
        gather(word, count, -PMID) 
    
    rep_labels <- nrow(df_matrix) 
    
    df_matrix$database <- rep(tool, rep_labels)
    
    word_matrix_list[[k]] <- df_matrix
    
}

# Merge by row 
v <- 2
merged_word_matrix <- word_matrix_list[[1]]
limit2 <- length(word_matrix_list) + 1 

while (v < limit2) {
    
    print(v)
    
    merged_word_matrix <- rbind(merged_word_matrix, word_matrix_list[[v]])
    
    v <- v + 1
}


top_word <- c(pc1_loadings, pc2_loadings)

merged_word_matrix2 <- merged_word_matrix %>% 
    filter(word %in% top_word)

head(merged_word_matrix2)


## Prepare df for PCA analysis 
merged_word_matrix3 <- merged_word_matrix2 %>%
    tidyr::spread(word, count) 

# head(merged_word_matrix3)

## Replace NA with 0 
merged_word_matrix4 <- merged_word_matrix3 %>% replace(is.na(.), 0)

head(merged_word_matrix4)

str(merged_word_matrix4)

## compute PCA matrix 
tmic_word_pca <- prcomp(merged_word_matrix4[,c(3:16)],
                        center = T,
                        scale. = T) 

# Plot PCA 
autoplot(tmic_word_pca,
         data = merged_word_matrix4,
         colour = "database")

biplot(tmic_word_pca)


#-----------------------------END OF EXPERIMENT---------------------------------








