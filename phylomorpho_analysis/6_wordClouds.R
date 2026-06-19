#https://cran.r-project.org/web/packages/wordcloud2/vignettes/wordcloud.html#lettercloud-function
library(wordcloud)
library(tm)
library(dplyr)
library(tidyverse)

#rm(list=ls())
setwd("/Users/huiqing/Documents/aedes_coi/7genes")
biplot<-read.csv('mca_adult_biplot_imptChr_dim1-3_20250203.csv', header = T)

# Dimension 1 #####
dim1_pos <- biplot %>% filter(dim=="dim1") %>% filter(dim.direction=="positive")%>%select(characterFull)
dim1_neg <- biplot %>% filter(dim=="dim1") %>% filter(dim.direction=="negative")%>%select(characterFull)

# create a corpus
pos_corpus <- Corpus(VectorSource(dim1_pos))%>% tm_map(removePunctuation)%>% tm_map(stripWhitespace)%>% tm_map(removeWords, stopwords("english")) %>% tm_map(removeWords, c("caaa"))
inspect(pos_corpus)
neg_corpus <- Corpus(VectorSource(dim1_neg))

# clean the dataset. Remove special characters, numbers, punctuation
pos_clean <- pos_corpus %>% tm_map(removePunctuation)%>% tm_map(stripWhitespace)%>% tm_map(removeWords, stopwords("english")) %>% tm_map(removeWords, c("caaa"))
inspect(pos_clean)

pos_2 <- as.matrix(TermDocumentMatrix(pos_clean, control=list(tolower=F)))
matrix.pos <- sort(rowSums(pos_2),decreasing=TRUE) 
df.pos <- data.frame(word = names(matrix.pos),freq=matrix.pos)

# plot pos word cloud
set.seed(1234) # for reproducibility 
pdf(file = "wordCloud_dim1_pos.pdf", width = 8, height = 7)
wordcloud(words = df.pos$word, freq = df.pos$freq,
          max.words=200, random.order=FALSE, colors="#eddb00", scale=c(4,0.2), min.freq = 1)
dev.off()

# negative dataset
neg_clean <- neg_corpus %>% tm_map(removePunctuation) %>% tm_map(stripWhitespace)%>% tm_map(removeWords, stopwords("english")) %>% tm_map(removeWords, c("caaa"))
inspect(neg_clean)

neg_2 <- as.matrix(TermDocumentMatrix(neg_clean, control=list(tolower=F)))
matrix.neg <- sort(rowSums(neg_2),decreasing=TRUE) 
df.neg <- data.frame(word = names(matrix.neg),freq=matrix.neg)

# plot neg word cloud
set.seed(1234) # for reproducibility 
pdf(file = "wordCloud_dim1_neg.pdf", width = 8, height = 7)
wordcloud(words = df.neg$word, freq = df.neg$freq,
          max.words=200, random.order=FALSE,
          colors="#00306FFF", scale=c(4,0.2), min.freq = 1)
dev.off()

# Plot comparison word cloud #####
# Combine the 'text' columns from both data frames
combined_text <- c(dim1_pos, dim1_neg)

# Create a Corpus from the combined text and clean the dataset
combined_corpus <- Corpus(VectorSource(combined_text))%>% tm_map(removePunctuation) %>% tm_map(stripWhitespace)%>% tm_map(removeWords, stopwords("english")) %>% tm_map(removeWords, c("caaa"))
inspect(combined_corpus)

term.matrix <- as.matrix(TermDocumentMatrix(combined_corpus, control=list(tolower=F)))
term.matrix
# comparison.cloud() plots words even if they are present in the same numbers in both documents
# inspect term.matrix2 and manually remove those words
combined_corpus <- tm_map(combined_corpus, removeWords,
                          c("abdomenTergumSetae","headDecumbentScales",
                            "headInterocularScales","headMaxillaryPalpusSetae","headSetaeMales",
                            "headVertex", "setaeMales", "thoraxProepisternalSetae",
                            "thoraxScalesNarrow"))
term.matrix <- as.matrix(TermDocumentMatrix(combined_corpus, control=list(tolower=F)))

pdf(file = "wordCloud_dim1_diff.pdf", width = 8, height = 7)
comparison.cloud(term.matrix,max.words=250,random.order=FALSE,
                 title.size=0.5,scale=c(4,0.1),colors=c("#eddb00","#00306FFF"))
dev.off()


# Dimension 2 #####
dim2_pos <- biplot %>% filter(dim=="dim2") %>% filter(dim.direction=="positive")%>%select(characterFull)
dim2_neg <- biplot %>% filter(dim=="dim2") %>% filter(dim.direction=="negative")%>%select(characterFull)

# create a corpus and clean the dataset
pos_corpus <- Corpus(VectorSource(dim2_pos))%>% tm_map(removePunctuation)%>% tm_map(stripWhitespace)%>% tm_map(removeWords, stopwords("english")) %>% tm_map(removeWords, c("caaa"))
inspect(pos_corpus)

neg_corpus <- Corpus(VectorSource(dim2_neg))%>% tm_map(removePunctuation)%>% tm_map(stripWhitespace)%>% tm_map(removeWords, stopwords("english")) %>% tm_map(removeWords, c("caaa"))
inspect(neg_corpus)

# positive dataset
pos_2 <- as.matrix(TermDocumentMatrix(pos_corpus, control=list(tolower=F)))
matrix.pos <- sort(rowSums(pos_2),decreasing=TRUE) 
df.pos <- data.frame(word = names(matrix.pos),freq=matrix.pos)

# plot pos word cloud
set.seed(1234) # for reproducibility 
pdf(file = "wordCloud_dim2_pos.pdf", width = 8, height = 7)
wordcloud(words = df.pos$word, freq = df.pos$freq,
          max.words=200, random.order=FALSE, colors="#eddb00", scale=c(4,0.2), min.freq = 1)
dev.off()

# negative dataset
neg_2 <- as.matrix(TermDocumentMatrix(neg_corpus, control=list(tolower=F)))
matrix.neg <- sort(rowSums(neg_2),decreasing=TRUE) 
df.neg <- data.frame(word = names(matrix.neg),freq=matrix.neg)

# plot neg word cloud
set.seed(1234) # for reproducibility 
pdf(file = "wordCloud_dim2_neg.pdf", width = 8, height = 7)
wordcloud(words = df.neg$word, freq = df.neg$freq,
          max.words=200, random.order=FALSE,
          colors="#00306FFF", scale=c(4,0.2), min.freq = 1)
dev.off()


# Plot comparison word cloud #####
# Combine the 'text' columns from both data frames
combined_text <- c(dim2_pos, dim2_neg)

# Create a Corpus from the combined text and clean the dataset
combined_corpus <- Corpus(VectorSource(combined_text))%>% tm_map(removePunctuation) %>% tm_map(stripWhitespace)%>% tm_map(removeWords, stopwords("english")) %>% tm_map(removeWords, c("caaa"))
inspect(combined_corpus)
term.matrix <- as.matrix(TermDocumentMatrix(combined_corpus, control=list(tolower=F)))
term.matrix
# comparison.cloud() plots words even if they are present in the same numbers in both documents
# inspect term.matrix2 and manually remove those words
combined_corpus <- tm_map(combined_corpus, removeWords,
                          c("headMaxillaryPalpusScales","thoraxProepisternalScales",
                            "thoraxScalesPresent"))
term.matrix <- as.matrix(TermDocumentMatrix(combined_corpus, control=list(tolower=F)))

pdf(file = "wordCloud_dim2_diff.pdf", width = 8, height = 7)
comparison.cloud(term.matrix,max.words=250,random.order=FALSE,
                 title.size=0.5,scale=c(4,0.1),colors=c("#eddb00","#00306FFF"))
dev.off()


# Dimension 3 #####
dim3_pos <- biplot %>% filter(dim=="dim3") %>% filter(dim.direction=="positive")%>%select(characterFull)
dim3_neg <- biplot %>% filter(dim=="dim3") %>% filter(dim.direction=="negative")%>%select(characterFull)

# create a corpus and clean the dataset
pos_corpus <- Corpus(VectorSource(dim3_pos))%>% tm_map(removePunctuation)%>% tm_map(stripWhitespace)%>% tm_map(removeWords, stopwords("english")) %>% tm_map(removeWords, c("caaa"))
inspect(pos_corpus)

neg_corpus <- Corpus(VectorSource(dim3_neg))%>% tm_map(removePunctuation)%>% tm_map(stripWhitespace)%>% tm_map(removeWords, stopwords("english")) %>% tm_map(removeWords, c("caaa"))
inspect(neg_corpus)

# positive dataset
pos_2 <- as.matrix(TermDocumentMatrix(pos_corpus, control=list(tolower=F)))
matrix.pos <- sort(rowSums(pos_2),decreasing=TRUE) 
df.pos <- data.frame(word = names(matrix.pos),freq=matrix.pos)

# plot pos word cloud
set.seed(1234) # for reproducibility 
pdf(file = "wordCloud_dim3_pos.pdf", width = 8, height = 7)
wordcloud(words = df.pos$word, freq = df.pos$freq,
          max.words=200, random.order=FALSE, colors="#eddb00", scale=c(4,0.2), min.freq = 1)
dev.off()

# negative dataset
neg_2 <- as.matrix(TermDocumentMatrix(neg_corpus, control=list(tolower=F)))
matrix.neg <- sort(rowSums(neg_2),decreasing=TRUE) 
df.neg <- data.frame(word = names(matrix.neg),freq=matrix.neg)

# plot neg word cloud
set.seed(1234) # for reproducibility 
pdf(file = "wordCloud_dim3_neg.pdf", width = 8, height = 7)
wordcloud(words = df.neg$word, freq = df.neg$freq,
          max.words=200, random.order=FALSE,
          colors="#00306FFF", scale=c(2,0.2), min.freq = 1)
dev.off()

# Plot comparison word cloud #####
# Combine the 'text' columns from both data frames
combined_text <- c(dim3_pos, dim3_neg)

# Create a Corpus from the combined text and clean the dataset
combined_corpus <- Corpus(VectorSource(combined_text))%>% tm_map(removePunctuation) %>% tm_map(stripWhitespace)%>% tm_map(removeWords, stopwords("english")) %>% tm_map(removeWords, c("caaa"))
inspect(combined_corpus)
term.matrix <- as.matrix(TermDocumentMatrix(combined_corpus, control=list(tolower=F)))
term.matrix
# comparison.cloud() plots words even if they are present in the same numbers in both documents
# inspect term.matrix2 and manually remove those words
#combined_corpus <- tm_map(combined_corpus, removeWords,
#                          c("headMaxillaryPalpusScales","thoraxProepisternalScales","thoraxScalesPresent"))
term.matrix<- as.matrix(TermDocumentMatrix(combined_corpus, control=list(tolower=F)))

pdf(file = "wordCloud_dim3_diff.pdf", width = 8, height = 7)
comparison.cloud(term.matrix,max.words=250,random.order=FALSE,
                 title.size=0.5,scale=c(3.5,0.1),colors=c("#eddb00","#00306FFF"))
dev.off()
