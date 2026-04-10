setwd("/path/to/working/directory")
rm(list = ls())
#library(devtools)
#install_github("ropensci/rentrez")
#https://stackoverflow.com/questions/66495221/how-to-retrieve-data-using-the-rentrez-package-by-giving-a-list-of-query-names-i
#https://github.com/ropensci/rentrez

#Load required packages
library(rentrez)
library(data.table)
library(dplyr)
library(Biostrings)
library(tidyr)
library(phylotools)
library(tibble)

#Set paths and attributes
mainDir<-"/path/to/working/directory"
fastaPath<-"./fasta/" #output directory for fasta sequences
dat<-read.csv("txid.csv") #file for genus/subgenus and txid
dat$id<-as.character(dat$id)
txid.list<-dat$id

#Target genes. Unhash only relevant line to retrieve one gene at a time
#targetGene <- "AND (((COI[Gene] OR CO1[Gene] OR COXI[Gene] OR COX1[Gene]) AND (500[SLEN]:3000[SLEN])) OR complete genome[All Fields] OR mitochondrial genome[All Fields])"
targetGene <- "AND ((28S[Gene]) OR large subunit[All Fields] OR 28S[All Fields] NOT whole genome[All Fields])"

#For renaming file names. Unhash only relevant line
#gene<-"COI"
gene<-"28S"

# Set retmax
# Do a quick search of the number of sequences from target gene and taxa on genbank to decide if need to increase this number
retmax=99999

#NOTE: entrez search DOES include unverified sequences
#entrez_dbs(), finds list of available database. Using nucleotide in this case

count<-list()
#txid.list=cbind("300144","149531","1548035") #for testing loop code
#taxa="53541" #stegomyia
#taxa="1806171"
for (taxa in txid.list){
  taxaName=paste0(dat[dat[,2]==taxa,1])
  term<-paste("txid",taxa,"[Orgn]",targetGene, sep='',collapse = NULL)
  search<-entrez_search(db="nucleotide",term=term,retmax=retmax,use_history = TRUE)
  output<- tryCatch({entrez_fetch(db="nucleotide",web_history=search$web_history,rettype="fasta")},
                    error = function(e){NA})
  write(output, paste0(fastaPath,gene,"_",taxaName,"_",taxa,".fasta"))
  count[[taxa]]<-length(search$ids) #count number of IDs and store that in list
  print(paste("***done for",taxaName))
}  
 
# convert list to dataframe, left join to includ taxa names
count_data = do.call(rbind, count)
count.df<-as.data.frame(count_data)
df <- tibble::rownames_to_column(count.df, "id")
colnames(df)=c("id","count.ids")
count.df<-left_join(df,dat,by="id")
write.csv(count.df, "28S_count.csv")

# exclude taxa with 0 ids
avail<-count.df[!(count.df$count.ids %in% "0"),]
avail.list<-as.character(avail$id)
sum(avail$count.ids) #2082

##codes below are for standardizing headers and getting a tab delimited file for ease of sorting
csv<-list()
for (taxa in avail.list){
  taxaName=paste0(avail[avail[,1]==taxa,3])
  fastaFile <- readDNAStringSet(paste0(fastaPath,gene,"_",taxaName,"_",taxa,".fasta"))
  df<-data.frame(name=names(fastaFile), seq=paste(fastaFile)) #fasta to tab delimited
  #manipulate header names
  names<-as.data.frame(df$name)
  names(names)[1] <- "header"
  rename<-separate(names, header, c("V1","V2","V3","V4"), extra = "merge", sep=" ") #splits fasta names into three columns (acc numbers, genus, species) and the rest of the headers
  csv[[taxa]]<-data.frame(name=rename, seq=paste(fastaFile))
  }

# convert to dataframe and export as csv 
big_data = do.call(rbind, csv)
write.csv(big_data, "28S_all.csv")
nrow(big_data) #2082
