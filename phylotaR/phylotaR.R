setwd("/path/to/working/directory")
rm(list = ls())

library(phylotaR)
library(ggplot2)

# helpful tutorials
#https://docs.ropensci.org/phylotaR/articles/phylotaR.html#running
#https://docs.ropensci.org/phylotaR/reference/drop_by_rank.html

# Set paths
mainDir<-"/path/to/working/directory"
ncbi<-'/opt/homebrew/bin'#directory to ncbi

# The working directory should have the following subdirectories:
cluster.plots<-"./_plots/" #output directory for cluster plots
fasta<-"./_fasta/" #output directory for fasta sequences

# Prepare a list of taxids for taxa of interest
# Taxids can be obtained from NCBI taxonomy
# In this csv, the first column (taxa) containes taxon names, second column (id) contains taxids, and third column (rank) is the taxonomic rank of the taxa (e.g. genus, subgenus)
dat<-read.csv("txid_subset.csv")
txid.list<-dat$id

# Setup for each taxa
# Identify and download taxonomic information
# Get sequences down the taxonomic hierarchy
for (taxa in txid.list){
    wd=paste0(dat[dat[,2]==taxa,1]) # wd = corresponding taxon name
    dir.create(file.path(mainDir, wd)) # create folder for each taxon
    txid=taxa
    setup(wd = wd, txid = txid, ncbi_dr = ncbi, v = TRUE, ncps=2, overwrite = TRUE)
    print(paste0("***** Setup done for:",wd))
}

# Generate clusters down taxonomic hierarchy
for (taxa in txid.list){
    wd=paste0(dat[dat[,2]==taxa,1])
    run(wd = wd)
    print(paste0("***** Run done for:",wd))
}

# Summarize results and plot to visualise taxa vs gene (cluster) completeness
df_total = data.frame()
for (taxa in txid.list){
  tryCatch({
    wd=paste0(dat[dat[,2]==taxa,1])
    all_clusters <- read_phylota(wd)
    smmry <- data.frame(summary(all_clusters)) #summarize phylota results
    smmry$taxa <- taxa  #add taxa id
    smmry$wd <- wd #add taxa name
    df <- data.frame(smmry)
    df_total <- rbind(df_total,df) #append it to dataframe
    print(paste0("***** Summary done for:",wd))
  }, error=function(e){cat("ERROR :",conditionMessage(e), wd, "\n")}) #added tryCatch function to skip taxa with only one cluster or no sequences

     species_txids <- get_txids(all_clusters, txids = all_clusters@txids, rnk = 'species') # get species-level taxonomic names
     species_txids <- unique(species_txids) #obtain unique txids 
     species_txids <- species_txids[species_txids !=  ''] #dropping missing
     species_nms <- get_tx_slot(all_clusters, species_txids, slt_nm = 'scnm') #get species names
     species_nms <- sort(species_nms, decreasing = TRUE) #sort alphabetically for plotting
     tryCatch({
     p <- plot_phylota_pa(phylota = all_clusters, cids = all_clusters@cids, txids = species_txids,
                         txnms = species_nms) # generate geom_object
    }, error=function(e){cat("ERROR :",conditionMessage(e),wd, "\n")}) #added tryCatch function to skip taxa without sequences
     ggsave(paste0(cluster.plots,wd,".jpg"), p, width=12, height=8, units="in")
     print(paste0("***** Plot done for:",wd))
}

# Save dataframe
# Information here contains summary of all clusters
# E.g. number of clusters in each taxa, number of sequences etc.
write.csv(df_total,"smmry_all.csv")

# Obtain best fasta sequences per species for cluster 1
for (taxa in txid.list){
  wd=paste0(dat[dat[,2]==taxa,1])
  all_clusters <- read_phylota(wd)
  tryCatch({
  smmry <- data.frame(summary(all_clusters)) #summarize phylota results
  cid_keep <- smmry[1, 'ID'] #select cluster 1
  selected <- drop_clstrs(phylota = all_clusters, cid = cid_keep) #drop clusters
  reduced <- drop_by_rank(phylota = selected, rnk = 'species', n = 1, #choose best seq per species
                          choose_by = c('pambgs', 'nncltds'), #criteria: ambiguous bases, length of sequence. Excluded: age of sequence
                          greatest=c(FALSE,TRUE)) #select seq with fewest number of ambgs bases, and longest seq length
  txids <- get_txids(phylota = reduced, cid = cid_keep, rnk = 'species') # get txids at the species level for each sequence
  scientific_names <- get_tx_slot(phylota = reduced, txid = txids, slt_nm = 'scnm') # look up species names for txids
  scientific_names <- gsub('\\.', '', scientific_names) #clean names
  scientific_names <- gsub('\\s+', '_', scientific_names)
  sids <- reduced@clstrs[[cid_keep]]@sids   # Add sequence IDs
  scientific_names_cluster<-paste(scientific_names,sep="_",sids)
  write_sqs(phylota = reduced, sid = sids, sq_nm = scientific_names_cluster, 
            outfile = paste0(fasta,wd,'_c1_bestSeqs.fasta')) #export fasta
  print(paste0("***** Exported best sequences for:",wd))
  }, error=function(e){cat("ERROR :",conditionMessage(e), wd, "\n")}) #added tryCatch function to skip taxa without sequences
}

# Note: errors when extracting clusters 1-3 (cluster id 0-2) can be safely ignored.
# e.g., ERROR : [NA] not in records albuginosus - refers to groups that only have one cluster
# Modifying the code (phylotaR_fork) ensures that taxa and sequences forming one cluster are still retained.
# e.g., ERROR : 'data' must be of a vector type, was 'NULL' cancraedes - refers to groups where the txid was specified, but no sequences were available in genbank

# Obtain best fasta sequences per species for cluster 2
for (taxa in txid.list){
  wd=paste0(dat[dat[,2]==taxa,1])
  all_clusters <- read_phylota(wd)
  tryCatch({
    smmry <- data.frame(summary(all_clusters)) #summarize phylota results
    cid_keep <- smmry[2, 'ID'] #select cluster 2
    selected <- drop_clstrs(phylota = all_clusters, cid = cid_keep) #drop clusters
    reduced <- drop_by_rank(phylota = selected, rnk = 'species', n = 1, #choose best seq per species
                            choose_by = c('pambgs', 'nncltds'), #criteria: ambiguous bases, length of sequence. Excluded: age of sequence
                            greatest=c(FALSE,TRUE)) #select seq with fewest number of ambgs bases, and longest seq length
    txids <- get_txids(phylota = reduced, cid = cid_keep, rnk = 'species') # get txids at the species level for each sequence
    scientific_names <- get_tx_slot(phylota = reduced, txid = txids, slt_nm = 'scnm') # look up species names for txids
    scientific_names <- gsub('\\.', '', scientific_names) #clean names
    scientific_names <- gsub('\\s+', '_', scientific_names)
    sids <- reduced@clstrs[[cid_keep]]@sids   # Add sequence IDs
    scientific_names_cluster<-paste(scientific_names,sep="_",sids)
    write_sqs(phylota = reduced, sid = sids, sq_nm = scientific_names_cluster, 
              outfile = paste0(fasta,wd,'_c2_bestSeqs.fasta')) #export fasta
    print(paste0("***** Exported best sequences for:",wd))
  }, error=function(e){cat("ERROR :",conditionMessage(e), wd, "\n")}) #added tryCatch function to skip taxa without sequences
}

# Obtain best fasta sequences per species for cluster 3
for (taxa in txid.list){
  wd=paste0(dat[dat[,2]==taxa,1])
  all_clusters <- read_phylota(wd)
  tryCatch({
    smmry <- data.frame(summary(all_clusters)) #summarize phylota results
    cid_keep <- smmry[3, 'ID'] #select cluster 3
    selected <- drop_clstrs(phylota = all_clusters, cid = cid_keep) #drop clusters
    reduced <- drop_by_rank(phylota = selected, rnk = 'species', n = 1, #choose best seq per species
                            choose_by = c('pambgs', 'nncltds'), #criteria: ambiguous bases, length of sequence. Excluded: age of sequence
                            greatest=c(FALSE,TRUE)) #select seq with fewest number of ambgs bases, and longest seq length
    txids <- get_txids(phylota = reduced, cid = cid_keep, rnk = 'species') # get txids at the species level for each sequence
    scientific_names <- get_tx_slot(phylota = reduced, txid = txids, slt_nm = 'scnm') # look up species names for txids
    scientific_names <- gsub('\\.', '', scientific_names) #clean names
    scientific_names <- gsub('\\s+', '_', scientific_names)
    sids <- reduced@clstrs[[cid_keep]]@sids   # Add sequence IDs
    scientific_names_cluster<-paste(scientific_names,sep="_",sids)
    write_sqs(phylota = reduced, sid = sids, sq_nm = scientific_names_cluster, 
              outfile = paste0(fasta,wd,'_c3_bestSeqs.fasta')) #export fasta
    print(paste0("***** Exported best sequences for:",wd))
  }, error=function(e){cat("ERROR :",conditionMessage(e), wd, "\n")}) #added tryCatch function to skip taxa without sequences
}

# Concatenate best sequences from each cluster together and carried out blast.

# To retrieve COI sequences from cluster 3 onwards, I first exported a list of accessions from smmry_all.csv
# I used batch entrez (https://www.ncbi.nlm.nih.gov/sites/batchentrez) to query the list of accessions and exported the results
# Identified accession numbers which matched to the gene COI
# Created a csv file (taxids_add.csv) for phylotaR
# This step can also be done with entrez in command line.
dat<-read.csv("txids_add.csv")
list<-dat$sn

for (taxa in list){
  wd=paste0(unique(dat[dat[,1]==taxa,4])) ###
  c.id=as.numeric(paste0(dat[dat[,1]==taxa,2]))+1 ### +1 here because cluster x = c.id+1 row in the smmry table
  all_clusters <- read_phylota(wd)
  tryCatch({
    smmry <- data.frame(summary(all_clusters)) #summarize phylota results
    cid_keep <- smmry[(c.id), 'ID'] #select cluster id
    selected <- drop_clstrs(phylota = all_clusters, cid = cid_keep) #drop clusters
    reduced <- drop_by_rank(phylota = selected, rnk = 'species', n = 1, #choose best seq per species
                            choose_by = c('pambgs', 'nncltds'), #criteria: ambiguous bases, length of sequence. Excluded: age of sequence
                            greatest=c(FALSE,TRUE)) #select seq with fewest number of ambgs bases, and longest seq length
    txids <- get_txids(phylota = reduced, cid = cid_keep, rnk = 'species') # get txids at the species level for each sequence
    scientific_names <- get_tx_slot(phylota = reduced, txid = txids, slt_nm = 'scnm') # look up species names for txids
    scientific_names <- gsub('\\.', '', scientific_names) #clean names
    scientific_names <- gsub('\\s+', '_', scientific_names)
    sids <- reduced@clstrs[[cid_keep]]@sids   # Add sequence IDs
    scientific_names_cluster<-paste0(scientific_names,"_",sids,"_c",c.id)
    write_sqs(phylota = reduced, sid = sids, sq_nm = scientific_names_cluster, 
              outfile = paste0(fasta,wd,"_c",c.id,'_bestSeqs.fasta')) #export fasta
    print(paste0("***** Exported best sequences for:",wd))
}, error=function(e){cat("ERROR :",conditionMessage(e), wd, c.id,"\n")}) #added tryCatch function to skip taxa without sequences
}

# Error comes from Stegomyia because there were some NAs introduced (rows 7 and 8 in the txid_add.csv)
# Identified the seeds as KT358463.1, AJ971003.1
# Checked that the fasta output looks okay
