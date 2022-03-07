library(tidyverse)
library(phylotools)

setwd(paste0("/Users/ikaialeleiwi/Desktop/Lab/Salmonella_NIH/Lactobacillus/",
             "Omics/Metagenome/rRNA/16S_From_All_Bins/MQHQ_Picrust2/",
             "MQHQ_Picrust2/"))

#Get fasta of unique ASVs from Rory's ASV_bin_match pipeline for input into
#Picrust2 and filter feature table to include only those asvs

#data
match_statistics <- read_tsv("Data/match_statistics.tsv")
asvs <- read.fasta("Data/r1-r5_16S_ASV_sequences.fasta")
feature_table  <- read_tsv("Data/feature_table_r1-r5.tsv")


#unique ASV names that are matched to bins
unique_asvs <- match_statistics %>%
  pull(ASV_header) %>%
  unique()

#filter feature_table to include only unique asvs
unique_feature_table <- feature_table %>%
  filter(`#OTU ID` %in% unique_asvs)
  
which((unique_asvs %in% feature_table[,1]))
  
#write out new unique asv fasta
dat2fasta(uniqe_asvs_seqs, outfile = "Clean_Data/uniqe_asvs_with_match.fasta")
