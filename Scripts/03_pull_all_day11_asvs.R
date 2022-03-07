library(tidyverse)
library(phylotools)

setwd(paste0("/Users/ikaialeleiwi/Desktop/Lab/Salmonella_NIH/Lactobacillus/Omics/",
      "Metagenome/rRNA/16S_From_All_Bins/ASVs_to_Bins_MMSEQS2"))

#Get fasta of unique ASVs from day11 for input into
#Picrust2 and filter feature table to include only those asvs

#data
asvs <- read.fasta("Data/r1-r5_16S_ASV_sequences.fasta")

feature_table  <- read_tsv("Data/feature_table_r1-r5.tsv") %>%
  rename("ASV_id" = "#OTU ID")

hr_ctrl_samples <- read_csv("Clean_Data/hr_ctrl_samples.csv")

#filter feature table to hr_ctrl samples from day 11
feature_table_11 <- feature_table %>%
  select(ASV_id, all_of(hr_ctrl_samples$Sample))

#check for any ASVs with 0 in every sample
keep <- rowSums(feature_table_11[,2:ncol(feature_table_11)]) >= 0
sum(keep) == length(keep)

#check for duplicate ASVs
length(unique(feature_table_11$ASV_id)) == nrow(feature_table_11)

#filter rep_seqs to ASVs in feature_table_11

asvs_filtered <- asvs %>%
  filter(seq.name %in% feature_table_11$ASV_id)

#write out new unique asv fasta
dat2fasta(asvs_filtered, outfile = "Clean_Data/asvs_in_day11_hr_ctrl.fasta")

#write out filtered feature table
write_tsv(feature_table_11, "Clean_Data/feature_table_11.tsv")
