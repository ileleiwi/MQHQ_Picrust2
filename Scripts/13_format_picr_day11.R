library(tidyverse)

setwd(paste0("/Users/ikaialeleiwi/Desktop/Lab/Salmonella_NIH/Lactobacillus/",
             "Omics/Metagenome/rRNA/16S_From_All_Bins/MQHQ_Picrust2/",
             "MQHQ_Picrust2/"))

#data
ec <- read_tsv("Data/EC_predicted_day11.tsv") %>%
  as.data.frame() %>%
  pivot_longer(cols = starts_with("EC"),
               values_to = "counts",
               names_to = "gene_id") %>%
  mutate(gene_id = str_remove(gene_id, ":")) %>%
  rename("bin.id" = "sequence")

ko <- read_tsv("Data/KO_predicted_day11.tsv") %>%
  as.data.frame() %>%
  pivot_longer(cols = starts_with("K"),
               values_to = "counts",
               names_to = "gene_id") %>%
  rename("bin.id" = "sequence")

#combine ec and ko predicetd counts
picr_day11 <- rbind(ko, ec)

#write out data
write_tsv(picr_day11, "Clean_Data/picr_day11.tsv")