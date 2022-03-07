library(tidyverse)

setwd(paste0("/Users/ikaialeleiwi/Desktop/Lab/Salmonella_NIH/Lactobacillus/",
             "Omics/Metagenome/rRNA/16S_From_All_Bins/MQHQ_Picrust2/",
             "MQHQ_Picrust2/"))

#Data

dbcan <- read_tsv("Data/dbcan_db_out_sorted.tsv")

#Clean Data

dbcan_clean <- dbcan %>%
  mutate(target_name = str_replace(target_name, ".hmm", ""),
         target_name = str_replace(target_name, "_.*", ""),
         target_name = str_replace(target_name, "CBM35inCE17", "CBM35"),
         query_name = str_replace(query_name, "_k.*", ""),
         query_name = str_replace(query_name, "_scaffold.*", "")) %>%
  group_by(target_name, query_name) %>%
  summarise(count = n()) %>%
  ungroup() %>%
  pivot_wider(names_from = "query_name",
              values_from = "count",
              values_fill = 0)


#Write out df
write_tsv(dbcan_clean, "Clean_Data/GH_in_MQHQ_bins.tsv")
