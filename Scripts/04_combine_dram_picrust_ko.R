library(tidyverse)
library(readxl)

setwd(paste0("/Users/ikaialeleiwi/Desktop/Lab/Salmonella_NIH/Lactobacillus/",
             "Omics/Metagenome/rRNA/16S_From_All_Bins/MQHQ_Picrust2/"))

##Read in data
#dram gene counts
metab_summary <- "Data/metabolism_summary_deduped.xlsx" %>%
  excel_sheets() %>%
  set_names() %>%
  map_dfr(read_excel, path = "Data/metabolism_summary_deduped.xlsx", 
          .id = "SheetName") 

#deduplicated metab_summary gene_id column
metab_summary_dedup <- metab_summary %>%
  distinct(gene_id, .keep_all = TRUE)

#picrust2 predicted KO and EC counts per ASV
#pictrust2 was run on uniuqe ASV sequences from ASV_bin match
picrust_ko <- read_tsv("Data/KO_predicted.tsv")
picrust_ec <- read_tsv("Data/EC_predicted.tsv")

#ASV to Bins match statistics
asv_bins <- read_tsv("Data/match_statistics.tsv") %>%
  mutate(bin.id = sub(x = bin_scaffold_header, 
                      pattern = "_scaffold.*|_k.*",
                      replacement = "")) # <---- create column with bin name only

#Bins stats MQHQ
bin_stats <- read_tsv("Data/V5_Bin_stats.tsv")

##Filter data
#get unique bin.id's that match to an ASV and are medium or high quality
asv_bins_bin.id <- asv_bins %>%
  filter(bin.id %in% bin_stats$bin.id) %>%
  pull(bin.id) %>%
  unique()

#get number of ASVs per unique MQHQ bin for normalizaiton
asv_per_bin <- asv_bins %>%
  group_by(bin.id) %>%
  summarise(asvs_per_bin = n()) %>%
  filter(bin.id %in% asv_bins_bin.id)

##PICRUST2 predicted gene copy numbers
#summarize to total copy number for each KO and bin
#filter to include only MQHQ bins that match to an ASV

picrust_ko_long <- picrust_ko %>%
  pivot_longer(cols = starts_with("K"),
               values_to = "copy_number",
               names_to = "gene_id") %>%
  left_join(asv_bins, by = c("sequence" = "ASV_header")) %>%
  group_by(bin.id, gene_id) %>%
  summarise(total_copy_number = sum(copy_number)) %>%
  filter(bin.id %in% asv_bins_bin.id) %>%
  arrange(bin.id, gene_id)


picrust_ec_long <- picrust_ec %>%
  pivot_longer(cols = starts_with("EC"),
               values_to = "copy_number",
               names_to = "gene_id") %>%
  left_join(asv_bins, by = c("sequence" = "ASV_header")) %>%
  group_by(bin.id, gene_id) %>%
  summarise(total_copy_number = sum(copy_number)) %>%
  filter(bin.id %in% asv_bins_bin.id)


#normalize picrust data to # of ASVs contributing to each bin
picrust_ko_long_norm <- picrust_ko_long %>%
  left_join(asv_per_bin, by = "bin.id") %>%
  mutate(total_copy_number = ifelse(total_copy_number > 0,
                                    total_copy_number/asvs_per_bin,
                                    total_copy_number))

picrust_ec_long_norm <- picrust_ec_long %>%
  left_join(asv_per_bin, by = "bin.id") %>%
  mutate(total_copy_number = ifelse(total_copy_number > 0,
                                    total_copy_number/asvs_per_bin,
                                    total_copy_number))

##DRAM KO gene counts
#summarize to total copy number for each KO and bin
#filter to include only MQHQ bins that match to an ASV

metab_summary_long <- metab_summary_dedup %>%
  select(SheetName, gene_id, gene_description, module, header, subheader,
         all_of(asv_bins_bin.id)) %>%
  pivot_longer(cols = c(contains("bin"), contains("KE")),
               values_to = "copy_number",
               names_to = "bin.id") %>%
  filter(str_starts(gene_id, "K")) %>%
  select(bin.id, gene_id, copy_number) %>%
  rename("total_copy_number" = "copy_number") %>%
  arrange(bin.id, gene_id)
  

#combine DRAM and picrust dataframes
picr_dram_ko <- picrust_ko_long_norm %>%
  left_join(metab_summary_long, 
            by = c("bin.id", "gene_id"), 
            suffix = c("_picrust", "_dram"))

#fill NA's in total_copy_number_dram with 0's
picr_dram_ko <- picr_dram_ko %>%
  mutate(total_copy_number_dram = replace_na(total_copy_number_dram, 0)) 

write_tsv(picr_dram_ko, "Clean_Data/picr_dram_ko.tsv")