library(tidyverse)
library(readxl)

setwd(paste0("/Users/ikaialeleiwi/Desktop/Lab/Salmonella_NIH/Lactobacillus/",
             "Omics/Metagenome/rRNA/16S_From_All_Bins/MQHQ_Picrust2/",
             "MQHQ_Picrust2/"))

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
picrust_ec_dbcan <- read_tsv("Clean_Data/picr_ec_dbcan.tsv")

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


picrust_dbcan_long <- picrust_ec_dbcan %>%
  left_join(asv_bins, by = c("sequence" = "ASV_header")) %>%
  group_by(bin.id, dbcan_column) %>%
  summarise(total_copy_number = sum(count)) %>%
  filter(bin.id %in% asv_bins_bin.id)


#normalize picrust data to # of ASVs contributing to each bin
picrust_ko_long_norm <- picrust_ko_long %>%
  left_join(asv_per_bin, by = "bin.id") %>%
  mutate(total_copy_number = ifelse(total_copy_number > 0,
                                    total_copy_number/asvs_per_bin,
                                    total_copy_number)) %>%
  select(bin.id, gene_id, total_copy_number) %>%
  ungroup() %>%
  rename("counts" = "total_copy_number")

picrust_dbcan_long_norm <- picrust_dbcan_long %>%
  left_join(asv_per_bin, by = "bin.id") %>%
  mutate(total_copy_number = ifelse(total_copy_number > 0,
                                    total_copy_number/asvs_per_bin,
                                    total_copy_number)) %>%
  select(bin.id, dbcan_column, total_copy_number) %>%
  ungroup()


#function to turn string into character vector based on separator
str_to_vect <- function(string, sep = ", "){
  string_list <- strsplit(string, split = sep)
  return(unlist(string_list))
}

#pull value from df
pullVal <- function(dataframe, column = "total_copy_number", row = 1){
  y <- dataframe %>%
    slice(row) %>%
    pull(column)
  return(y)
}


#build wider df from a tibble row
buildWide <- function(dataframe){
  df_out <- dataframe
  
  vect <- pullVal(dataframe, column = "dbcan_column")
  string <- str_to_vect(vect)
  copy_num <- pullVal(dataframe)
  
  temp_df <- data.frame(matrix(data = copy_num, 
                               nrow = 1,
                               ncol = length(string)))
  colnames(temp_df) <- string
 
  return(cbind(df_out, temp_df))
}


#combine wide rows with different column names
combineRows <- function(dataframe){
 
  for(i in 1:nrow(dataframe)){
    if(i == 1){
    working_df <- buildWide(dataframe[i,])
    }else{
    new_row <- dataframe %>%
      slice(i)
    row_to_bind <- buildWide(new_row)
    working_df <- bind_rows(working_df, row_to_bind)
    }
  }
  return(working_df)
}

#make a column with count data for each dbcan id in each row
picrust_dbcan_wide <- combineRows(picrust_dbcan_long_norm)

#transform df to long format and keep the highest count value for
#each bin:dbcan_id pair when there are duplicate dbcan_ids for a bin
picrust_dbcan_single_ids <- picrust_dbcan_wide %>%
  pivot_longer(values_to = "counts",
               names_to = 'gene_id',
               cols = 4:205) %>%
  replace_na(list(counts = 0)) %>%
  group_by(bin.id, gene_id) %>%
  summarise(max_count = max(counts)) %>%
  rename("counts" = "max_count")


#Combine picrust KO and picrust ec/dbcanid data

picrust_ko_dbcan <- rbind(picrust_ko_long_norm, picrust_dbcan_single_ids) 
  

##DRAM KO and Cazy gene counts
#summarize to total copy number for each gene and bin
#filter to include only MQHQ bins that match to an ASV

metab_summary_long <- metab_summary_dedup %>%
  select(SheetName, gene_id, gene_description, module, header, subheader,
         all_of(asv_bins_bin.id)) %>%
  pivot_longer(cols = c(contains("bin"), contains("KE")),
               values_to = "counts",
               names_to = "bin.id") %>%
  filter(str_starts(gene_id, "K") | header == "CAZY") %>%
  select(bin.id, gene_id, counts) %>%
  arrange(bin.id, gene_id)


#combine DRAM and picrust dataframes
picr_dram <- picrust_ko_dbcan %>%
  left_join(metab_summary_long, 
            by = c("bin.id", "gene_id"), 
            suffix = c("_picrust", "_dram"))

#fill NA's in total_copy_number_dram with 0's
picr_dram <- picr_dram %>%
  mutate(counts_dram = replace_na(counts_dram, 0)) 

#write out datframe
write_tsv(picr_dram, "Clean_Data/picr_dram.tsv")








