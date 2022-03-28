## ---------------------------
##
## Script: 23_format_GH_picrust_day11
##
## Purpose: Create feature table of GH genes and picrust predicted genomes from day11
##
## Author: Ikaia Leleiwi
##
## Date Created: 2022-03-28
##
## Copyright (c) Ikaia Leleiwi, 2022
## Email: ileleiwi@gmail.com
##
## ---------------------------
##
## Notes: requires picr_ec_dbcan_day11.tsv
##   
##
## ---------------------------

## set working directory

setwd(paste0("/Users/ikaialeleiwi/Desktop/Lab/Salmonella_NIH/Lactobacillus/",
             "Omics/Metagenome/rRNA/16S_From_All_Bins/MQHQ_Picrust2/",
             "MQHQ_Picrust2/"))

## ---------------------------

library(tidyverse)

#data
picrust_ec_dbcan <- read_tsv("Clean_Data/picr_ec_dbcan_day11.tsv")


###functions###

str_to_vect <- function(string, sep = ", "){
  #function to turn string into character vector based on separator
  string_list <- strsplit(string, split = sep)
  return(unlist(string_list))
}


pullVal <- function(dataframe, column = "count", row = 1){
  #pull value from df
  y <- dataframe %>%
    slice(row) %>%
    pull(column)
  return(y)
}


buildWide <- function(dataframe){
  #build wider df from a tibble row
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



combineRows <- function(dataframe){
  #combine wide rows with different column names
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
picrust_dbcan_wide <- combineRows(picrust_ec_dbcan)

#transform df to long format and keep the highest count value for
#each bin:dbcan_id pair when there are duplicate dbcan_ids for a bin
picrust_dbcan_single_ids <- picrust_dbcan_wide %>%
  pivot_longer(values_to = "counts",
               names_to = 'gene_id',
               cols = 5:206) %>%
  replace_na(list(counts = 0)) %>%
  group_by(sequence, gene_id) %>%
  summarise(max_count = max(counts)) %>%
  rename("counts" = "max_count")


#make wide format dbcan single ids with bins as columns and dbcan ids as rows
GH_in_picrust_day11 <- picrust_dbcan_single_ids %>%
  pivot_wider(names_from = sequence,
              values_from = counts) %>%
  rename("target_name" = "gene_id")

write_tsv(GH_in_picrust_bins, "Clean_Data/GH_in_picrust_day11.tsv")

