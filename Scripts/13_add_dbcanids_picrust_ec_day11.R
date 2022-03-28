## ---------------------------
##
## Script: add_dbcanids_picrust_ec_day11
##
## Purpose: Get DBCan IDs from picrust ec data for day 11 ASVs
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
## Notes: 
##   
##
## ---------------------------

## set working directory

setwd(paste0("/Users/ikaialeleiwi/Desktop/Lab/Salmonella_NIH/Lactobacillus/",
             "Omics/Metagenome/rRNA/16S_From_All_Bins/MQHQ_Picrust2/",
             "MQHQ_Picrust2/"))

## ---------------------------
library(tidyverse)
library(lubridate)


##Read in data
#DRAM cazy/EC matches from genome summary form
cazyid_EC <- read_tsv("Clean_Data/cazyid_EC.tsv") 

#predicted ECs from picrust2
picr_ec <- read_tsv("Data/EC_predicted_day11.tsv")

#extracts EC #'s from string and returns a list of matches
PullStr <- function(string, type=1){
  if(type == 1){
    #match EC number format EC -.-.-.-
    x <- str_extract_all(string, 
                         "EC [\\d-]{1,2}\\.[\\d-]{1,2}\\.[\\d-]{1,2}\\.[\\d-]{1,2}")
  }else if(type == 2){
    #match EC number format EC:-.-.-.-
    x <- str_extract_all(string, 
                         "EC:[\\d-]{1,2}\\.[\\d-]{1,2}\\.[\\d-]{1,2}\\.[\\d-]{1,2}")
  }
  return(x)
}


#collaps all elements of character vector into one string
collapseVect <- function(char_vect){
  out_string <- ""
  for (i in char_vect){
    out_string <- paste(out_string, i, sep = ", ")
  }
  out_string <- str_remove(out_string, "^, ?")
  out_string <- str_replace_all(out_string, "EC ", "EC:")
  return(out_string)
}


#pull EC string from dataframe
pullECcazyid_EC <- function(dataframe = cazyid_EC, in_row){
  ec_string <- dataframe %>%
    slice(in_row) %>%
    pull(EC)
  return(ec_string)
}

#pull dbcan string from dataframe
pullDBCANcazyid_EC <- function(dataframe = cazyid_EC, in_row){
  ec_string <- dataframe %>%
    slice(in_row) %>%
    pull(cazyid)
  return(ec_string)
}

#pull picrust EC number from dataframe
pullECpicr_EC <- function(dataframe = picr_ec_long, in_row){
  ec_string <- dataframe %>%
    slice(in_row) %>%
    pull(EC)
  return(ec_string)
}


#collects all relevant dbCAN id's that match a particular EC number
buildDBCANstring <- function(row, df.picr=picr_ec_long, df.cazyid=cazyid_EC){
  
  dbcan_vect <- c()
  query_string <- pullECpicr_EC(dataframe = df.picr, in_row = row)
  query_string_b <- paste0("\\b", query_string, "\\b")
  for(j in 1:nrow(df.cazyid)){
    if(pullDBCANcazyid_EC(dataframe = df.cazyid, in_row = j) %in% dbcan_vect){
      next
    }else{
      subject_string <- pullECcazyid_EC(dataframe = df.cazyid, in_row = j)
      if(str_detect(subject_string, query_string_b)){
        dbcan_vect <- c(dbcan_vect, pullDBCANcazyid_EC(dataframe = df.cazyid, 
                                                       in_row = j))
      }
    }
  }
  
  if(length(dbcan_vect) > 0){
    dbcan_string <- collapseVect(dbcan_vect)
  }else{
    dbcan_string <- "-"
  }
  
  return(dbcan_string)
}


#makes column with all relevant dbCAN id's for each EC in picrust output
makeDBCANcolumn <- function(df.picr=picr_ec_long, df.cazyid=cazyid_EC){
  n_iter <- nrow(df.picr)
  pb <- txtProgressBar(min = 0,
                       max = n_iter,
                       style = 3,
                       width = n_iter, 
                       char = "=") 
  
  init <- numeric(n_iter)
  end <- numeric(n_iter)
  
  out_column <- c()
  
  for(i in 1:nrow(df.picr)){
    init[i] <- Sys.time()
    
    #code
    out_column <- c(out_column, buildDBCANstring(df.picr=df.picr,
                                                 df.cazyid=df.cazyid,
                                                 row = i))
    
    end[i] <- Sys.time()
    
    #progress bar
    setTxtProgressBar(pb, i)
    time <- round(seconds_to_period(sum(end - init)), 0)
    est <- n_iter * (mean(end[end != 0] - init[init != 0])) - time
    remaining <- round(seconds_to_period(est), 0)
    cat(paste(" // Execution time:", time,
              " // Estimated time remaining:", remaining), "")
  }
  return(out_column)
}


#get list of all EC's in cazyid_EC
all_ECs <- c()
for(i in cazyid_EC$EC){
  all_ECs <- c(all_ECs, unlist(PullStr(i, type = 2)))
}
#remove duplicates
all_ECs_deduped <- unique(all_ECs)

#filter picr_ec to only include relevant EC's in all_EC vector
#remove any ASVs with 0 in every column

col_ids <- colnames(picr_ec)[-1]
cols_keep <- all_ECs_deduped[all_ECs_deduped %in% col_ids]

picr_ec_f<- picr_ec %>%
  select(sequence, all_of(cols_keep)) 
total <- rowSums(picr_ec_f[,-1]) 
picr_ec_f <- cbind(picr_ec_f, total)

picr_ec_f <- picr_ec_f %>%
  filter(total > 0) %>%
  select(-total)


picr_ec_long <- picr_ec_f %>%
  pivot_longer(cols = starts_with("EC:"),
               names_to = "EC",
               values_to = "count")

##make dbCAN column
dbcan_column <- makeDBCANcolumn()

#add dbcan_column to picr counts
picr_ec_dbcan <- cbind(picr_ec_long, dbcan_column)

write_tsv(picr_ec_dbcan, "Clean_Data/picr_ec_dbcan_day11.tsv")
