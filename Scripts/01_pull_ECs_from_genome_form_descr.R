library(tidyverse)

setwd(paste0("/Users/ikaialeleiwi/Desktop/Lab/Salmonella_NIH/Lactobacillus/",
             "Omics/Metagenome/rRNA/16S_From_All_Bins/MQHQ_Picrust2/",
             "MQHQ_Picrust2/"))

##Read in data

#DRAM genome summary form
genome_form_cazy <- read_tsv("Data/genome_summary_form.tsv") %>%
  filter(header == "CAZY") 

##Create df of dbCAN/EC matches fromm genome_form

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

#deduplicates list elements and returns character vector of unique values
getUniqueECs <- function(list_in){
  un <- unlist(list_in)
  return(unique(un))
}

#function to handle wildcard in EC number
wcHandle <- function(vect_in){
  vect_out <- c()
  for(i in vect_in){
    if(!str_detect(i, "-")){
      vect_out <- c(vect_out, i)
    }else{
      next
    }
  }
  return(vect_out)
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

#Create dataframe
MakeCazyECdf <- function(){
  cazyid_EC <- data.frame(cazyid = character(),
                          EC = character())
  for (i in 1:nrow(genome_form_cazy)){
    ls <- PullStr(genome_form_cazy[i,2])
    ECs <- getUniqueECs(ls)
    ECs_no_wc <- wcHandle(ECs)
    EC_string <- collapseVect(ECs_no_wc)
    if(length(ECs_no_wc) > 0){
      gene_id <- genome_form_cazy %>%
        slice(i) %>%
        pull(gene_id)
      df <- data.frame(cazyid = gene_id,
                       EC = EC_string)
      cazyid_EC <- rbind(cazyid_EC, df)
  
    } else {
      next
    }
  }
  return(cazyid_EC)
}

cazyid_EC <- MakeCazyECdf()

write_tsv(cazyid_EC, "Clean_Data/cazyid_EC.tsv")
