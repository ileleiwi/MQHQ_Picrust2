library(tidyverse)
library(readxl)

setwd(paste0("/Users/ikaialeleiwi/Desktop/Lab/Salmonella_NIH/Lactobacillus/",
             "Omics/Metagenome/rRNA/16S_From_All_Bins/MQHQ_Picrust2/",
             "MQHQ_Picrust2/"))


#Data 
resp_rules <- read_tsv("Data/Respiration_rules.tsv") %>%
  as.data.frame()

ko_rules <- read_xlsx("Data/KOs_for_functions.xlsx") %>%
  select(1:3)

#make nested list

#get function name
getName <- function(df, row, column){
  name <- df %>%
    slice(row) %>%
    pull(column)
  return(name)
}

#get the main rule for each function
pullMainRule <- function(df){
  list_parent <- list()
  for(i in 1:nrow(df)){
    if(df[i,1] == df[i,2]){
      x <- list(main_rule = getName(df, i, 3))
      names(x) <- getName(df, i, 1)
      list_parent <- append(list_parent, x)
    }
  }
  return(list_parent)
}


main_rules <- pullMainRule(resp_rules)

#get subrules
pullSubRules <- function(df){
  list_sub <- list()
  
  for(i in 2:nrow(df)){
    if(df[i,1] != df[i,2] & df[i-1,1] == df[i,1]){
      x <- list(sub_rule = getName(df, i, 3))
      names(x) <- getName(df,i,2)
      list_sub <- append(list_sub, x)
    }
  }
  return(list_sub)
}

sub_rules <- pullSubRules(resp_rules)
