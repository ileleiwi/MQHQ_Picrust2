library(tidyverse)
library(readxl)

setwd(paste0("/Users/ikaialeleiwi/Desktop/Lab/Salmonella_NIH/Lactobacillus/",
             "Omics/Metagenome/rRNA/16S_From_All_Bins/MQHQ_Picrust2/",
             "MQHQ_Picrust2/"))


#Data 
picr_dram <- read_tsv("Clean_Data/picr_dram_ko.tsv")

ko_rules <- read_xlsx("Data/KOs_for_functions.xlsx") %>%
  select(1:3)

resp_rules <- read_tsv("Data/Respiration_rules.tsv")


#functions to parse through resp_rules

  

