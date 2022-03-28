## ---------------------------
##
## Script: 23_call_functions_picr_Day11
##
## Purpose: Provides T/F for function presence in picrust2 predicted genomes
##
## Author: Ikaia Leleiwi
##
## Date Created: 2022-03-22
##
## Copyright (c) Ikaia Leleiwi, 2022
## Email: ileleiwi@gmail.com
##
## ---------------------------
##
## Notes: 
## Requires function rules tsv and properly formatted gene counts table
##   
##
## ---------------------------

## set working directory

setwd(paste0("/Users/ikaialeleiwi/Desktop/Lab/Salmonella_NIH/Lactobacillus/",
             "Omics/Metagenome/rRNA/16S_From_All_Bins/MQHQ_Picrust2/",
             "MQHQ_Picrust2/"))

## ---------------------------

#load libraries
library(tidyverse)

#load functions for function calling
source("Scripts/rule_functions.R")

#load data
picr_day11 <- read_tsv("Clean_Data/picr_day11.tsv") %>%
  as.data.frame()
rules <- read_tsv("Clean_Data/picr_rules.tsv")

#create checkRule function
checkRule <- makecheckRule(rules)

#call functions
picr_day11_functions <- evaluateCountsDf(picr_day11, rules, omit = "ETC50")

#write out function dataframe
write_tsv(picr_day11_functions, "Clean_Data/picr_day11_functions.tsv")
