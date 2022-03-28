## ---------------------------
##
## Script: 33_day11_16S_relative_abundace
##
## Purpose: calculate relative abundance for asvs from day11 samples
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
## Notes: requires feature table filtered to day11 samples
##   
##
## ---------------------------

## set working directory

setwd("/Users/ikaialeleiwi/Desktop/Lab/Salmonella_NIH/Lactobacillus/",
      "Omics/Metagenome/rRNA/16S_From_All_Bins/MQHQ_Picrust2/",
      "MQHQ_Picrust2/")

## ---------------------------

library(tidyverse)
library(edgeR) #GeTMM
library(compositions) #clr
library(ComplexHeatmap) #heatmap

#functions
relabund <- function(df, columns = c(NA)) 
  #takes feature table and calculates relative abundance of columns, omitting NA's
  #considers only columns listed in columns argument (character vector of column names). 
  #columns (default) = all columns
{
  if (NA %in% columns){
    df <- sweep(df, 2, colSums(df, na.rm = TRUE), '/')
  }
  else {
    df_relabund <- df %>% select(all_of(columns)) %>%
      sweep(., 2, colSums(., na.rm = TRUE), '/')
    df[,columns] <- df_relabund
  }
  
  return(df)
}

#read in data
ft_day11 <- read_tsv("Clean_Data/feature_table_11.tsv")


