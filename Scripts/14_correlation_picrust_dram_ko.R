library(tidyverse)
library(GGally)
library(ggcorrplot)

setwd(paste0("/Users/ikaialeleiwi/Desktop/Lab/Salmonella_NIH/Lactobacillus/",
             "Omics/Metagenome/rRNA/16S_From_All_Bins/MQHQ_Picrust2/",
             "MQHQ_Picrust2/"))

#Data
picr_dram <- read_tsv("Clean_Data/function_pa_dram_picrust.tsv")


##Correlations
pcor <- cor(picr_dram$dr_prop, 
            picr_dram$pi_prop,
            method = "pearson")

scor <- cor(picr_dram$dr_prop, 
            picr_dram$pi_prop,
            method = "spearman")


##correlation matrix
corr <- round(cor(pdk_mtx), 1)
head(corr[, 1:6])

#####do this with function groups
pdk_mtx <- picr_dram %>%
  select(ends_with("prop")) %>%
  as.matrix() %>%
  t()
colnames(pdk_mtx) <- picr_dram$functions

pdk_mtx %>% 
  ggcorr(method = c("everything", "spearman")) 

