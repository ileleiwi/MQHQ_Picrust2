library(tidyverse)
library(GGally)

setwd(paste0("/Users/ikaialeleiwi/Desktop/Lab/Salmonella_NIH/Lactobacillus/",
             "Omics/Metagenome/rRNA/16S_From_All_Bins/MQHQ_Picrust2/",
             "MQHQ_Picrust2/"))

#Data
picr_dram <- read_tsv("Clean_Data/picr_dram.tsv")

#remove functions with 0 in both picrust and dram
picr_dram <- picr_dram %>%
  filter(counts_picrust > 0 & counts_dram > 0)

##Correlations
pcor <- cor(picr_dram$counts_dram, 
            picr_dram$counts_picrust,
            method = "pearson")

scor <- cor(picr_dram$counts_dram, 
            picr_dram$counts_picrust,
            method = "spearman")


#####do this with function groups
pdk_mtx <- picr_dram %>%
  select(starts_with("counts")) %>%
  as.matrix() %>%
  t()
colnames(pdk_mtx) <- picr_dram$gene_id

pdk_mtx %>% 
  ggcorr(method = c("everything", "pearson")) 

