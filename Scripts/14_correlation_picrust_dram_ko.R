library(tidyverse)
library(GGally)

setwd(paste0("/Users/ikaialeleiwi/Desktop/Lab/Salmonella_NIH/Lactobacillus/",
             "Omics/Metagenome/rRNA/16S_From_All_Bins/ASVs_to_Bins_MMSEQS2/"))

#Data
picr_dram_ko <- read_tsv("Clean_Data/picr_dram_ko.tsv")

#remove functions with 0 in both picrust and dram
picr_dram_ko <- picr_dram_ko %>%
  filter(total_copy_number_picrust > 0 & total_copy_number_dram > 0)

##Correlations
pcor <- cor(picr_dram_ko$total_copy_number_dram, 
            picr_dram_ko$total_copy_number_picrust,
            method = "pearson")

scor <- cor(picr_dram_ko$total_copy_number_dram, 
            picr_dram_ko$total_copy_number_picrust,
            method = "spearman")


#####do this with function groups
pdk_mtx <- picr_dram_ko %>%
  select(starts_with("Total")) %>%
  as.matrix() %>%
  t()
colnames(pdk_mtx) <- picr_dram_ko$gene_id

pdk_mtx %>% 
  ggcorr(method = c("everything", "pearson")) 

