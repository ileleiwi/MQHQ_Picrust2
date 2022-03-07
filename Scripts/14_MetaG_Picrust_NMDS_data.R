library(tidyverse)
library(readxl)
library(vegan)

setwd(paste0("/Users/ikaialeleiwi/Desktop/Lab/Salmonella_NIH/Lactobacillus/",
             "Omics/Metagenome/rRNA/16S_From_All_Bins/MQHQ_Picrust2/",
             "MQHQ_Picrust2/"))

#Data
picr_dram_ko <- read_tsv("Clean_Data/picr_dram_ko.tsv")

#remove functions with 0 in both picrust and dram
picr_dram_ko <- picr_dram_ko %>%
  filter(total_copy_number_picrust > 0 & total_copy_number_dram > 0)

#dataframe for nmds
picr_dram_long <- picr_dram %>%
  pivot_longer(cols = starts_with("total"),
               values_to = "total_copy_number",
               names_to = "method") %>%
  mutate(method = ifelse(str_detect(method, "dram"),
                         "DRAM",
                         "Picrust2"),
         unique_name = paste(gene_id, method, sep = "_")) 


picr_dram_matrix <- picr_dram_long %>%
  ungroup() %>%
  select(unique_name, bin.id, total_copy_number) %>%
  pivot_wider(names_from = "unique_name",
              values_from = "total_copy_number") %>%
  mutate(
    across(everything(), ~replace_na(.x, 0))
  ) %>%
  column_to_rownames(var = "bin.id") %>%
  as.matrix() %>%
  t()

ord_bray <- metaMDS(picr_dram_matrix, 
                    distance = "bray", 
                    autotransform = FALSE, 
                    noshare = 0.1, 
                    trace = 1)

# #getting scores for NMDS object and saving as df
ord_bray_scrs <- as.data.frame(scores(ord_bray),display="sites")

#make plot dataframe
nmds_plot_df <- ord_bray_scrs %>%
  rownames_to_column(var = "unique_id") %>%
  separate(unique_id,
           c("gene_id", "method"),
           sep = "_")
#perform Shepards goodness test on NMDS and anosim on control vs salmonella
Shepards_goodness_test_results_dram <- goodness(ord_bray)

anosim_results <- anosim(picr_dram_matrix, 
                         nmds_plot_df$method, 
                         permutations = 999, 
                         distance = "bray", 
                         strata = NULL,
                         parallel = getOption("mc.cores"))

stress_plot <- stressplot(ord_bray)

#mrpp
#getting grouping information
grouping_variable <- nmds_plot_df$method

mrpp_out <- mrpp(picr_dram_matrix,
                 grouping = grouping_variable,
                 distance = "bray")

write_tsv(nmds_plot_df, "Clean_Data/nmds_plot_df.tsv")
