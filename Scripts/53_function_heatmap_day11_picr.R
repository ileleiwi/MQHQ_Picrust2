## ---------------------------
##
## Script: 53_function_heatmap_day11_picr
##
## Purpose: Make heatmap of picrust2 predicted genome functions from 16S relative abundance 
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
## Notes: requires T/F function dataframe and feature_table from qiime2 and taxonomy from qiime2
## 
## Summed relative abundance of each ASV in each sample with given function normalized within each function
##   
##
## ---------------------------

## set working directory

setwd(paste0("/Users/ikaialeleiwi/Desktop/Lab/Salmonella_NIH/Lactobacillus/",
             "Omics/Metagenome/rRNA/16S_From_All_Bins/MQHQ_Picrust2/",
             "MQHQ_Picrust2/"))

## ---------------------------

library(tidyverse)
library(ComplexHeatmap) #heatmap
library(circlize) #colorRamp2
library(ggpubr)
library(cowplot)
library(RColorBrewer)


#read in data
metaG_samples <- c("13-(11)-S-F-KE", "20-(11)-S-F-KL", "27-(11)-S-F-KS", 
                   "49-(11)-C-F-LO", "47-(11)-C-F-LM", "50-(11)-C-F-LP")

design <- list(control = c("49-(11)-C-F-LO", "47-(11)-C-F-LM", "50-(11)-C-F-LP"),
               infected = c("13-(11)-S-F-KE", "20-(11)-S-F-KL", "27-(11)-S-F-KS"))

taxonomy <- read_tsv("Clean_Data/taxonomy.tsv") %>%
  rename("FeatureID" = "Feature ID")

dram_product <- read_tsv("Clean_Data/dram_product_picrust_day11.tsv")

picr_func <- read_tsv("Clean_Data/picr_day11_functions.tsv") %>%
  left_join(taxonomy, by = c("genome" = "FeatureID")) %>%
  filter(!str_detect(Taxon, "D_4__Mitochondria|D_3__Chloroplast")) %>%# remove mitochondria and chloroplast from picrust2
  select(-Taxon, -Confidence) %>%
  left_join(dram_product, by = c("genome" = "bins")) %>%
  column_to_rownames(var = "genome")
  
  
ft <- read_tsv("Clean_Data/feature_table_clean_r1-r5.tsv") %>%
  select(asv_id, all_of(metaG_samples)) 

################functions################

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


FillSampleMatrices <- function(sample, function_df, feature_table){
  #Function that makes df for each sample filled with 
  #relative abundance of mapped reads from each function present
  #sample = column name of feature table
  for(asv_id in rownames(function_df)){
    for(funct in colnames(function_df)){
      ifelse(function_df[asv_id,funct],
             function_df[asv_id,funct] <- feature_table[asv_id, sample],
             function_df[asv_id,funct] <- 0)
    }
  }
  return(function_df)
}

generateRelabundFunctList <- function(samples = colnames(feature_table), 
                                      function_df, 
                                      feature_table){
  #Creates list of relabund dfs for each sample
  #samples = column names of feature table (samples)
  #function_df = T/F dataframe with functions as columns and asvs/genomes as rows
  #feature_table = qimme2 feature table (samples as columns rownames as asv_ids)
  working_list <- list()
  for(samp in samples){
    working_list[[length(working_list) + 1]] <- FillSampleMatrices(samp, 
                                                                   function_df, 
                                                                   feature_table)
  }
  names(working_list) <- samples
  return(working_list)
}


SumAcrossCols <- function(samples_list){
  #Sum relative abundances within each sample for each column (function)
  for(sample in names(samples_list)){
    samples_list[[sample]] <- apply(samples_list[[sample]], 2, sum)
  }
  return(samples_list)
}


MeanWithinTrt <- function(sample_function_list, design){
  #Mean relative abundances within each treatment for each row (function)
  #sample_function_list = list of sample_function dataframes
  #design = named list of treatment samples and controls
  working_list <- list()
  for(trt in 1:length(design)){
    trts <- sample_function_list[design[[trt]]]
    working_list[[length(working_list)+1]] <- Reduce('+', trts)
    working_list[[length(working_list)]] <- map_dbl(working_list[[length(working_list)]], ~.x/length(trts))
  }
  names(working_list) <- names(design)
  return(working_list)
}

SumWithinTrt <- function(sample_function_list, design){
  #Sum relative abundances within each treatment for each row (function)
  #sample_function_list = list of sample_function dataframes
  #design = named list of treatment samples and controls
  working_list <- list()
  for(trt in 1:length(design)){
    trts <- sample_function_list[design[[trt]]]
    working_list[[length(working_list)+1]] <- t( Reduce('+', trts))
  }
  names(working_list) <- names(design)

  return(working_list)
}


CountAcrossCols <- function(working_list){
  #count number of genomes within each sample that are more than 0.5% relative abundance
  #for each function
  for(i in names(working_list)){
    sample_matrix <- ifelse(working_list[[i]] >= 0.5, 1, 0)
    sample_matrix <- apply(sample_matrix, 2, sum)
    sample_df <- as.data.frame(sample_matrix)
    colnames(sample_df) <- i
    working_list[[i]] <- sample_df
  }
  return(working_list)
}
################^functions^################

#calculate relative abundance
ft_rel <- ft %>%
  column_to_rownames(var = "asv_id") %>%
  relabund() %>%
  mutate_all(~.*100) %>%
  rownames_to_column(var = "asv_id") %>%
  filter(asv_id %in% rownames(picr_func)) %>% #filter feature table to picr2 predicted asvs
  column_to_rownames(var = "asv_id")

#filter function dataframe to ASVs present in feature table
picr_func <- picr_func[rownames(ft_rel),] 

#create list of summed relative abundance for each sample and function
summed_relabund <- generateRelabundFunctList(function_df = picr_func, feature_table = ft_rel) %>%
  SumAcrossCols()

#create dataframe where columns are sample names and rows are functions
samples_by_function <- do.call(cbind.data.frame, summed_relabund)

#modify rownames
samples_by_function <- samples_by_function %>%
  rownames_to_column(var = "functions") %>%
  mutate(functions = case_when(str_detect(functions, "NonSpecificHexoseSugar") ~ "Hexose",
                               str_detect(functions, "GalacturonicAcid") ~ "Galacturonic Acid",
                               str_detect(functions, "Non.SpecificAlcohols") ~ "Alcohol",
                               str_detect(functions, "SulfateReduction") ~ "Sulfate Reduction",
                               str_detect(functions, "FumarateReduction") ~ "Fumarate Reduction",
                               str_detect(functions, "TetrathionateReduction") ~ "Tetrathionate Reduction",
                               str_detect(functions, "DenitrificationNotETC") ~ "Denitrification (Not ETC)",
                               str_detect(functions, "Alpha-galactans") ~ "Alpha-Galactans",
                               str_detect(functions, "Beta-galactan") ~ "Beta-Galactan (Pectic Galactan)",
                               str_detect(functions, "Mixed-Linkage glucans") ~ "Mixed-Linkage Glucans",
                               str_detect(functions, "Beta-mannan") ~ "Beta-Mannan",
                               str_detect(functions, "Alpha-mannan") ~ "Alpha-Mannan",
                               str_detect(functions, "Rhamnose") ~ "Rhamnose Cleavage",
                               str_detect(functions, "Arabinose") ~ "Arabinose Cleavage",
                               functions == "Fucose" ~ "Fucose Cleavage",
                               T ~ functions),
         functions = str_remove(functions, "Degradation")) %>%
  column_to_rownames(var = "functions")

#reorder samples_by_function
polymers <- c("Arabinan", "Pectin", "Arabinose Cleavage", "Polyphenolics", "Xyloglucan",
              "Amorphous Cellulose", "Xylans", "Mixed-Linkage Glucans", 
              "Crystalline Cellulose", "Alpha-Galactans", "Chitin", 
              "Beta-Galactan (Pectic Galactan)", "Alpha-Mannan", "Starch", 
              "Fucose Cleavage", "Rhamnose Cleavage", "Sulf-Polysachharides", "Mucin",
              "Beta-Mannan")


sugars <- c("Galactose", "Galacturonic Acid", "Lactose", "Hexose", "Fructose",
            "Mannose", "Fucose", "Sucrose", "Xylose")

scfas <- c("Acetate", "Butyrate", "Propionate", "Alcohol", "Lactate")

resp <- c("Denitrification (Not ETC)", "Sulfate Reduction", 
          "Fumarate Reduction", "TMAO", "Aerobic", "Tetrathionate Reduction",
          "Microaerophillic")
new_order <- c(polymers, sugars, scfas, resp)

samples_by_function <- samples_by_function[match(new_order,
                                                 rownames(samples_by_function)),]

#scale by function
samples_by_function_scaled <- t(scale(t(samples_by_function)))
samples_by_function_scaled[is.nan(samples_by_function_scaled)] <- 0
samples_by_function_scaled <- as.data.frame(samples_by_function_scaled)

#order samples
samples_by_function_scaled <- samples_by_function_scaled[,c(6,4,5,1,2,3)] 

#write out for stats script
samples_by_function_out <- samples_by_function %>%
  rownames_to_column(var = "function")
write_tsv(samples_by_function_out, "Clean_Data/samples_by_function_out_day11.tsv")

#run stats on samples by function
source("Scripts/43_stats_for_function_heatmap_day11.R")

#create dataframe where columns are sample names and rows are functions for healthy and infected
treatment_by_function <- SumWithinTrt(generateRelabundFunctList(function_df = picr_func, feature_table = ft_rel), design)
#combine treatments into one dataframe with samples labeled from which treatment they originated
treatment_by_function <- do.call(cbind.data.frame, treatment_by_function)

treatment_by_function <- treatment_by_function %>%
  rownames_to_column(var = "functions") %>%
  pivot_longer(cols = -functions,
               names_to = "trt_bin",
               values_to = "summed_relative_abundance") %>%
  mutate(treatment = case_when(str_starts(trt_bin, "control.") ~ "Control",
                               str_starts(trt_bin, "infected.") ~ "Infected"),
         trt_bin = str_remove(trt_bin, "control.|infected."),
         functions = as.factor(functions),
         treatment = as.factor(treatment),
         trt_bin = as.factor(trt_bin)) %>%
  rename("bins" = "trt_bin") %>%
  select(bins, treatment, functions, summed_relative_abundance) %>%
  left_join(taxonomy, by = c("bins" = "FeatureID"))


#Split dataframes by healthy and infected for plotting
healthy <- samples_by_function_scaled %>%
  rownames_to_column(var = "functions") %>%
  select(matches("LP|LO|LM|functions")) %>%
  arrange(match(functions, c("Polyphenolics", "Crystalline Cellulose","Amorphous Cellulose", 
                             "Mixed-Linkage Glucans","Xylans","Xyloglucan", "Alpha-Mannan", 
                             "Beta-Mannan", "Mucin", "Starch", "Pectin", "Chitin", "Arabinan",
                             "Alpha-Galactans", "Beta-Galactan (Pectic Galactan)",
                             "Sulf-Polysachharides", "Arabinose Cleavage","Rhamnose Cleavage",
                             "Fucose Cleavage", "Lactose", "Sucrose", "Galactose", "Mannose", 
                             "Galacturonic Acid", "Hexose", "Xylose", "Fucose", 
                             "Fructose","Butyrate", "Propionate", "Acetate", "Lactate", 
                             "Alcohol","Sulfate Reduction", "Fumarate Reduction", "Aerobic", 
                             "Tetrathionate Reduction", "Microaerophillic",
                             "Denitrification (Not ETC)", "TMAO"))) %>%
  column_to_rownames(var = "functions")


infected <- samples_by_function_scaled %>%
  rownames_to_column(var = "functions") %>%
  select(matches("KE|KL|KS|functions")) %>%
  arrange(match(functions, c("Polyphenolics", "Crystalline Cellulose","Amorphous Cellulose", 
                             "Mixed-Linkage Glucans","Xylans","Xyloglucan", "Alpha-Mannan", 
                             "Beta-Mannan", "Mucin", "Starch", "Pectin", "Chitin", "Arabinan",
                             "Alpha-Galactans", "Beta-Galactan (Pectic Galactan)",
                             "Sulf-Polysachharides", "Arabinose Cleavage","Rhamnose Cleavage",
                             "Fucose Cleavage", "Lactose", "Sucrose", "Galactose", "Mannose", 
                             "Galacturonic Acid", "Hexose", "Xylose", "Fucose", 
                             "Fructose","Butyrate", "Propionate", "Acetate", "Lactate", 
                             "Alcohol","Sulfate Reduction", "Fumarate Reduction", "Aerobic", 
                             "Tetrathionate Reduction", "Microaerophillic",
                             "Denitrification (Not ETC)", "TMAO"))) %>%
  column_to_rownames(var = "functions")

healthy_unscaled <- samples_by_function %>%
  select(matches("LM|LO|LP")) %>%
  rownames_to_column(var = "functions") %>%
  arrange(match(functions, c("Polyphenolics", "Crystalline Cellulose","Amorphous Cellulose", 
                             "Mixed-Linkage Glucans","Xylans","Xyloglucan", "Alpha-Mannan", 
                             "Beta-Mannan", "Mucin", "Starch", "Pectin", "Chitin", "Arabinan",
                             "Alpha-Galactans", "Beta-Galactan (Pectic Galactan)",
                             "Sulf-Polysachharides", "Arabinose Cleavage","Rhamnose Cleavage",
                             "Fucose Cleavage", "Lactose", "Sucrose", "Galactose", "Mannose", 
                             "Galacturonic Acid", "Hexose", "Xylose", "Fucose", 
                             "Fructose","Butyrate", "Propionate", "Acetate", "Lactate", 
                             "Alcohol","Sulfate Reduction", "Fumarate Reduction", "Aerobic", 
                             "Tetrathionate Reduction", "Microaerophillic",
                             "Denitrification (Not ETC)", "TMAO"))) %>%
  column_to_rownames(var = "functions")

infected_unscaled <- samples_by_function %>%
  select(matches("KE|KL|KS")) %>%
  rownames_to_column(var = "functions") %>%
  arrange(match(functions, c("Polyphenolics", "Crystalline Cellulose","Amorphous Cellulose", 
                             "Mixed-Linkage Glucans","Xylans","Xyloglucan", "Alpha-Mannan", 
                             "Beta-Mannan", "Mucin", "Starch", "Pectin", "Chitin", "Arabinan",
                             "Alpha-Galactans", "Beta-Galactan (Pectic Galactan)",
                             "Sulf-Polysachharides", "Arabinose Cleavage","Rhamnose Cleavage",
                             "Fucose Cleavage", "Lactose", "Sucrose", "Galactose", "Mannose", 
                             "Galacturonic Acid", "Hexose", "Xylose", "Fucose", 
                             "Fructose","Butyrate", "Propionate", "Acetate", "Lactate", 
                             "Alcohol","Sulfate Reduction", "Fumarate Reduction", "Aerobic", 
                             "Tetrathionate Reduction", "Microaerophillic",
                             "Denitrification (Not ETC)", "TMAO"))) %>%
  column_to_rownames(var = "functions")


#create list of number of genomes more than 0.5% relative abundancefor each sample and function
genome_count <- CountAcrossCols(generateRelabundFunctList(function_df = picr_func, feature_table = ft_rel))

genome_count_df <- do.call(cbind.data.frame, genome_count)

#modify rownames
genome_count_df <- genome_count_df %>%
  rownames_to_column(var = "functions") %>%
  mutate(functions = case_when(str_detect(functions, "NonSpecificHexoseSugar") ~ "Hexose",
                               str_detect(functions, "GalacturonicAcid") ~ "Galacturonic Acid",
                               str_detect(functions, "Non.SpecificAlcohols") ~ "Alcohol",
                               str_detect(functions, "SulfateReduction") ~ "Sulfate Reduction",
                               str_detect(functions, "FumarateReduction") ~ "Fumarate Reduction",
                               str_detect(functions, "TetrathionateReduction") ~ "Tetrathionate Reduction",
                               str_detect(functions, "DenitrificationNotETC") ~ "Denitrification (Not ETC)",
                               str_detect(functions, "Alpha-galactans") ~ "Alpha-Galactans",
                               str_detect(functions, "Beta-galactan") ~ "Beta-Galactan (Pectic Galactan)",
                               str_detect(functions, "Mixed-Linkage glucans") ~ "Mixed-Linkage Glucans",
                               str_detect(functions, "Beta-mannan") ~ "Beta-Mannan",
                               str_detect(functions, "Alpha-mannan") ~ "Alpha-Mannan",
                               str_detect(functions, "Rhamnose") ~ "Rhamnose Cleavage",
                               str_detect(functions, "Arabinose") ~ "Arabinose Cleavage",
                               functions == "Fucose" ~ "Fucose Cleavage",
                               T ~ functions),
         functions = str_remove(functions, "Degradation")) %>%
  column_to_rownames(var = "functions")

#dataframes for heatmap barplot annotations
genomes_healthy <- genome_count_df %>%
  select(matches("LM|LO|LP")) %>%
  mutate(healthy = rowSums(.)) %>%
  rownames_to_column(var = "functions") %>%
  arrange(match(functions, c("Polyphenolics", "Crystalline Cellulose","Amorphous Cellulose", 
                             "Mixed-Linkage Glucans","Xylans","Xyloglucan", "Alpha-Mannan", 
                             "Beta-Mannan", "Mucin", "Starch", "Pectin", "Chitin", "Arabinan",
                             "Alpha-Galactans", "Beta-Galactan (Pectic Galactan)",
                             "Sulf-Polysachharides", "Arabinose Cleavage","Rhamnose Cleavage",
                             "Fucose Cleavage", "Lactose", "Sucrose", "Galactose", "Mannose", 
                             "Galacturonic Acid", "Hexose", "Xylose", "Fucose", 
                             "Fructose","Butyrate", "Propionate", "Acetate", "Lactate", 
                             "Alcohol","Sulfate Reduction", "Fumarate Reduction", "Aerobic", 
                             "Tetrathionate Reduction", "Microaerophillic",
                             "Denitrification (Not ETC)", "TMAO"))) %>%
  select(functions, healthy) %>%
  column_to_rownames(var = "functions")

genomes_infected <- genome_count_df %>%
  select(matches("KE|KL|KS")) %>%
  mutate(infected = rowSums(.)) %>%
  rownames_to_column(var = "functions") %>%
  arrange(match(functions, c("Polyphenolics", "Crystalline Cellulose","Amorphous Cellulose", 
                             "Mixed-Linkage Glucans","Xylans","Xyloglucan", "Alpha-Mannan", 
                             "Beta-Mannan", "Mucin", "Starch", "Pectin", "Chitin", "Arabinan",
                             "Alpha-Galactans", "Beta-Galactan (Pectic Galactan)",
                             "Sulf-Polysachharides", "Arabinose Cleavage","Rhamnose Cleavage",
                             "Fucose Cleavage", "Lactose", "Sucrose", "Galactose", "Mannose", 
                             "Galacturonic Acid", "Hexose", "Xylose", "Fucose", 
                             "Fructose","Butyrate", "Propionate", "Acetate", "Lactate", 
                             "Alcohol","Sulfate Reduction", "Fumarate Reduction", "Aerobic", 
                             "Tetrathionate Reduction", "Microaerophillic",
                             "Denitrification (Not ETC)", "TMAO"))) %>%
  select(functions, infected) %>%
  column_to_rownames(var = "functions")


#Split dataframes by function type for plotting
cazy_sbf <- samples_by_function_scaled %>%
  rownames_to_column(var = "functions") %>%
  filter(functions %in% c("Polyphenolics", "Crystalline Cellulose","Amorphous Cellulose", 
                          "Mixed-Linkage Glucans","Xylans","Xyloglucan", "Alpha-Mannan", 
                          "Beta-Mannan", "Mucin", "Starch", "Pectin", "Chitin", "Arabinan",
                          "Alpha-Galactans", "Beta-Galactan (Pectic Galactan)",
                          "Sulf-Polysachharides", "Arabinose Cleavage","Rhamnose Cleavage",
                          "Fucose Cleavage")) %>%
  arrange(match(functions, c("Polyphenolics", "Crystalline Cellulose","Amorphous Cellulose", 
                             "Mixed-Linkage Glucans","Xylans","Xyloglucan", "Alpha-Mannan", 
                             "Beta-Mannan", "Mucin", "Starch", "Pectin", "Chitin", "Arabinan",
                             "Alpha-Galactans", "Beta-Galactan (Pectic Galactan)",
                             "Sulf-Polysachharides", "Arabinose Cleavage","Rhamnose Cleavage",
                             "Fucose Cleavage"))) %>%
  column_to_rownames(var = "functions")

sug_sbf <- samples_by_function_scaled %>%
  rownames_to_column(var = "functions") %>%
  filter(functions %in% c("Lactose", "Sucrose", "Galactose", "Mannose", 
                          "Galacturonic Acid", "Hexose", "Xylose", "Fucose", 
                          "Fructose")) %>%
  arrange(match(functions, c("Lactose", "Sucrose", "Galactose", "Mannose", 
                             "Galacturonic Acid", "Hexose", "Xylose", "Fucose", 
                             "Fructose"))) %>%
  column_to_rownames(var = "functions")

scfa_sbf <- samples_by_function_scaled %>%
  rownames_to_column(var = "functions") %>%
  filter(functions %in% c("Butyrate", "Propionate", "Acetate", "Lactate", 
                          "Alcohol")) %>%
  arrange(match(functions, c("Butyrate", "Propionate", "Acetate", "Lactate", 
                             "Alcohol"))) %>%
  column_to_rownames(var = "functions")

resp_sbf <- samples_by_function_scaled %>%
  rownames_to_column(var = "functions") %>%
  filter(functions %in% c("Sulfate Reduction", "Fumarate Reduction",
                          "Aerobic", "Tetrathionate Reduction", "Microaerophillic",
                          "Denitrification (Not ETC)", "TMAO")) %>%
  arrange(match(functions, c("Sulfate Reduction", "Fumarate Reduction",
                             "Aerobic", "Tetrathionate Reduction", "Microaerophillic",
                             "Denitrification (Not ETC)", "TMAO"))) %>%
  column_to_rownames(var = "functions")

## plots
healthy_annotation <- rowAnnotation(genomes = anno_barplot(genomes_healthy$healthy,
                                                           which = "row",
                                                           border = FALSE,
                                                           axis_param = list(direction = "reverse"),
                                                           ylim = c(0,60),
                                                           height = unit(8, "cm")),
                                    annotation_label = "# of Genomes\n(RelAbund >=0.5%)",
                                    annotation_name_gp= gpar(fontsize = 8))

infected_annotation <- rowAnnotation(genomes = anno_barplot(genomes_infected$infected,
                                                            which = "row",
                                                            border = FALSE,
                                                            axis_param = list(direction = "reverse"),
                                                            ylim = c(0,60),
                                                            height = unit(8, "cm")),
                                     annotation_label = "# of Genomes\n(RelAbund >=0.5%)",
                                     annotation_name_gp= gpar(fontsize = 8))

ht_healthy <- Heatmap(as.matrix(healthy),
                      col = c("#2A3479", "#FCF8CC"),
                      border = TRUE,
                      row_split = factor(c(rep("Polymers", 19),
                                           rep("Sugars", 9),
                                           rep("Alcohol &\nFatty Acids", 5),
                                           rep("Respiration", 7)), 
                                         levels = c("Polymers", "Sugars",
                                                    "Alcohol &\nFatty Acids",
                                                    "Respiration")),
                      column_split = rep("Healthy", 3),
                      column_gap = unit(2, "mm"),
                      row_gap = unit(4, "mm"),
                      rect_gp = gpar(col = "black", lwd = 2),
                      cluster_rows = TRUE,
                      cluster_columns = TRUE,
                      cluster_column_slices = FALSE,
                      cluster_row_slices = FALSE,
                      clustering_method_rows = "complete",
                      left_annotation = healthy_annotation)

ht_infected <- Heatmap(as.matrix(infected),
                       col = c("#2A3479", "#FCF8CC"),
                       border = TRUE,
                       row_split = factor(c(rep("Polymers", 19),
                                            rep("Sugars", 9),
                                            rep("Alcohol &\nFatty Acids", 5),
                                            rep("Respiration", 7)), 
                                          levels = c("Polymers", "Sugars",
                                                     "Alcohol &\nFatty Acids",
                                                     "Respiration")),
                       column_split = rep("Infected", 3),
                       column_gap = unit(2, "mm"),
                       row_gap = unit(4, "mm"),
                       rect_gp = gpar(col = "black", lwd = 2),
                       row_names_gp = gpar(fontsize = 8), 
                       row_names_side = "right", 
                       row_title_side = "right",
                       cluster_rows = TRUE,
                       cluster_columns = TRUE,
                       cluster_column_slices = FALSE,
                       cluster_row_slices = FALSE,
                       clustering_method_rows = "complete",
                       left_annotation = infected_annotation)

ht_both <- ht_healthy + ht_infected

#pdf("day11_polymer_heatmap_scaled_largest_to_smallest.pdf", width = 10, height = 10)
ht_both
#dev.off()


## Heatmap Scaled by Row, Clustered hclust complete
genomes <- cbind(genomes_healthy, genomes_infected) %>%
  rownames_to_column(var = "functions")

polymers <- c("Arabinan", "Pectin", "Arabinose Cleavage", "Polyphenolics", "Xyloglucan",
  "Amorphous Cellulose", "Xylans", "Mixed-Linkage Glucans", 
  "Crystalline Cellulose", "Alpha-Galactans", "Chitin", 
  "Beta-Galactan (Pectic Galactan)", "Alpha-Mannan", "Starch", 
  "Fucose Cleavage", "Rhamnose Cleavage", "Sulf-Polysachharides", "Mucin",
  "Beta-Mannan")
cazy_g <- genomes[match(polymers, genomes$functions),] 
rownames(cazy_g) <- cazy_g$functions
cazy_g <- cazy_g %>%
  select(-functions)

sugars <- c("Galactose", "Galacturonic Acid", "Lactose", "Hexose", "Fructose",
            "Mannose", "Fucose", "Sucrose", "Xylose")
sugar_g <- genomes[match(sugars, genomes$functions),] 
rownames(sugar_g) <- sugar_g$functions
sugar_g <- sugar_g %>%
  select(-functions)

scfas <- c("Acetate", "Butyrate", "Propionate", "Alcohol", "Lactate")
scfa_g <- genomes[match(scfas, genomes$functions),] 
rownames(scfa_g) <- scfa_g$functions
scfa_g <- scfa_g %>%
  select(-functions)

resp <- c("Denitrification (Not ETC)", "Sulfate Reduction", 
          "Fumarate Reduction", "TMAO", "Aerobic", "Tetrathionate Reduction",
          "Microaerophillic")
resp_g <- genomes[match(resp, genomes$functions),] 
rownames(resp_g) <- resp_g$functions
resp_g <- resp_g %>%
  select(-functions)

genomes_reordered <- rbind(cazy_g, sugar_g, scfa_g, resp_g)
index <- match(rownames(genomes_reordered), significance_df$functions)

sig_df_reorderd <- significance_df[index,] 

rownames(sig_df_reorderd) <- NULL
sig_df_reorderd <- sig_df_reorderd %>%
  column_to_rownames(var = "functions") %>%
  select(pval) %>%
  mutate(pval = if_else(condition = pval <= 0.05,
                        true = pval,
                        false = NaN))

colors_sig <- c("#E54F6D", "#22162B")

pvalue_col_fun = colorRamp2(c(0, 2, 3), c("#E54F6D", "#22162B", "black")) 

significance_annoatation <- rowAnnotation(`p-value >= 0.05` = anno_simple(sig_df_reorderd$pval, col = pvalue_col_fun))

healthy_annotation <- rowAnnotation(genomes = anno_barplot(genomes_reordered$healthy,
                                                           which = "row",
                                                           border = FALSE,
                                                           axis_param = list(direction = "reverse"),
                                                           ylim = c(0,60),
                                                           height = unit(8, "cm")),
                                    annotation_label = "# of Genomes\n(RelAbund >=0.5%)",
                                    annotation_name_gp= gpar(fontsize = 8))

infected_annotation <- rowAnnotation(genomes = anno_barplot(genomes_reordered$infected,
                                                            which = "row",
                                                            border = FALSE,
                                                            axis_param = list(direction = "reverse"),
                                                            ylim = c(0,60),
                                                            height = unit(8, "cm")),
                                     annotation_label = "# of Genomes\n(RelAbund >=0.5%)",
                                     annotation_name_gp= gpar(fontsize = 8))





ht_both_clust <- Heatmap(as.matrix(samples_by_function_scaled),
                         col = c("#2A3479", "#FCF8CC"),
                         border = TRUE,
                         row_split = factor(c(rep("Polymers", 19),
                                              rep("Sugars", 9),
                                              rep("Alcohol &\nFatty Acids", 5),
                                              rep("Respiration", 7)), 
                                            levels = c("Polymers", "Sugars",
                                                       "Alcohol &\nFatty Acids",
                                                       "Respiration")),
                         column_split = c(rep("Healthy", 3), rep("Infected", 3)),
                         column_gap = unit(2, "mm"),
                         row_gap = unit(4, "mm"),
                         rect_gp = gpar(col = "black", lwd = 2),
                         row_names_gp = gpar(fontsize = 8), 
                         row_names_side = "right", 
                         row_title_side = "left",
                         cluster_rows = TRUE,
                         cluster_columns = TRUE,
                         cluster_column_slices = FALSE,
                         cluster_row_slices = FALSE,
                         clustering_method_rows = "complete",
                         left_annotation = significance_annoatation)

lgd_pvalue = Legend(title = "p-value", col_fun = pvalue_col_fun, at = c(0, 1, 2, 3), 
                    labels = c("1", "0.1", "0.01", "0.001"))
lgd_sig = Legend(pch = "*", type = "points", labels = "< 0.05")

#pdf("day11_polymer_heatmap_clustered_significance.pdf", width = 10, height = 10)
#ht_both_clust
#draw(ht_both_clust, annotation_legend_list = list(lgd_pvalue, lgd_sig))
#dev.off()


