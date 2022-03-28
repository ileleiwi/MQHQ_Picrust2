## ---------------------------
##
## Script: 33_function_heatmap_day11_picr
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

picr_func <- read_tsv("Clean_Data/picr_day11_functions.tsv") %>%
  left_join(taxonomy, by = c("genome" = "FeatureID")) %>%
  filter(!str_detect(Taxon, "D_4__Mitochondria|D_3__Chloroplast")) %>%# remove mitochondria and chloroplast from picrust2
  select(-Taxon, -Confidence) %>%
  column_to_rownames(var = "genome")
  
  
ft <- read_tsv("Clean_Data/feature_table_clean_r1-r5.tsv") %>%
  select(asv_id, all_of(metaG_samples)) 

dram_product <- read_tsv("Clean_Data/dram_product_picrust_day11.tsv")

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


CountAcrossRows <- function(working_list){
  #count number of genomes within each sample that are more than 0.5% relative abundance
  #for each function
  for(i in names(working_list)){
    sample_matrix <- ifelse(working_list[[i]] >= 0.5, 1, 0)
    sample_matrix <- apply(sample_matrix, 1, sum)
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

#create list of number of genomes more than 0.5% relative abundancefor each sample and function
genome_count <- CountAcrossRows(generateRelabundFunctList(function_df = picr_func, feature_table = ft_rel))

genome_count_df <- do.call(cbind.data.frame, genome_count)



#scale by function
samples_by_function_scaled <- t(scale(t(samples_by_function)))
samples_by_function_scaled[is.nan(samples_by_function_scaled)] <- 0
samples_by_function_scaled <- as.data.frame(samples_by_function_scaled)

#modify rownames
samples_by_function_scaled <- samples_by_function_scaled %>%
  rownames_to_column(var = "functions") %>%
  mutate(case_when())

#Split dataframes by healthy and infected for plotting
healthy <- samples_by_function_scaled %>%
  select(LM, LO, LP) %>%
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

infected <- samples_by_function_scaled %>%
  select(KE, KL, KS) %>%
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

healthy_unscaled <- samples_by_function %>%
  select(LM, LO, LP) %>%
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
  select(KE, KL, KS) %>%
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

#dataframes for heatmap barplot annotations
genomes_healthy <- genome_count_df %>%
  select(LM, LO, LP) %>%
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
  select(KE, KL, KS) %>%
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
  filter(functions %in% colnames(cazy)) %>%
  arrange(match(functions, c("Polyphenolics", "Crystalline Cellulose","Amorphous Cellulose", 
                             "Mixed-Linkage Glucans","Xylans","Xyloglucan", "Alpha-Mannan", 
                             "Beta-Mannan", "Mucin", "Starch", "Pectin", "Chitin", "Arabinan",
                             "Alpha-Galactans", "Beta-Galactan (Pectic Galactan)",
                             "Sulf-Polysachharides", "Arabinose Cleavage","Rhamnose Cleavage",
                             "Fucose Cleavage"))) %>%
  column_to_rownames(var = "functions")

sug_sbf <- samples_by_function_scaled %>%
  rownames_to_column(var = "functions") %>%
  filter(functions %in% colnames(sugar)) %>%
  arrange(match(functions, c("Lactose", "Sucrose", "Galactose", "Mannose", 
                             "Galacturonic Acid", "Hexose", "Xylose", "Fucose", 
                             "Fructose"))) %>%
  column_to_rownames(var = "functions")

scfa_sbf <- samples_by_function_scaled %>%
  rownames_to_column(var = "functions") %>%
  filter(functions %in% colnames(scfa)) %>%
  arrange(match(functions, c("Butyrate", "Propionate", "Acetate", "Lactate", 
                             "Alcohol"))) %>%
  column_to_rownames(var = "functions")

resp_sbf <- samples_by_function_scaled %>%
  rownames_to_column(var = "functions") %>%
  filter(functions %in% colnames(respiration)) %>%
  arrange(match(functions, c("Sulfate Reduction", "Fumarate Reduction",
                             "Aerobic", "Tetrathionate Reduction", "Microaerophillic",
                             "Denitrification (Not ETC)", "TMAO"))) %>%
  column_to_rownames(var = "functions")
