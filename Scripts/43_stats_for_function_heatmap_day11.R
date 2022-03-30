## ---------------------------
##
## Script: 43_stats_for_function_heatmap_day11
##
## Purpose: run statistics on different functions and groups
##
## Author: Ikaia Leleiwi
##
## Date Created: 2022-03-29
##
## Copyright (c) Ikaia Leleiwi, 2022
## Email: ileleiwi@gmail.com
##
## ---------------------------
##
## Notes: to be run in 53_function_heatmap_day11_picr.R
##   
##
## ---------------------------

## set working directory

setwd(paste0("/Users/ikaialeleiwi/Desktop/Lab/Salmonella_NIH/Lactobacillus/",
             "Omics/Metagenome/rRNA/16S_From_All_Bins/MQHQ_Picrust2/",
             "MQHQ_Picrust2/"))

## ---------------------------

library(tidyverse)

PATH <- paste0("Clean_Data/samples_by_function_out_day11.tsv")

df <- read_tsv(PATH) %>%
  rename("KE" = "13-(11)-S-F-KE",
         "KL" = "20-(11)-S-F-KL",
         "KS" = "27-(11)-S-F-KS",
         "LO" = "49-(11)-C-F-LO",
         "LM" = "47-(11)-C-F-LM",
         "LP" = "50-(11)-C-F-LP") %>%
  mutate(category = c(rep(1, 19), rep(2, 9), rep(3, 5), rep(4, 7)),
         category = case_when(category == 1 ~ "Polymers",
                              category == 2 ~ "Sugars",
                              category == 3 ~ "SCFA",
                              category == 4 ~ "Respiration",
                              TRUE ~ "other"),
         category = as.factor(category),
         `function` = as.factor(`function`)) %>%
  pivot_longer(cols = c("LM", "LO", "LP", "KE", "KL", "KS"),
               names_to = "sample",
               values_to = "relative_abundance") %>%
  mutate(treatment = case_when(sample %in% c("LM", "LO", "LP") ~ "Control",
                               sample %in% c("KE", "KL", "KS") ~ "Infected"),
         sample = as.factor(sample),
         treatment = as.factor(treatment)) %>%
  select(treatment, sample, category, `function`, relative_abundance) %>%
  rename("functions" = "function")

#anova on each function between control and treatment
anova_list <- list()
for( i in unique(df$functions)){
  x <- df[df$functions == i,]
  fit <- lm(relative_abundance ~ treatment, data = x)
  anova <- anova(fit)
  anova_list[[length(anova_list) + 1]] <- anova
}
names(anova_list) <- unique(df$functions)

#anova on each category between control and treatment
anova_list_cat <- list()
for( i in unique(df$category)){
  x <- df[df$category == i,]
  fit <- lm(relative_abundance ~ treatment, data = x)
  anova <- anova(fit)
  anova_list_cat[[length(anova_list_cat) + 1]] <- anova
}
names(anova_list_cat) <- unique(df$category)

#dataframe of significance between treatments for each function
significant <- c()
pval <- c()
for (i in names(anova_list)){
  x <- as.data.frame(anova_list[i])
  if(is.nan(x["treatment",5])){
    significant <- c(significant, NaN)
    pval <- c(pval, NaN)
  }
  else if(x["treatment",5] <=0.05){
    significant <- c(significant, as.character(x["treatment",5]))
    pval <- c(pval, x["treatment",5])
  }
  else{
    significant <- c(significant, "ns")
    pval <- c(pval, x["treatment",5])
  }
}

#dataframe of significance between treatments for each function
significant_cat <- c()
pval_cat <- c()
for (i in names(anova_list_cat)){
  x <- as.data.frame(anova_list_cat[i])
  if(is.nan(x["treatment",5])){
    significant_cat <- c(significant_cat, NaN)
    pval_cat <- c(pval_cat, NaN)
  }
  else if(x["treatment",5] <=0.05){
    significant_cat <- c(significant_cat, as.character(x["treatment",5]))
    pval_cat <- c(pval_cat, x["treatment",5])
  }
  else{
    significant_cat <- c(significant_cat, "ns")
    pval_cat <- c(pval_cat, x["treatment",5])
  }
}


significance_df <- data.frame(functions = names(anova_list),
                              significance = significant,
                              pval = pval)

significance_df_cat <- data.frame(category = names(anova_list_cat),
                                  significance = significant_cat,
                                  pval = pval_cat)


rm(anova, anova_list, df, fit, x, PATH, i, significant, pval,
   anova_list_cat, significant_cat, pval_cat)

