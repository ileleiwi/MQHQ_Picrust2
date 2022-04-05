library(tidyverse)
library(readxl)
library(RColorBrewer)


setwd(paste0("/Users/ikaialeleiwi/Desktop/Lab/Salmonella_NIH/Lactobacillus/",
             "Omics/Metagenome/rRNA/16S_From_All_Bins/MQHQ_Picrust2/",
             "MQHQ_Picrust2/"))

#Data 
#DRAM genome summary form
genome_form_cazy <- read_tsv("Data/genome_summary_form.tsv") %>%
  filter(header == "CAZY") 

#DRAM function heatmap form
function_heatmap_form_cazy <- read_tsv("Data/function_heatmap_form.tsv") %>%
  filter(category == "CAZy") %>%
  select(-gene_symbol)

#function rules
rules <- read_tsv("Clean_Data/picr_rules.tsv")

#picrust ASV gene counts and DRAM bin gene counts combined df
picr_dram <- read_tsv("Clean_Data/picr_dram.tsv") 

picr_dram_present_absent <- picr_dram %>%
  mutate(picrust_pa = ifelse(counts_picrust > 0,
                              yes = 1,
                              no = 0),
         dram_pa = ifelse(counts_dram > 0,
                          yes = 1,
                          no = 0))

expected_value <- factor(picr_dram_present_absent$dram_pa, 
                         levels = c(1, 0))
predicted_value <- factor(picr_dram_present_absent$picrust_pa, 
                          levels = c(1, 0))

#ASV to Bins match statistics
asv_bins <- read_tsv("Data/match_statistics.tsv") %>%
  mutate(bin.id = sub(x = bin_scaffold_header, 
                      pattern = "_scaffold.*|_k.*",
                      replacement = "")) %>% # <---- create column with bin name only
  filter(bin.id %in% picr_dram_present_absent$bin.id)

#picrust2 predicted ECs
picr_ecs <- read_tsv("Data/EC_predicted.tsv") %>%
  filter(sequence %in% asv_bins$ASV_header)

#DRAM metabolism summaray (for EC's)
#dram gene counts
metab_summary <- "Data/metabolism_summary_deduped.xlsx" %>%
  excel_sheets() %>%
  set_names() %>%
  map_dfr(read_excel, path = "Data/metabolism_summary_deduped.xlsx", 
          .id = "SheetName") 


#functions
#match KO number format K----- 
matchKO <- function(string){
  x <- str_extract_all(string, 
                       "K[\\d]{5}")
  return(unlist(x))
}

#match EC number format EC-.-.-.- or EC:-.-.-.-
matchEC <- function(string, format="", tf=F){
  if(format == ":"){
    x <- str_extract_all(string, 
                         "EC:[\\d-]{1,2}\\.[\\d-]{1,2}\\.[\\d-]{1,2}\\.[\\d-]{1,2}")
  }else{
    x <- str_extract_all(string, 
                         "EC[\\d-]{1,2}\\.[\\d-]{1,2}\\.[\\d-]{1,2}\\.[\\d-]{1,2}")
  }
  
  if(tf & format == ":"){
    x <- str_detect(string, 
                    "EC:[\\d-]{1,2}\\.[\\d-]{1,2}\\.[\\d-]{1,2}\\.[\\d-]{1,2}")
  }else if(tf){
    x <- str_detect(string, 
                    "EC[\\d-]{1,2}\\.[\\d-]{1,2}\\.[\\d-]{1,2}\\.[\\d-]{1,2}")
  }
  return(unlist(x))
}



#pull all KO's and EC's from rules dataframes
KO_keep <- c()
EC_keep <- c()
for(row in 1:nrow(function_heatmap_form_cazy)){
  KO_keep <- c(KO_keep, matchKO(function_heatmap_form_cazy[row,"function_ids"]))
}

for(row in 1:nrow(genome_form_cazy)){
  KO_keep <- c(KO_keep, matchKO(genome_form_cazy[row,"gene_description"]))
}

for(row in 1:nrow(rules)){
  KO_keep <- c(KO_keep, matchKO(rules[row,"child"]))
}

for(row in 1:nrow(rules)){
  EC_keep <- c(EC_keep, matchEC(rules[row,"child"]))
}

for(row in 1:nrow(rules)){
  print(matchEC(EC_keep, format = ":"))
}


KO_keep <- unique(KO_keep)
EC_keep <- unique(EC_keep) %>%
  str_replace(pattern = "EC",
              replacement = "EC:")

#filter df to only include KO's relevant to functions and dbcan ids
picr_dram_present_absent <- picr_dram_present_absent %>%
  mutate(gene_id = ifelse(!startsWith(gene_id, "K") | gene_id %in% KO_keep,
                          yes = gene_id,
                          no = NA)) %>%
  drop_na(gene_id)

#filter metabolism summary to include only rows with EC's and relevent bin columns
metab_ecs <- c()
for(row in 1:nrow(metab_summary)){
  ecs <- matchEC(string = as.character(metab_summary[row,"gene_description"]), format = ":")
  if(length(ecs) > 0){
    metab_ecs <- c(metab_ecs, ecs)
  }else{
    metab_ecs <- c(metab_ecs, "ecs")
  }
}

row_keep <- metab_ecs %in% EC_keep %>% which()

metab_filtered <- metab_summary[row_keep,] %>% 
  as.data.frame() %>%
  mutate(gene_id = matchEC(gene_description, ":")) %>%
  select(gene_id, all_of(unique(picr_dram_present_absent$bin.id))) %>%
  distinct() %>%
  pivot_longer(cols = -gene_id,
               names_to = "bin.id",
               values_to = "counts_dram") %>%
  group_by(bin.id, gene_id) %>%
  summarise(counts_dram = sum(counts_dram))

#filter picr_ecs to include only relevant ecs
picr_ecs_filtered <- picr_ecs %>%
  select(sequence, any_of(EC_keep))

#all 0 so just adding a counts_picrust column to metab_filtered with all 0's
metab_filtered$counts_picrust <-  rep(0, nrow(metab_filtered))
metab_filtered$picrust_pa <- rep(0, nrow(metab_filtered))
metab_filtered$dram_pa <- rep(0, nrow(metab_filtered))

metab_filtered <- metab_filtered %>%
  select(colnames(picr_dram_present_absent))

#add metab_filteredto picr_dram_present_absent
picr_dram_present_absent <- rbind(picr_dram_present_absent, metab_filtered)


#make new df with functions as rows and  columns for ids and # of ids
funct_df <- data.frame(functions = unique(function_heatmap_form_cazy$function_name),
                       ids = rep(NA, length(unique(function_heatmap_form_cazy$function_name))),
                       num_ids = rep(NA, length(unique(function_heatmap_form_cazy$function_name))))

for(row in 1:nrow(function_heatmap_form_cazy)){
  for(funct in 1:nrow(funct_df)){
    if(function_heatmap_form_cazy[row,"function_name"] == funct_df[funct, "functions"]){
      funct_df[funct, "ids"] <- paste(funct_df[funct, "ids"], function_heatmap_form_cazy[row,"function_ids"], sep = ", ")
    }
  }
}

funct_df <- funct_df %>%
  mutate(ids = str_remove(ids, "NA, "))

for(row in 1:nrow(funct_df)){
  x <- unlist(strsplit(funct_df[row, "ids"], ", "))
  funct_df[row, "ids"] <- paste(unique(x), collapse = ", ")
  funct_df[row, "num_ids"] <- length(x)
}


#collapse rule dataframe
ETC50 <- rules %>%
  filter(name == "ETC50") %>%
  pull(child) %>%
  matchKO() %>%
  unique()

TMAO <- rules %>%
  filter( name == "TMAO") %>%
  pull(child) %>%
  matchKO() %>%
  unique()

DenitrificationNotETCKO <- rules %>%
  filter( name == "DenitrificationNotETC") %>%
  pull(child) %>%
  matchKO() %>%
  unique()

DenitrificationNotETCEC <- rules %>%
  filter( name == "DenitrificationNotETC") %>%
  pull(child) %>%
  matchEC() %>%
  unique()

DenitrificationNotETC <- paste(unique(c(DenitrificationNotETCKO, DenitrificationNotETCEC)), collapse = ", ")

MicroaerophillicKO <- rules %>%
  filter( name == "Microaerophillic") %>%
  pull(child) %>%
  matchKO() %>%
  unique()

Microaerophillic <- paste(unique(c(MicroaerophillicKO, ETC50)), collapse = ", ")

TetrathionateReduction <- rules %>%
  filter( name == "TetrathionateReduction") %>%
  pull(child) %>%
  matchKO() %>%
  unique() %>%
  paste(collapse = ', ')

AerobicKO <- rules %>%
  filter( name == "Aerobic") %>%
  pull(child) %>%
  matchKO() %>%
  unique()

Aerobic <- paste(unique(c(AerobicKO, ETC50)), collapse = ", ")


FumarateReduction <- rules %>%
  filter( name == "FumarateReduction") %>%
  pull(child) %>%
  matchKO() %>%
  unique() %>%
  paste(collapse = ", ")

SulfateReductionKO <- rules %>%
  filter( name == "SulfateReduction") %>%
  pull(child) %>%
  matchKO() %>%
  unique()

SulfateReductionEC <- rules %>%
  filter( name == "SulfateReduction") %>%
  pull(child) %>%
  matchEC() %>%
  unique()

SulfateReduction <- paste(unique(c(SulfateReductionKO, SulfateReductionEC)), collapse = ", ")


rules_df <- data.frame(functions = unique(rules$name),
                       ids = rep(NA, length(unique(rules$name))),
                       num_ids = rep(NA, length(unique(rules$name)))) %>%
  filter(functions != "ETC50")

for(row in 15:nrow(rules_df)){
  rules_df[row,"ids"] <- get(rules_df[row,"functions"])
}

for(row in 1:14){
  funct <- rules_df[row, "functions"]
  x <- matchKO(rules[rules$name == funct, "child"])
  
  rules_df[row, "ids"] <- paste(unique(x), collapse = ", ")
  rules_df[row, "num_ids"] <- length(x)
}

for(row in 15:nrow(rules_df)){
  string <- unlist(strsplit(rules_df[row,"ids"], ", "))
  rules_df[row,"num_ids"] <- length(string)
}

#dataframe with all functions and gene_ids for each function
funct_df <- rbind(funct_df, rules_df)

#calculate the % present of ids in picrust_pa and dram_pa
#need to account for each bin and each gene
#sum all gene calls from each bin

picr_dram_sum_bins <- picr_dram %>%
  group_by(gene_id) %>%
  summarise(sum_counts_picrust = sum(counts_picrust),
            sum_counts_dram = sum(counts_dram))

metab_filtered_sum_bins <- metab_filtered %>%
  group_by(gene_id) %>%
  summarise(sum_counts_picrust = sum(counts_picrust),
            sum_counts_dram = sum(counts_dram))

summed_df <- rbind(picr_dram_sum_bins, metab_filtered_sum_bins)

summed_df_pa <- summed_df %>%
  mutate(picrust_pa = ifelse(sum_counts_picrust == 0,
                             yes = 0,
                             no = 1),
         dram_pa = ifelse(sum_counts_dram == 0,
                             yes = 0,
                             no = 1),
         gene_id = case_when(startsWith(gene_id, "EC") ~ str_replace(gene_id,
                                                                     "EC:",
                                                                     "EC"),
                             TRUE ~ gene_id))


#query funct_df and summed_df_pa for gene id's presence and absence
funct_df <- funct_df %>%
  mutate(picrust_pa = NA,
         pi_prop = NA,
         dram_pa = NA,
         dr_prop = NA)

for(row in 1:nrow(funct_df)){
  ids <- funct_df[row, "ids"]
  id_vect <- unlist(strsplit(ids, ", "))
  
  id_count_pi <- 0
  id_count_dr <- 0
  for(id in id_vect){
    if(id %in% summed_df_pa$gene_id){
      pi_pa <- summed_df_pa %>%
        filter(gene_id == id) %>%
        pull(picrust_pa)
      
      dr_pa <- summed_df_pa %>%
        filter(gene_id == id) %>%
        pull(dram_pa)
      
      if(pi_pa == 1){
        id_count_pi <- id_count_pi + 1
      }
      
      if(dr_pa == 1){
        id_count_dr <- id_count_dr + 1
      }
    }else{
      next
    }
  }
  
  funct_df[row, "picrust_pa"] <- id_count_pi
  funct_df[row, "dram_pa"] <- id_count_dr
  funct_df[row, "pi_prop"] <- (id_count_pi/funct_df[row, "num_ids"])*100
  funct_df[row, "dr_prop"] <- (id_count_dr/funct_df[row, "num_ids"])*100
}


write_tsv(funct_df, "Clean_Data/function_pa_dram_picrust.tsv")
#confustion matrix
m <- table(expected_value, predicted_value)
TP <- m[1,1] #true positives
TN <- m[2,2] #true negatives
FP <- m[2,1] #false positives
FN <- m[1,2] #false negatives
P <- TP + FN #real positives
N <- FP + TN #real neagatives


TPR <- TP/P #Sensitivity (1 - FDR) (TP/TP+FN) (True positive rate) #how good at detecting positive events
FDR <- 1-TPR #False discovery rate (1 - TPR) 
NPV <- TN/N #Specificity (1 - FOR) (TN/FP+TN) #how how exact the assignment to the positive class is
FOR <- 1-NPV #False omission rate (1 - NPV)



RECALL <- TP/P #measures how good the model is in detecting positive events
PERCISION <- TP/(TP+FP) # how good the model is at assigning positive events
ACCURACY <- (TP + TN)/(TP + TN + FP + FN)
FMEASURE <- (2*RECALL*PERCISION)/(RECALL+PERCISION)

#match output to metabolisms
#picharts

plot_df <- funct_df %>%
  mutate(p_not = 100 - pi_prop,
         d_not = 100 - dr_prop) %>%
  pivot_longer(cols = c("pi_prop", "dr_prop", "p_not", "d_not"),
               values_to = "prop",
               names_to = "method") %>%
  mutate(method = case_when(str_detect(method, "p_not") ~ "Picrust2 Not Predicted",
                            str_detect(method, "d_not") ~ "DRAM Not Predicted",
                            str_detect(method, "pi") ~ "Picrust2",
                            str_detect(method, "dr") ~ "DRAM"),
         label = paste0(as.character(round(prop, 1)), "/", as.character(num_ids)))

pal <- brewer.pal(n=11,name = "RdBu")[c(1,5,11,7)]

plot_out <- plot_df %>%
  ggplot(aes(x = "", y = prop, fill = method)) +
  geom_bar(position = "fill", stat = "identity", color = "black") +
  coord_polar("y", start = 0) +
  scale_fill_manual(values = pal) +
  facet_wrap(~functions, ncol = 8) +
  theme_void() +
  theme(legend.position = "bottom") 


# pdf("Figs/function_piecharts.pdf", width = 10, height = 10)
# plot_out
# dev.off()

#heatmap for confusion matrix
hm_out <- m %>%
  as.data.frame() %>%
  ggplot(aes(x = expected_value, y = predicted_value, fill = Freq)) +
  geom_tile(color = "black", size = 3) +
  scale_fill_distiller(palette="Greens", direction=1) +
  geom_text(aes(label = Freq), color = "black") +
  theme_void() +
  theme(legend.position = "none")

# pdf("Figs/confustion_matrix.pdf", width = 10, height = 10)
# hm_out
# dev.off()
