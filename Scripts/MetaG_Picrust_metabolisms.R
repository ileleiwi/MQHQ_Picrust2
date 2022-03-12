library(tidyverse)


setwd(paste0("/Users/ikaialeleiwi/Desktop/Lab/Salmonella_NIH/Lactobacillus/",
             "Omics/Metagenome/rRNA/16S_From_All_Bins/MQHQ_Picrust2/",
             "MQHQ_Picrust2/"))


#Data 
picr_ko <- read_tsv("Clean_Data/picr_dram.tsv") %>%
  filter(str_detect(gene_id, "K\\d{5}")) %>%
  select(-counts_dram) %>%
  rename("counts" = "counts_picrust")

picr_ec <- read_tsv("Clean_Data/picrust_ec_mqhq_match.tsv") 

picr <- rbind(picr_ko, picr_ec)

#rules
source("Scripts/change_Rules_to_list.R")

#make wide dataframes
picr <- picr %>%
  pivot_wider(values_from = "counts",
              names_from = "gene_id") %>%
  column_to_rownames(var = "bin.id") %>%
  mutate(across(everything(),
                ~ifelse(.x > 0, 1, .x)))


#functions to parse through rules
matchEC <- function(string){
    #match EC number format EC-.-.-.-
    x <- str_extract_all(string, 
                         "EC[\\d-]{1,2}\\.[\\d-]{1,2}\\.[\\d-]{1,2}\\.[\\d-]{1,2}")
  return(x)
}

matchKO <- function(string){
  #match KO number format K-----
  x <- str_extract_all(string, 
                       "K[\\d]{5}")
  return(x)
}

parseOr <- function(string){
 if(str_detect(string, "|")){
   return(TRUE)
 }else{
   return(FALSE)
 }
}

parseAnd <- function(string){
  if(str_detect(string, "&")){
    return(TRUE)
  }else{
    return(FALSE)
  }
}

parse50 <- function(string){
  if(str_detect(string, "percent50")){
    return(TRUE)
  }else{
    return(FALSE)
  }
}



checkMainRule <- function(rules = main_rules, funct){
  return(rules[[funct]])
}

checkSubRule <- function(rules = sub_rules, funct){
  return(rules[[funct]])
}






priorityString <- function(string){
  strlist <- str_match_all(string, "\\[.+?\\]")
  #strlist_trimmed <- str_sub(strlist, 2, -2)
  if(length(strlist) == 0){
    return(FALSE)
  }else{
    return(TRUE) 
  }
}


compIA <- "K00330,[K00331&K00332&[K00333|K00331]&[K13378|K13380]],K00334,K00335,K00336,K00337,K00338,K00339,K00340,[K00341&K00342|K15863],K00343"
test_str <- "[percent50CytochromeCOxidase|percent50CytochromeAa3600MenaquinolOxidase|percent50CytochromeOUbiquinolOxidase]&ETC50"
compIB <- "K05574,K05582,K05581,K05579,K05572,K05580,K05578,K05576,K05577,K05575,K05573,K05583,K05584,K05585"

test_p_list <- buildPriorityList(buildIndex(compIA), compIA)

buildIndex <- function(rule_string){
  index <- c()
  for(i in 1:nchar(rule_string)){
    if(str_sub(rule_string, i, i) == "["){
      index <- c(index, i)
    }else if(str_sub(rule_string, i, i) == "]"){
      index <- c(index, -i)
    }
  }
  return(index)
}


buildPriorityList <- function(index, rule_string){
  priority_list <- list()
  dbl_idx_strt <- c()
  dbl_idx_end <- c()
  
  if(is.null(index)){
    priority_list[length(priority_list)+1] <- rule_string
    names(priority_list) <- "org"
    return(priority_list)
  }else{
    for(i in 2:length(index)){
      if(index[i] > 0 & index[i-1] > 0){
        dbl_idx_strt <- c(dbl_idx_strt, index[i-1])
      }else if(index[i] < 0 & index[i-1] > 0){
        x <- unlist(str_sub(rule_string, index[i-1], -index[i]))
        priority_list[length(priority_list)+1] <- str_remove_all(x, "\\[|\\]")
      }else if(index[i] < 0 & index[i-1] < 0){
        dbl_idx_end <- c(dbl_idx_end, -index[i-1])
      }
    }
    names(priority_list) <- letters[1:length(priority_list)]
    first_names <- names(priority_list)
    
    dbl_names <- c()
    if(length(dbl_idx_strt) != 0 & length(dbl_idx_end) != 0){
      for(i in 1:length(dbl_idx_strt)){
        priority_list[length(priority_list)+1] <- str_sub(compIA, dbl_idx_strt[i], dbl_idx_end[i])
        dbl_names <- c(dbl_names, paste0("dbl", as.character(i)))
      }
    }
    
    priority_list[length(priority_list)+1] <- rule_string
    
    names(priority_list) <- c(first_names, dbl_names, "org")
    return(priority_list)
  }
}




evaluateRule <- function(string){
  if(priorityString(string)){
    FIRST <- str_sub(str_match_all(string, "\\[.+?\\]"), 2, -2)
    FIRST_ANDs <- 
  }
  if(parseAnd(string)){
    ANDs <- str_split(string, "&")
  }
  if(parseOr(string)){
    ORs <- str_split(string, "|")
  }
}



checkSulfRed <- function(dataframe){
  main <- checkMainRule("Sulfate Reduction")
  sub <- checkSubRule(main)
  parseOr(sub)
}

####EC50 Functions####

compIA50 <- function(dataframe, row){
  df <- as.data.frame(dataframe) %>%
    slice(row) %>%
    select(all_of(c("K00330", "K00331", "K00332", "K00333", "K00331", "K13378",
                    "K13380", "K00334", "K00335", "K00336", "K00337", "K00338", 
                    "K00339", "K00340", "K00341", "K00342", "K15863", 
                    "K00343")))
  #first condition
  if(df[1, "K13378"] == 1 | df[1, "K13380"] == 1){
    A <- 1
  }else{
    A <- 0
  }
  
  #second condition
  if(df[1, "K00331"] == 1 & df[1, "K00332"] == 1 & df[1, "K00333"] == 1){
    B <- 1
  }else{
    B <- 0
  }
  
  ##third condition
  if(A & B){
    C <- 1
  }else{
    C <- 0
  }
  
  #fourth condition
  if(df[1, "K00342"] == 1 | df[1, "K15863"] == 1){
    D <- TRUE
  }else{
    D <- 0
  }
  
  ##fifth condition
  if(df[1, "K00341"] == 1 & D){
    E <- 1
  }else{
    E <- 0
  }

  if(sum(df[1,c("K00330","K00334","K00335","K00336","K00337","K00338","K00339","K00340","K00343")], C, E) >= 6){
    return(TRUE)
  }else{
    return(FALSE)
  }
}

compIB50 <- function(dataframe, row){
  df <- as.data.frame(dataframe) %>%
    slice(row) %>%
    select(all_of(c("K05574","K05582","K05581","K05579","K05572","K05580",
                    "K05578","K05576","K05577","K05575",'K05573',"K05583",
                    "K05584","K05585")))
  
  if(sum(df[1,]) >= 7){
    return(TRUE)
  }else{
    return(FALSE)
  }
}
  
compIC50 <- function(dataframe, row){
  df <- as.data.frame(dataframe) %>%
    slice(row) %>%
    select(all_of(c("K03945","K03946","K03947","K03948","K03949","K03950",
                    "K03951","K03952","K03953","K03954","K03955","K03956",
                    "K11352","K11353")))
  
  if(sum(df[1,]) >= 7){
    return(TRUE)
  }else{
    return(FALSE)
  }
}

checkEC50 <- function(dataframe, row){
  if(compIA50(dataframe, row) | 
     compIB50(dataframe, row) | 
     compIC50(dataframe, row)){
       return(TRUE)
     }else{
       return(FALSE)
     }
}







#Rules

# Sulfate Reduction 
#   must have Dsr
#   Dsr - K11180|EC1.8.99.5|K11181|EC1.8.99.5
#   Asr - K00380|EC1.8.1.2|K00381|EC1.8.1.2|K00392|EC1.8.7.1
# 
# Fumarate Reduction
#   must have one of the following KOs
#   K00244,K00245,K00246,K00247
# 
# Aerobic
#   must have [50% CytochromeCOxidase OR 50% CytochromeAa3600MenaquinolOxidase OR 50% CytochromeOUbiquinolOxidase] AND ETC50
# 
#   ETC50
#     must have 50% ComplexIA OR 50% ComplexIB OR 50% ComplexIC
#     ComplexIA - K00330,[K00331&K00332&[K00333|K00331]&[K13378|K13380]],K00334,K00335,K00336,K00337,K00338,K00339,K00340,[K00341&K00342|K15863],K00343
#     ComplexIB - K05574,K05582,K05581,K05579,K05572,K05580,K05578,K05576,K05577,K05575,K05573,K05583,K05584,K05585
#     ComplexIC - K03945,K03946,K03947,K03948,K03949,K03950,K03951,K03952,K03953,K03954,K03955,K03956,K11352,K11353
#   
#   CytochromeCOxidase - K02275,K02274,K02276|K15408,K02277
#   CytochromeAa3600MenaquinolOxidase - K02827,K02826,K02828,K02829
#   CytochromeOUbiquinolOxidase - K02297,K02298,K02299,K02300
# 
# 
# Tetrathionate Reduction
#   must have [1 of ttrrs AND 1 of ttrabc] OR 2 of ttrabc
#   ttrabc - K08359,K08358,K08357
#   ttrrs - K13040|K13041
# 
# Microaerophillic
#   must have [50% CytochromeBDUbiquinolOxidase OR 50% CytochromeCOxidaseCBB3] AND ETC50
#     
#     CytochromeBDUbiquinolOxidase - K00425,K00426,K00424|K22501
#     CytochromeCOxidaseCBB3 - K00404|K00405|K15862,K00407,K00406
#     
# Denitrification (Not ETC)
#   must have [NitrateReductase OR NitriteReductase OR NitricOxideReductase OR NitrousOxideReductase] AND does not have ETC50
#     NitrateReductase - K00370|K00371|K00374|K02567|K02568|EC1.7.5.1|EC1.9.6.1
#     NitriteReductase - K00368|K15864|EC1.7.2.1
#     NitricOxideReductase - K04561|K02305|EC1.7.2.5
#     NitrousOxideReductase - K00376|EC1.7.2.4