library(tidyverse)


setwd(paste0("/Users/ikaialeleiwi/Desktop/Lab/Salmonella_NIH/Lactobacillus/",
             "Omics/Metagenome/rRNA/16S_From_All_Bins/MQHQ_Picrust2/",
             "MQHQ_Picrust2/"))


#Data 
picr_ko <- read_tsv("Clean_Data/picr_dram.tsv") %>%
  filter(str_detect(gene_id, "K\\d{5}")) %>%
  select(-counts_dram) %>%
  rename("counts" = "counts_picrust")

picr_ec <- read_tsv("Clean_Data/picrust_ec_mqhq_match.tsv") %>%
  mutate(gene_id = str_remove(gene_id, ":"))

picr <- rbind(picr_ko, picr_ec) %>%
  as.data.frame()

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
  return(unlist(x))
}

matchEC("K00376|EC1.7.2.4")
str_extract_all("K00376|EC1.7.2.4", 
                "EC[\\d-]{1,2}\\.[\\d-]{1,2}\\.[\\d-]{1,2}\\.[\\d-]{1,2}")

matchKO <- function(string){
  #match KO number format K-----
  x <- str_extract_all(string, 
                       "K[\\d]{5}")
  return(unlist(x))
}

#match if 50% marker is present
check50 <- function(string){
  if(str_detect(string, "percent50")){
    loc_50 <- unlist(str_match_all(string, "percent50(.+?)[\\|,\\]]|percent50(.+?)$"))
    loc_50_1 <- loc_50[!str_detect(loc_50, "percent50")]
    return(loc_50_1[!is.na(loc_50_1)])
  }else{
    return(NA)
  }
}


#match if 1of  marker is present
one <- function(string){
  if(str_detect(string, "1of")){
    loc_1 <- unlist(str_match_all(string, "1of(.+?)[\\|,\\]].*$"))
    return(loc_1[!str_detect(loc_1, "1of")])
  }else{
    return(NA)
  }
}


#match if 2of  marker is present
two <- function(string){
  if(str_detect(string, "2of")){
    loc_2 <- unlist(str_match_all(string, "^.*2of(.+?)$"))
    return(loc_2[!str_detect(loc_2, "2of")])
  }else{
    return(NA)
  }
}

checkMainRule <- function(funct){
  switch(funct,
         SulfateReduction = "Dsr",
         FumarateReduction = "FumarateReduction",
         Aerobic = "[percent50CytochromeCOxidase|percent50CytochromeAa3600MenaquinolOxidase|percent50CytochromeOUbiquinolOxidase]&ETC50",
         TetrathionateReduction = "[ttrrs&1ofttrabc]|2ofttrabc",
         Microaerophillic = "[percent50CytochromeBDUbiquinolOxidase|percent50CytochromeCOxidaseCBB3]&ETC50",
         DenitrificationNotETC = "notETC50&[NitrateReductase|NitriteReductase|NitricOxideReductase|NitrousOxideReductase]", 
         ETC50 = "percent50ComplexIA|percent50ComplexIB|percent50ComplexIC",
         stop("Unknown input"))
}

pullIds("[ttrrs&1ofttrabc]|2ofttrabc")

checkSubRule <- function(funct){
  switch(funct,
         Dsr = "K11180|EC1.8.99.5|K11181|EC1.8.99.5",
         FumarateReduction = "K00244|K00245|K00246|K00247",
         CytochromeCOxidase = "K02275,K02274,K02276|K15408,K02277",
         CytochromeAa3600MenaquinolOxidase = "K02827,K02826,K02828,K02829",
         CytochromeOUbiquinolOxidase = "K02297,K02298,K02299,K02300",
         ETC50 = "percent50ComplexIA|percent50ComplexIB|percent50ComplexIC",
         ttrrs = "K13040|K13041",
         ttrabc = "K08359,K08358,K08357",
         CytochromeBDUbiquinolOxidase = "K00425,K00426,K00424|K22501",
         CytochromeCOxidaseCBB3 = "K00404|K00405|K15862,K00407,K00406",
         NitrateReductase = "K00370|K00371|K00374|K02567|K02568|EC1.7.5.1|EC1.9.6.1",
         NitriteReductase = "K00368|K15864|EC1.7.2.1",
         NitricOxideReductase = "K04561|K02305|EC1.7.2.5",
         NitrousOxideReductase = "K00376|EC1.7.2.4",
         ComplexIA = "K00330,[K00331&K00332&[K00333|K00331]&[K13378|K13380]],K00334,K00335,K00336,K00337,K00338,K00339,K00340,[K00341&K00342|K15863],K00343",
         ComplexIB = "K05574,K05582,K05581,K05579,K05572,K05580,K05578,K05576,K05577,K05575,K05573,K05583,K05584,K05585",
         ComplexIC = "K03945,K03946,K03947,K03948,K03949,K03950,K03951,K03952,K03953,K03954,K03955,K03956,K11352,K11353",
         stop("Unknown input"))
}

setETC50 <- function(dataframe, genome){
  df <- dataframe %>%
    filter(bin.id == genome)
  
  main_ETC50 <- checkMainRule("ETC50")
  main_ETC50 <- check50(main_ETC50)
  sub_ETC50 <- c()
  for(i in main_ETC50){
    sub_ETC50 <- c(sub_ETC50, checkSubRule(i))
  }
  
  out_list <- list()
  for(i in sub_ETC50){
    out_list[length(out_list)+1] <- list(checkGenome(df, i, genome))
  }
  
  names(out_list) <- main_ETC50
  return(out_list)
}
setETC50(picr, "KE_lactobacillus_man_1")

main_ETC50 <- checkMainRule("ETC50")
check50(main_ETC50)

#pull index of evaluation priority
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

#pull individual Ids
pullIds <- function(string){
  out_vect <- c()
  for(i in resp_rules$parent){
    if(str_detect(string, i)){
      out_vect <- c(out_vect, i)
    }
  }
  if(is.null(out_vect)){
    out_vect <- c(out_vect, matchKO(string), matchEC(string))
    return(out_vect)
  }else{
    return(out_vect)
  }
}


#build list for rule evaluation priority
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
        priority_list[length(priority_list)+1] <- str_sub(rule_string, dbl_idx_strt[i], dbl_idx_end[i])
        dbl_names <- c(dbl_names, paste0("dbl", as.character(i)))
      }
    }
    
    priority_list[length(priority_list)+1] <- rule_string
    names(priority_list) <- c(first_names, dbl_names, "org")
    return(priority_list)
  }
}

evaluatePriorityExpr <- function(logic_list){
  logic_out <- logic_list
  for(id in seq_along(logic_list)){
    x <- logic_list[[id]]
    if(str_detect(x, "|")){ 
      x_or <- unlist(strsplit(x, split = "\\|")) #splits or expressions, returns chr vect of length 1 if no | present
      and_index <- c()
      for(ex1 in str_detect(x_or, "&")){ #record positions of & if present
        ifelse(ex1, and_index <- c(and_index, TRUE),and_index <- c(and_index, FALSE))
      }
      if(is.null(and_index)){
        ifelse(sum(as.logical(x_or)) > 0, #evaluate or expression and return T or F to logic list
               logic_out[[id]] <- TRUE,
               logic_out[[id]] <- FALSE)
      }else{
        x_and <- c()
        for(ex2 in x_or[and_index]){
          x_and <- c(x_and, unlist(strsplit(ex2, split = "\\&"))) #split & expressions
        }
        ifelse(sum(as.logical(x_and)) > 0, 
               x_or[and_index] <- TRUE, 
               x_or[and_index] <- FALSE) #evaluate and expressions and replace in x_or with T or F
        
        ifelse(sum(as.logical(x_or)) > 0,
               logic_out[[id]] <- TRUE,
               logic_out[[id]] <- FALSE) #evaluate or expression and return T or F to logic list
      }
    }
  }
  return(logic_out)
}

evaluateFunct <- function(genome_priority_list){
  id_list <- genome_priority_list[[1]]
  logic_list <- genome_priority_list[[2]]
  
  if(sum(names(id_list) %in% letters) > 0){ #checks if priority list has expressions to evaluate first
    geonome_priority_list_1 <- evaluatePriorityExpr(id_list, logic_list)
    id_list <- genome_priority_list_1[[1]]
    logic_list <- genome_priority_list_1[[2]]
  }
}



test <- checkGenome(as.data.frame(picr), checkSubRule("CytochromeCOxidaseCBB3"), "KE_lactobacillus_man_1")
logic_list <- test[[2]]
id_list <- test[[1]]




test_out <- evaluateExpr(id_list, logic_list)
test_out



checkGenome <- function(dataframe, rule_string, genome){
  rule_list <- buildPriorityList(buildIndex(rule_string), rule_string)
  
  ids_keep <- dataframe %>%
    filter(bin.id == genome) %>%
    pull(gene_id) %>%
    unique()
  
  df <- dataframe %>%
    filter(bin.id == genome,
           gene_id %in% ids_keep)
  
  rule_ids <- pullIds(rule_string)
  
  for(i in rule_ids){
    if(i %in% df$gene_id){
      next
    }else{
      rule_string <- str_replace_all(rule_string, i, "F")
    }
  }
  
  for(i in 1:nrow(df)){
    if(df[i ,"counts"] > 0){
      rule_string <- str_replace_all(rule_string, df[i, "gene_id"], "T") 
    }else if(df[i ,"counts"] == 0){
      rule_string <- str_replace_all(rule_string, df[i, "gene_id"], "F")
    }
  }
  priority_list <- buildPriorityList(buildIndex(rule_string), rule_string)
  return(list(rule_list, priority_list))
}



checkDbl <- function(priority_list){
  if(str_detect(names(priority_list), "Dbl")){
    return(TRUE)
  }else{
    return(FALSE)
  }
}



pullIds("[percent50CytochromeCOxidase|percent50CytochromeAa3600MenaquinolOxidase|percent50CytochromeOUbiquinolOxidase]&ETC50")
pullIds("K00376|EC1.7.2.4")

main_r <- checkMainRule("Aerobic")
main_idx <- buildIndex(main_r)
buildPriorityList(main_idx, main_r)


checkGenome(as.data.frame(picr), checkMainRule("Aerobic"), "KE_lactobacillus_man_1")
checkGenome(as.data.frame(picr), checkSubRule("CytochromeCOxidase"), "KE_lactobacillus_man_1")

for(i in unique(resp_rules$name)){
  main_rule <- checkMainRule(i)
  if(!is.na(check50(main_rule))){
    main_rule <- check50(main_rule)
    for(j in main_rule){
      genome <- checkGenome(as.data.frame(picr), checkSubRule(j), "KE_lactobacillus_man_1")
      print(genome)  
    }
  }else{
    for(j in pullIds(main_rule)){
      sub_rule <- checkSubRule(j)
      genome <- checkGenome(as.data.frame(picr), checkSubRule(j), "KE_lactobacillus_man_1")
      print(genome)  
    }
  }
}
check50(checkMainRule("SulfateReduction"))
test <- picr %>%
  filter(bin.id == "KE_lactobacillus_man_1")

"EC1.7.5.1" %in% test$gene_id

genome <- setRefClass("genome", fields = list(SulfateReduction = "ANY",
                                    FumarateReduction = "ANY",
                                    Aerobic = "ANY",
                                    TetrathionateReduction = "ANY",
                                    Microaerophillic = "ANY",
                                    DenitrificationNotETC = "ANY", 
                                    ETC50 = "ANY"))

as.name("t") <- genome()
rm(t)
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