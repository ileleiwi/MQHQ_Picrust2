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

#functions
source("Scripts/rule_functions.R")


#checks rule structure and gets id's or sub funct
checkRule <- function(funct){
  switch(funct,
         SulfateReduction = "Dsr",
         Dsr = "K11180|EC1.8.99.5|K11181|EC1.8.99.5",
         FumarateReduction = "K00244|K00245|K00246|K00247",
         Aerobic = "[percent50CytochromeCOxidase|percent50CytochromeAa3600MenaquinolOxidase|percent50CytochromeOUbiquinolOxidase]&ETC50",
         CytochromeCOxidase = "K02275,K02274,K02276|K15408,K02277",
         CytochromeAa3600MenaquinolOxidase = "K02827,K02826,K02828,K02829",
         CytochromeOUbiquinolOxidase = "K02297,K02298,K02299,K02300",
         TetrathionateReduction = "[ttrrs&1ofttrabc]|2ofttrabc",
         ttrrs = "K13040|K13041",
         ttrabc = "K08359,K08358,K08357",
         Microaerophillic = "[percent50CytochromeBDUbiquinolOxidase|percent50CytochromeCOxidaseCBB3]&ETC50",
         CytochromeBDUbiquinolOxidase = "K00425,K00426,K00424|K22501",
         CytochromeCOxidaseCBB3 = "K00404|K00405|K15862,K00407,K00406",
         DenitrificationNotETC = "notETC50&[NitrateReductase|NitriteReductase|NitricOxideReductase|NitrousOxideReductase]",
         NitrateReductase = "K00370|K00371|K00374|K02567|K02568|EC1.7.5.1|EC1.9.6.1",
         NitriteReductase = "K00368|K15864|EC1.7.2.1",
         NitricOxideReductase = "K04561|K02305|EC1.7.2.5",
         NitrousOxideReductase = "K00376|EC1.7.2.4",
         ETC50 = "percent50ComplexIA|percent50ComplexIB|percent50ComplexIC",
         ComplexIA = "K00330,[K00331&K00332&[K00333|K00331]&[K13378|K13380]],K00334,K00335,K00336,K00337,K00338,K00339,K00340,[K00341&K00342|K15863],K00343",
         ComplexIB = "K05574,K05582,K05581,K05579,K05572,K05580,K05578,K05576,K05577,K05575,K05573,K05583,K05584,K05585",
         ComplexIC = "K03945,K03946,K03947,K03948,K03949,K03950,K03951,K03952,K03953,K03954,K03955,K03956,K11352,K11353",
         NA)
}


#evaluates string from rule structure and outputs rule list
evaluateRule <- function(rule_df, rule_out){
  if(is.na(rule_out)){
    return(NA)
  }else{
    if(!hasIds(rule_out)){ 
      ids <- pullIds(rule_out)
      rule_list <- map(ids, checkRule)
      names(rule_list) <- ids
    }else{
      rule_list <- as.list(rule_out)
      names(rule_list) <- rule_df[rule_df$child == rule_out, "parent"]
    }
    
    #handle ETC50
    for(rule in names(rule_list)){
      if(hasIds(rule_list[rule])){
        next
      }else{
        rule_list[rule] <- list(evaluateRule(rule_df, rule_list[rule]))
      }
    }
    
    #add original rule to list
    original_names <- names(rule_list)
    rule_list[length(rule_list)+1] <- rule_out
    names(rule_list) <- c(original_names, "org")
    
  return(rule_list)
  }
}


#checks rule list against counts dataframe and provides T/F list for given genome
buildLogicList <- function(dataframe, rule_list, genome){

  logic_list <- rule_list
  #check for nestedness
  if(listOfLists(rule_list)){ 
    nested_index <- str_detect(map_chr(rule_list, class), "list")
    nested_list <- pluck(rule_list,which(nested_index))
    logic_list[[which(nested_index)]] <- buildLogicList(dataframe, 
                                                        nested_list, 
                                                        genome)
  }
  
  #get unique IDs
  ids <- unlist(map(unlist(rule_list), pullIds))
  names(ids) <- NULL
  ids <- unique(c(matchKO(ids), matchEC(ids)))
  
  df <- dataframe %>%
    filter(bin.id == genome)
  
  for(id in ids){
    if(!(id %in% df[,"gene_id"])){
      logic_list <- replaceIdInList(logic_list, id, repl = "F")
    }else if(df[which(df[,"gene_id"] == id),"counts"] > 0){
      logic_list <- replaceIdInList(logic_list, id, repl = "T")
    }else if(df[which(df[,"gene_id"] == id),"counts"] == 0){
      logic_list <- replaceIdInList(logic_list, id, repl = "F")
    }
  }

return(logic_list)
}




#build index of evaluation priority (run on individual elements of rule_list)
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

#replaces DBL instances in org string
replaceDBLs <- function(priority_list, dbl_names){
  out_list <- priority_list
  for(i in dbl_names){
    dbl_repl <- str_extract(string = out_list[["org"]], pattern = "\\[.*\\]" )
    out_list[[i]] <- dbl_repl
    out_list[[i]] <- str_remove_all(out_list[[i]], "\\[")
    out_list[[i]] <- str_remove_all(out_list[[i]], "\\]")
    
    out_list[["org"]] <- str_replace(out_list[["org"]], 
                                     out_list[[i]], i)
    
    out_list[["org"]] <-str_remove_all(out_list[["org"]], "\\[")
    out_list[["org"]] <-str_remove_all(out_list[["org"]], "\\]")
  }
  return(out_list)
}

#build list for rule evaluation priority (run on individual elements of rule_list)
buildPriorityList <- function(index, rule_string){
  priority_list <- list()
  dbl_idx_strt <- c()
  dbl_idx_end <- c()
  org <- rule_string
  
  if(is.null(index)){
    priority_list[length(priority_list)+1] <- rule_string
    names(priority_list) <- "org"
    return(priority_list)
  }else{
    let_cnt <- 1
    repl_idx <- abs(index)
    for(i in 2:length(index)){
      if(index[i] > 0 & index[i-1] > 0){
        dbl_idx_strt <- c(dbl_idx_strt, index[i-1])
      }else if(index[i] < 0 & index[i-1] > 0){
        x <- unlist(str_sub(rule_string, index[i-1], -index[i]))
        priority_list[length(priority_list)+1] <- str_remove_all(x, "\\[|\\]")
        str_sub(string = org, 
                start = repl_idx[i-1], 
                end = repl_idx[i]) <- letters[let_cnt]
        let_cnt <- let_cnt + 1
        repl_idx <- repl_idx - nchar(str_sub(string = org, 
                                           start = repl_idx[i-1], 
                                           end = repl_idx[i])) + 1 
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
        dbl_names <- c(dbl_names, paste0("DBL", as.character(i)))
      }
    }
    
    
    priority_list[length(priority_list)+1] <- org
    names(priority_list) <- c(first_names, dbl_names, "org")
    priority_list <- replaceDBLs(priority_list, dbl_names)
    
    return(priority_list)
  }
}



#collapses priority logic list: used in evaluateLogicList()
collapsePriority <- function(logic_list){
  collapsed_list <- logic_list
  elems_rm <- c()
  for(elem in names(collapsed_list)){
    if(is.logical(collapsed_list[[elem]])){
      for(elem_idx in seq_along(collapsed_list)){
        if(str_detect(collapsed_list[[elem_idx]], elem)){
          repl_char <- ifelse(collapsed_list[[elem]], "T", "F")
          collapsed_list[[elem_idx]] <- str_replace(collapsed_list[[elem_idx]], 
                                             elem, 
                                             repl_char)
          elems_rm <- c(elems_rm, elem)
        }
      }
    }else{
      next
      }
  }  
  for(i in elems_rm){
    collapsed_list[[i]] <- NULL
  }
  return(collapsed_list)
}



#evaluates | and & logic
evalLogic <- function(logic_in){
  logic_out <- logic_in
  for(idx in seq_along(logic_out)){
    if(str_detect(logic_out[[idx]], "|")){ 
      x_or <- unlist(strsplit(logic_out[[idx]], split = "\\|")) #splits or expressions, returns chr vect of length 1 if no | present
      and_index <- c()
      for(ex1 in str_detect(x_or, "&")){ #record positions of & if present
        ifelse(ex1, and_index <- c(and_index, TRUE),and_index <- c(and_index, FALSE))
      }
      if(is.null(and_index)){
        ifelse(sum(as.logical(x_or)) > 0, #evaluate or expression and return T or F to logic list
               logic_out[[idx]] <- TRUE,
               logic_out[[idx]] <- FALSE)
      }else{
        x_and <- c()
        for(ex2 in x_or[and_index]){
          x_and <- c(x_and, unlist(strsplit(ex2, split = "\\&"))) #split & expressions
        }
        ifelse(sum(as.logical(x_and)) == length(x_and), 
               x_or[and_index] <- TRUE, 
               x_or[and_index] <- FALSE) #evaluate and expressions and replace in x_or with T or F
        
        ifelse(sum(as.logical(x_or)) > 0,
               logic_out[[idx]] <- TRUE,
               logic_out[[idx]] <- FALSE) #evaluate or expression and return T or F to logic list
      }
    }else if(str_detect(logic_out[[idx]], "&")){
      x_and <- unlist(strsplit(logic_out[[idx]], split = "\\&")) #splits and expressions, returns chr vect of length 1 if no & present
      ifelse(sum(as.logical(x_and)) == length(x_and),
               logic_out[[idx]] <- TRUE,
               logic_out[[idx]] <- FALSE) #evaluate and expression and return T or F to logic list
    }
  }
  return(logic_out)
}



#returns percent true over percent false in comma separated logic string
commaPercent <- function(string){
  vect <- unlist(strsplit(string, split = ","))
  log_vect <- as.logical(evalLogic(vect))
  trues <- length(which(log_vect))
  falses <- length(which(!log_vect))
  
  percT <- round(trues/length(log_vect)*100, 0)
  percF <- round(falses/length(log_vect)*100, 0)
  return(paste0(percT, '/', percF))
}

t <- "F,x,T,T,T,F,F,F,F,F,F"
commaPercent(t)

#evaluates a logic list to level of original logic string
evaluateLogicList <- function(priority_list_logic){
  logic_out <- priority_list_logic
  #tests for complete list evaluation  
  if(length(logic_out) > 1){
    logic_out <- collapsePriority(evalLogic(logic_out))
    
  }else if(!is.logical(logic_out[["org"]])){
    logic_out[["org"]] <- commaPercent(logic_out[["org"]])
  }
  return(logic_out)
}


evaluateFunction(test2)
collapsePriority(ETC50_test)

test <- evaluateRule(resp_rules, checkRule("TetrathionateReduction"))
test <- evaluateRule(resp_rules, checkRule("Aerobic"))
test <- evaluateRule(resp_rules, checkRule("DenitrificationNotETC"))
test <- evaluateRule(resp_rules, checkRule("Microaerophillic"))


test2 <- buildLogicList(picr,test,"LM_bin.44")

#for each logic element in logic_list need to build priority index and priority_list then evaluate logic
elem_index <- buildIndex(test2[[1]])
test3 <- buildPriorityList(elem_index, test2[[1]])
test4 <- evaluateLogicList(test3)

#replace element in test2 with evaluated logic
test2[[1]] <- test4$org

ETC_testorg <- test2[[which(map_lgl(test2, is.list))]]

ETC_test <- test2[[which(map_lgl(test2, is.list))]]
c <- 0
for(element in seq_along(ETC_test)){
  elem_index <- buildIndex(ETC_test[[element]])
  prir_list <- buildPriorityList(elem_index, ETC_test[[element]])
  ETC_test[[element]] <- evaluateLogicList(prir_list)
  c <- c+1
  print(c)
  print(prir_list)
  print(ETC_test[[element]])
  #print(ETC_test)
}
ETC_test

ETC50_test <- test2[[4]]
evalLogic(collapsePriority(evalLogic(ETC50_test)))
collapsePriority(evalLogic(collapsePriority(ETC50_test)))

elem_index <- buildIndex(ETC_test[[1]])
prir_list <- buildPriorityList(elem_index, ETC_test[[1]])
#out_elem <- evaluateLogicList(prir_list)
collapsePriority(evalLogic(collapsePriority(evalLogic(collapsePriority(prir_list)))))
collapsePriority(evalLogic(collapsePriority(prir_list)))

evalLogic(x)

x <- which(map_lgl(test2, is.list))
test2[[which(map_lgl(test2, is.list))]]
evaluateFunction(ETC_test)

evaluateFunction <- function(logic_list){
  out_list <- logic_list #create output list
  sub_list <- logic_list[[which(map_lgl(logic_list, is.list))]] #save sub list as separate list
  non_list_loop_idx <- which(!map_lgl(out_list, is.list)) #save original index of main list
  out_list[[which(map_lgl(out_list, is.list))]] <- NULL #remove sublist from main list
 
  
  for(element in seq_along(sub_list)){
    elem_index <- buildIndex(sub_list[[element]])
    prir_list <-buildPriorityList(elem_index, sub_list[[element]])
    sub_list[[element]] <- evaluateLogicList(prir_list)
    
    for(chr in names(sub_list[[element]])){
      if(str_detect(chr, "DBL")){
        sub_list[[element]][[chr]]  <- evalLogic(sub_list[[element]][[chr]])
        sub_list[[element]] <- evaluateLogicList(sub_list[[element]])
      }
    }
  }
  out_list[[which(map_lgl(logic_list, is.list))]] <- unlist(sub_list)
  
  for(element in non_list_loop_idx){
      elem_index <- buildIndex(out_list[[element]])
      prir_list <-buildPriorityList(elem_index, out_list[[element]])
      out_list[[element]] <- unlist(evaluateLogicList(prir_list))
  }
  return(out_list)
}

logic_list <- test2
out_list <- logic_list #create output list
sub_list <- logic_list[[which(map_lgl(logic_list, is.list))]] #save sub list as separate list
sub_list_idx <- which(map_lgl(logic_list, is.list)) #save sub list position in main list
non_list_loop_idx <- which(!map_lgl(out_list, is.list)) #save original index of main list
out_list[[which(map_lgl(out_list, is.list))]] <- NULL #remove sublist from main list


for(element in seq_along(sub_list)){
  elem_index <- buildIndex(sub_list[[element]])
  prir_list <-buildPriorityList(elem_index, sub_list[[element]])
  for(chr in names(prir_list)){
    if(str_detect(chr, "DBL")){
      prir_list[[chr]]  <- evalLogic(prir_list[[chr]])
      prir_list <- evaluateLogicList(evaluateLogicList(prir_list))
    }
  }
  sub_list[[element]] <- evaluateLogicList(prir_list)
}

out_list[[which(map_lgl(logic_list, is.list))]] <- unlist(sub_list)
#########need to rebuild sublist here by adding back ETC50 in position 4
##########need to make function that evaluates percent string and collapses ETC50 sublist to T or F



elem_index <- buildIndex(sub_list[["ComplexIA"]])
tp <- buildPriorityList(elem_index, sub_list[["ComplexIA"]])
tp1 <- evaluateLogicList(tp)


for(chr in names(tp1)){
  if(str_detect(chr, "DBL")){
    tp1[[chr]]  <- evalLogic(tp1[[chr]])
    tp1 <- evaluateLogicList(tp1)
  }
}



for(i in seq_along(sub_list)){
  print(i)
  for(chr in sub_list[[i]]){
    print(chr)
  }
}






genome <- setRefClass("genome", fields = list(genome_name = "",
                                              SulfateReduction = "ANY",
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