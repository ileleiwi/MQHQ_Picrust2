library(tidyverse)


setwd(paste0("/Users/ikaialeleiwi/Desktop/Lab/Salmonella_NIH/Lactobacillus/",
             "Omics/Metagenome/rRNA/16S_From_All_Bins/MQHQ_Picrust2/",
             "MQHQ_Picrust2/"))

#functions to parse through rules

#match EC number format EC-.-.-.-
matchEC <- function(string){
  x <- str_extract_all(string, 
                       "EC[\\d-]{1,2}\\.[\\d-]{1,2}\\.[\\d-]{1,2}\\.[\\d-]{1,2}")
  return(unlist(x))
}


#match KO number format K-----
matchKO <- function(string){
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
checkOne <- function(string){
  if(str_detect(string, "1of")){
    loc_1 <- unlist(str_match_all(string, "1of(.+?)[\\|,\\]].*$"))
    return(loc_1[!str_detect(loc_1, "1of")])
  }else{
    return(NA)
  }
}


#match if 2of  marker is present
checkTwo <- function(string){
  if(str_detect(string, "2of")){
    loc_2 <- unlist(str_match_all(string, "^.*2of(.+?)$"))
    return(loc_2[!str_detect(loc_2, "2of")])
  }else{
    return(NA)
  }
}

#check if string has KO or EC ids
hasIds <- function(string){
  return(length(matchEC(string)) + length(matchKO(string)) > 0)
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


#check if list has another list as an element
listOfLists <- function(list){
  return(sum(str_detect(map_chr(list, class), "list")) > 0)
}

#replaces all instances of an id in a list with repl
replaceIdInList <- function(list, id, repl){
  list_out <- list
  if(listOfLists(list)){
    idx <- str_detect(map_chr(list, class), "list")
    for(elem in idx){
      if(elem){
        new_nested <- list_out[[which(idx)]]
        for(elem1 in seq_along(new_nested)){
          new_nested[[elem1]] <- gsub(pattern = id, 
                                      x = new_nested[[elem1]], 
                                      replacement = repl)
        }
        list_out[[which(idx)]] <- new_nested
      }else{
        for(elem1 in which(!idx)){
          list_out[[elem1]] <- gsub(pattern = id, 
                                    x = list_out[[elem1]], 
                                    replacement = repl)
        }
      }
    }
  }else{
    for(elem1 in seq_along(list_out)){
      list_out[[elem1]] <- gsub(pattern = id, 
                                x = list_out[[elem1]], 
                                replacement = repl)
    }
  }
  return(list_out)
}

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

##buildPriorityList output types:
#list - 1 entry named org, string is subrule