library(tidyverse)

setwd(paste0("/Users/ikaialeleiwi/Desktop/Lab/Salmonella_NIH/Lactobacillus/",
             "Omics/Metagenome/rRNA/16S_From_All_Bins/MQHQ_Picrust2/",
             "MQHQ_Picrust2/"))

#Data 
picr_dram <- read_tsv("Clean_Data/picr_dram.tsv") %>%
  filter(total_copy_number_picrust > 0 & total_copy_number_dram > 0)

picr_dram_present_absent <- picr_dram %>%
  mutate(picrust_pa = ifelse(total_copy_number_picrust > 0,
                              yes = 1,
                              no = 0),
         dram_pa = ifelse(total_copy_number_dram > 0,
                          yes = 1,
                          no = 0))

expected_value <- factor(picr_dram_present_absent$dram_pa, 
                         levels = c(1, 0))
predicted_value <- factor(picr_dram_present_absent$picrust_pa, 
                          levels = c(1, 0))

#confustion matrix
m <- table(expected_value, predicted_value)
TP <- m[1,1] #true positives
TN <- m[2,2] #true negatives
FP <- m[2,1] #false positives
FN <- m[1,2] #false negatives
P <- TP + FN #real positives
N <- TN + FP #real neagatives
PP <- sum(picr_dram_present_absent$picrust_pa) #predicted positives
PN <- length(picr_dram_present_absent$picrust_pa) - PP #predicted negatives



PPV <- TP/PP #Positive predicted value (1 - FDR) (precision)
FDR <- FP/PP #False discovery rate (1 - PPV) 
FOR <- FN/PN #False omission rate (1 - NPV)
NPV <- TN/PN #Negative predictive value (1 - FOR)


RECALL <- TP/P
PRECISION <- TP/(TP+FP)
ACCURACY <- (TP + TN)/(TP + TN + FP + FN)
FMEASURE <- (2*RECALL*PRECISION)/(RECALL+PRECISION)

#match output to metabolisms


