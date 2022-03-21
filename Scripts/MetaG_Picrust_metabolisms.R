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