library(tidyverse)
library(vegan)
library(ggtext)

setwd(paste0("/Users/ikaialeleiwi/Desktop/Lab/Salmonella_NIH/Lactobacillus/",
             "Omics/Metagenome/rRNA/16S_From_All_Bins/ASVs_to_Bins_MMSEQS2/"))
##Data
nmds_plot_df <- read_tsv("Clean_Data/nmds_plot_df.tsv")


#create hulls for treatment polygons
find_hull <- function(df) df[chull(df$NMDS1, df$NMDS2), ]
hulls <- plyr::ddply(nmds_plot_df, "method", find_hull)

# Bray Curtis NMDS
bray_nmds <- nmds_plot_df %>%
  ggplot(aes(x = NMDS1, 
             y = NMDS2, 
             fill = method)) +
  geom_point(aes(fill = method),
             shape = 23) +
  geom_polygon(data = hulls,
               aes(x = NMDS1, 
                   y = NMDS2),
               alpha = 0.3,
               size = .5) +
  geom_textbox(data = data.frame(NMDS1 = .5,
                                 NMDS2 = -1.75,
                                 method = NA),
               aes(label = "Anosim
               \nR: 0.01592
                   \nSignificance: 0.001"),
               box.color = "black",
               box.size = 0.3,
               fill = "white",
               height = NULL,
               width = NULL) +
  theme_classic() +
  scale_alpha_identity() +
  scale_size_identity(guide = "legend") +
  labs(title = "MetaG KO's and ASV Predicted KO's",
       fill="Method") 


jpeg(file="Figs/picrust_dram_KO.jpeg")
bray_nmds
dev.off()
