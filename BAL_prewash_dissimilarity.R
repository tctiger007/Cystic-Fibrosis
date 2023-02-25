library(dplyr)
library(ggplot2)
library(phyloseq)
library(microbiome)
library(vegan)
library(wesanderson)

dirname = dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(dirname)
phy = readRDS("../Data/phy_noncontam.rds")
meta = meta(phy)
tax = tax_table(phy)
otu = otu_table(phy)

prewash.phy = subset_samples(phy, Sample_type == "Prewash")
BAL.phy = subset_samples(phy, Sample_type == "BAL")
prewash.meta = meta(prewash.phy)
BAL.meta = meta(BAL.phy)
BAL.BA2 = meta %>% dplyr::filter(Sample_type == "BAL" & Visit == "BA2") 
BAL.BA7 = meta %>% dplyr::filter(Sample_type == "BAL" & Visit == "BA7") 
BA2_BA7_ID = intersect(BAL.BA2$Subject_ID, BAL.BA7$Subject_ID)  #45
BAL.meta = BAL.meta %>% dplyr::filter(Subject_ID %in% BA2_BA7_ID)
prewash.meta = prewash.meta %>% dplyr::filter(Subject_ID %in% BA2_BA7_ID)
commonID = intersect(BAL.meta$Sample_ID, prewash.meta$Sample_ID)
phy = subset_samples(phy, Sample_ID %in% commonID)
# bray = phyloseq::distance(phy, method = "bray")
# bray = as.matrix(bray)
relab_genera = transform_sample_counts(phy, function(x) x / sum(x) * 100) 
ord = ordinate(relab_genera, method="PCoA", distance = "bray")

png("Figures/BAL_prewash_ordinate.png", width = 6, height = 6, unit = "in", res = 300)
ordinate = plot_ordination(relab_genera, ord, color = "Sample_type") + 
  geom_point(size=4) + 
  scale_color_manual(values = wes_palette("FantasticFox1")[c(4,3)]) +
  stat_ellipse() + 
  theme_classic()
dev.off()


