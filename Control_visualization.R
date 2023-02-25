library(dplyr)
library(ggplot2)
library(readxl)
library(stringr)
library(factoextra)
library(FactoMineR)
library(Rtsne)
library(phyloseq)
library(DESeq2)
library(MASS)
library(phyloseq)
library(scales)
library(reshape2)
library(forcats)
library(RColorBrewer)
library(FSA)
library(pheatmap)
library(cluster)
library(gridExtra)
library(grid)
library(microViz)
library(alluvial)
library(ggalluvial)
library(ape)

dirname = dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(dirname)
BAL.phy = readRDS("../../Data/phy_noncontam.rds")
# only keep the controls 
BAL.phy = subset_samples(BAL.phy, Sample_type == "Prewash")
meta = BAL.phy@sam_data
meta$Visit[which(meta$Visit == "")] = "IASM"
meta$Sequence_ID = rownames(meta)
## Remove IASM samples
BAL.phy = subset_samples(BAL.phy, Visit %in% c("BA2", "BA7"))
meta = BAL.phy@sam_data
orderID = sort(meta$Sample_ID)[c(1:14, 55:56, 15:54)]

######################## Part 1: visualization of relative abundance ########################
## i. genus level 
BAL.phy_genus <- tax_glom(BAL.phy, "Genus", NArm = TRUE)
# Get top 20 genera
orderTaxa = names(sort(taxa_sums(BAL.phy_genus), decreasing=TRUE))

# apply(BAL.phy_genus_relabun@otu_table, 2, summary)
# sum(rowSums(BAL.phy@otu_table)/sum(rowSums(BAL.phy@otu_table)) > 0.02)

top20_genera <- orderTaxa[1:20]
# Transform Taxa counts to relative abundance
BAL.phy_genus_relabun <- transform_sample_counts(BAL.phy_genus, function(OTU) OTU/sum(OTU) * 100)

# Extract the top 20 taxa 
BAL.phy_genus_top20 <- prune_taxa(top20_genera, BAL.phy_genus_relabun)

# load the color for the genus on ggplot
load("color_for_genus.RData")
load("color_for_phylum.RData")
cols_genus = as.data.frame(cols_genus)
cols_phylum = as.data.frame(cols_phylum)

#### showing the top 20 taxa separately, and the rest in gray 
BAL.phy.df = psmelt(BAL.phy_genus_relabun)  %>% 
  mutate(name = fct_relevel(Sample_ID, orderID)) 
BAL.phy.df$Genus[which(BAL.phy.df$OTU %in% top20_genera == FALSE)] = "Low Abundance"
top20_genera_name = BAL.phy.df$Genus[match(top20_genera, BAL.phy.df$OTU)]
BAL.phy.df = BAL.phy.df %>% mutate(Genusorder = fct_relevel(Genus, c(top20_genera_name, "Low Abundance")))
BAL.phy.df$Country = factor(BAL.phy.df$Country, levels = c("USA", "AUS"))
BAL.phy.df$Group = paste0(BAL.phy.df$Country, " ", BAL.phy.df$Visit)
BAL.phy.df$Group = factor(BAL.phy.df$Group, levels = c("USA BA2","AUS BA2",  
                                                       "USA BA7", "AUS BA7"))
BAL.phy_genus_relabun_BA2 = subset_samples(BAL.phy_genus_relabun, Visit == "BA2")
BAL.phy_genus_relabun_BA7 = subset_samples(BAL.phy_genus_relabun, Visit == "BA7")
# write.csv(as.matrix(BAL.phy_genus_relabun@otu_table), "../Data/BAL_genus_relabun.csv")
# write.csv(as.matrix(BAL.phy_genus_relabun_BA2@otu_table), "../Data/BAL_genus_relabun_BA2.csv")
# write.csv(as.matrix(BAL.phy_genus_relabun_BA7@otu_table), "../Data/BAL_genus_relabun_BA7.csv")
BAL.phy_genus_relabun_BA2_OTU = as.matrix(BAL.phy_genus_relabun_BA2@otu_table)


# try to keep the genus that are present in more than 15% of the samples 
otu = otu_table(BAL.phy_genus_relabun)
otu.df = as.data.frame(otu)
otu.df = as.data.frame(t(otu.df))
otu.df$`Study ID` = rownames(otu.df)
rownames(otu.df) = NULL
dim(otu) #352 109
otu = t(otu)
otu = as.data.frame(otu)
zero.df = as.data.frame(lapply(otu, function(x){ sum(x==0)}))
zero.df = as.data.frame(t(zero.df))
nonzero.df = data.frame(nonzero = dim(otu)[1]-zero.df$V1)
rownames(nonzero.df) = rownames(zero.df)
rownames(nonzero.df) = tax_table(BAL.phy)[match(rownames(nonzero.df), 
                                                rownames(otu_table(BAL.phy))), "Genus"]
summary(nonzero.df$nonzero)
hist(nonzero.df$nonzero)
# keep the genera that were prevalent in more than 10% of the samples 
keep = rownames(nonzero.df)[which(nonzero.df$nonzero > dim(otu)[1] * 0.1)]
keep

commonID.T1T2 = c("IU120", "IU121", "IU123", "IU124", "IU125", "IU127", "IU128", "IU129",
                  "IU130", "IU131", "M134",  "M135",  "M136", "M137",  "M140",  "M141", 
                  "M151",  "M153",  "M155" )

# View(BAL.phy.df)
BAL.phy.df$Color = cols_genus$Color[match(BAL.phy.df$Genus, cols_genus$Genus)]
# genus in prewash but not in samples 
sum(is.na(BAL.phy.df$Color))
levels(as.factor(BAL.phy.df$Genus[which(is.na(BAL.phy.df$Color))])) #11 genus in prewash not in samples 
# create 11 more colors for the 11 genus in prewash but not in samples 
set.seed(111)
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
cols1 = sample(col_vector, 31)
cols_genus1 = cbind(levels(as.factor(BAL.phy.df$Genus[which(is.na(BAL.phy.df$Color))])),
                    sample(cols1[cols1 %in% cols_genus$Color == FALSE], 11, replace = F))
colnames(cols_genus1) = c("Genus", "Color")
cols_genus1 = as.data.frame(cols_genus1)
cols_genus1 = rbind.data.frame(cols_genus, cols_genus1)
BAL.phy.df$Color = cols_genus1$Color[match(BAL.phy.df$Genus, cols_genus1$Genus)]
cols1 = cols_genus1$Color[match(levels(BAL.phy.df$Genusorder), cols_genus1$Genus)]

png("prewash_genus_relative_abundance.png", width = 8, height = 10, units = "in", res = 300)
BAL.phy.df %>% dplyr::filter(Subject_ID %in% commonID.T1T2) %>% 
  ggplot(aes(x =name, y = Abundance, fill = Genusorder)) +
  geom_bar(stat = "identity", width = .95) +
  # geom_col(position = position_stack(reverse = TRUE)) + 
  labs(x = "",
       y = "Relative Abundance",
       title = "Genus Relative Abundance") +
  scale_fill_manual(values = cols1) +
  facet_grid(Country ~ Visit, scales = "free_x", space = "free_x") +
  labs(fill = "Genus") + 
  theme_classic() + 
  theme(axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, hjust = 1),
        axis.text.y = element_text(size = 12),
        legend.text = element_text(size = 14),
        strip.text = element_text(size = 14),
        # strip.background = element_blank(),
        axis.title.y = element_text(size = 16),
        plot.title = element_text(size = 16)) 
dev.off()



## ii. phylum level relative abundance 
BAL.phy_phylum <- tax_glom(BAL.phy, "Phylum", NArm = TRUE)
# Get top 10 genera
orderPhylum = names(sort(taxa_sums(BAL.phy_phylum), decreasing=TRUE))
top10_phylum <- orderPhylum [1:10]
# Transform Taxa counts to relative abundance
BAL.phy_phylum_relabun <- transform_sample_counts(BAL.phy_phylum, function(OTU) OTU/sum(OTU) * 100)

# Extract the top 10 taxa 
BAL.phy_phylum_top10 <- prune_taxa(top10_phylum, BAL.phy_phylum_relabun)

#### showing the top 10 taxa separately, and the rest in gray 
BAL.phy.df1 = psmelt(BAL.phy_phylum_relabun)  %>% 
  mutate(name = fct_relevel(Sample_ID, orderID)) 
BAL.phy.df1$Phylum[which(BAL.phy.df1$OTU %in% top10_phylum == FALSE)] = "Low Abundance"
top10_phylum_name = BAL.phy.df1$Phylum[match(top10_phylum, BAL.phy.df1$OTU)]
BAL.phy.df1 = BAL.phy.df1 %>% mutate(phylumorder = fct_relevel(Phylum, c(top10_phylum_name, "Low Abundance")))
BAL.phy.df1$Country = factor(BAL.phy.df1$Country, levels = c("USA", "AUS"))
BAL.phy.df1$Group = paste0(BAL.phy.df1$Country, " ", BAL.phy.df1$Visit)
BAL.phy.df1$Group = factor(BAL.phy.df1$Group, levels = c("USA BA2","AUS BA2",  
                                                         "USA BA7", "AUS BA7"))

cols_phylum = rbind.data.frame(cols_phylum, 
                               c("Chloroflexi", "#F0027F"))
colnames(cols_phylum) = c("Phylum", "Color")
cols2 = cols_phylum$Color[match(levels(BAL.phy.df1$phylumorder), 
                                cols_phylum$Phylum)]

png("prewash_phylum_relative_abundance.png", width = 8, height = 10, units = "in", res = 300)
BAL.phy.df1 %>% dplyr::filter(Subject_ID %in% commonID.T1T2) %>% 
  ggplot(aes(x =name, y = Abundance, fill = phylumorder)) +
  geom_bar(stat = "identity") +
  # geom_col(position = position_stack(reverse = TRUE)) + 
  labs(x = "",
       y = "Relative Abundance",
       title = "Phylum Relative Abundance") +
  scale_fill_manual(values = cols2) + 
  facet_grid(Country ~ Visit, scales = "free_x", space = "free_x") +
  labs(fill = "Phylum") + 
  theme_classic() + 
  theme(axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, hjust = 1),
        axis.text.y = element_text(size = 12),
        legend.text = element_text(size = 10),
        strip.text = element_text(size = 12)) 
dev.off()
# 
# BAL.phy_species <- tax_glom(BAL.phy, "Species", NArm = TRUE)
# # Get top 20 genera
# orderSpecies = names(sort(taxa_sums(BAL.phy_species), decreasing=TRUE))
# top20_species <- orderSpecies[1:20]
# # Transform Taxa counts to relative abundance
# BAL.phy_species_relabun <- transform_sample_counts(BAL.phy_species, function(OTU) OTU/sum(OTU) * 100)
# 
# # Extract the top 20 taxa 
# BAL.phy_species_top20 <- prune_taxa(top20_species, BAL.phy_species_relabun)
# 
# 
# #### showing the top 20 taxa separately, and the rest in gray 
# BAL.phy.df2 = psmelt(BAL.phy_species_relabun)  %>% 
#   mutate(name = fct_relevel(Sample_ID, orderID))
# BAL.phy.df2$Species = paste0(BAL.phy.df2$Genus, " ", BAL.phy.df2$Species)
# BAL.phy.df2$Species[which(BAL.phy.df2$OTU %in% top20_species == FALSE)] = "Low Abundance"
# top20_species_name = BAL.phy.df2$Species[match(top20_species, BAL.phy.df2$OTU)]
# BAL.phy.df2 = BAL.phy.df2 %>% mutate(Speciesorder = fct_relevel(Species, c(top20_species_name, "Low Abundance")))
# BAL.phy.df2$Country = factor(BAL.phy.df2$Country, levels = c("USA", "AUS"))
# BAL.phy.df2$Group = paste0(BAL.phy.df2$Country, " ", BAL.phy.df2$Visit)
# BAL.phy.df2$Group = factor(BAL.phy.df2$Group, levels = c("USA BA2","AUS BA2",  
#                                                          "USA BA7", "AUS BA7"))
# 
# # png("../Results/species_relative_abundance.png", width = 20, height = 10, units = "in", res = 300)
# BAL.phy.df2 %>% 
#   # mutate(`Species Name` = paste0(BAL.phy.df2$Genus, BAL.phy.df2$Species)) %>%
#   ggplot(aes(x =name, y = Abundance, fill = Speciesorder)) +
#   geom_bar(stat = "identity") +
#   geom_col(position = position_stack(reverse = TRUE)) + 
#   labs(x = "",
#        y = "Relative Abundance",
#        title = "Species Relative Abundance") +
#   scale_fill_manual(values = c(cols, "gray57")) + 
#   facet_grid(Country ~ Visit, scales = "free") +
#   labs(fill = "Species") + 
#   theme_classic() + 
#   theme(axis.text.x = element_text(size = 7, angle = 90, vjust = 0.5, hjust = 1),
#         axis.text.y = element_text(size = 12),
#         legend.text = element_text(size = 10),
#         strip.text = element_text(size = 12)) 
# # dev.off()
