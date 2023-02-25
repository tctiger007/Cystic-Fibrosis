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
library(microbiome)
library(wesanderson)


dirname = dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(dirname)
phy = readRDS("../../Data/phy_noncontam.rds")
meta = meta(phy)
tax = tax_table(phy)
otu = otu_table(phy)

# meta$Visit = ifelse(meta$Visit == "BA2", "T1", 
#                     ifelse(meta$Visit == "BA7", "T2", ""))
BAL.BA2 = meta %>% dplyr::filter(Sample_type == "BAL" & Visit == "BA2") #51
BAL.BA7 = meta %>% dplyr::filter(Sample_type == "BAL" & Visit == "BA7") #58
dim(meta %>% dplyr::filter(Study == "ECF" & Sample_type == "BAL" & Country == "AUS" ))
dim(meta %>% dplyr::filter(Study == "ECF" & Sample_type == "BAL" & Country == "USA" ))

BA2_BA7_ID = intersect(BAL.BA2$Subject_ID, BAL.BA7$Subject_ID) 
# 45; we keep these 45 matched BAL samples. 
# Now we find the matched prewash samples 

BA2_BA7_ID1 = c(paste0(BA2_BA7_ID, "BA2"), 
                paste0(BA2_BA7_ID, "BA7"))
prewash = meta %>% dplyr::filter(Sample_type == "Prewash")
sum(prewash$Sample_ID %in% paste0(BA2_BA7_ID, "BA2")) #21
sum(prewash$Sample_ID %in% paste0(BA2_BA7_ID, "BA7")) #30
  
prewash.BA2 = prewash %>% dplyr::filter(Sample_ID %in% paste0(BA2_BA7_ID, "BA2"))
prewash.BA7 = prewash %>% dplyr::filter(Sample_ID %in% paste0(BA2_BA7_ID, "BA7"))
# prewash_ID = c(intersect(BAL.BA2$Subject_ID, prewash.BA2$Subject_ID), #22
#                intersect(BAL.BA7$Subject_ID, prewash.BA7$Subject_ID)) #33
prewash_ID = c(prewash.BA2$Sample_ID, prewash.BA7$Sample_ID)

# BAL phyloseq object that contains matched BA2 and BA7 samples 
BAL.phy = subset_samples(phy, Subject_ID %in% BA2_BA7_ID)
BAL.phy = subset_samples(BAL.phy, Sample_type == "BAL")
BAL.meta = meta(BAL.phy)
# prewash phyloseq objects 
# prewash.phy = subset_samples(phy, Subject_ID %in% prewash_ID)
# prewash.phy = subset_samples(prewash.phy, Sample_type == "Prewash")
# prewash.meta = meta(prewash.phy)

# 1. visualization of the microbiome in BAL samples and matched prewash samples
#### also plot qPCR results
## 1). we plot BAL samples 
## i. genus level 
BAL.phy_genus <- tax_glom(BAL.phy, "Genus", NArm = TRUE)
# Get top 20 genera
orderTaxa = names(sort(taxa_sums(BAL.phy_genus), decreasing=TRUE))

top20_genera <- orderTaxa[1:20]
# Transform Taxa counts to relative abundance
BAL.phy_genus_relabun <- transform_sample_counts(BAL.phy_genus, function(OTU) OTU/sum(OTU) * 100)

# Extract the top 20 taxa 
BAL.phy_genus_top20 <- prune_taxa(top20_genera, BAL.phy_genus_relabun)

set.seed(111)
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
cols = sample(col_vector, 20)
cols = cols[-3]
cols = c("#E41A1C", cols[c(1:2)], cols[c(3:19)])
cols[8] = "#E6AB02"
cols[9] = "#A6CEE3"
cols[19] = "#620001"
cols = cols[c(8, 2:7, 1, 9:20)]

orderID = sort(BAL.meta$Sample_ID)[c(1:22,83:90,23:82)]

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
BAL.phy_genus_relabun_BA2_OTU = as.matrix(BAL.phy_genus_relabun_BA2@otu_table)

genera_summary = BAL.phy.df %>% group_by(Group, Genusorder) %>% 
  summarise(mean = mean(Abundance))
genera_summary1 = BAL.phy.df %>% group_by(Group, Genusorder) %>% 
  summarise(sd = sd(Abundance))
BAL.phy.df %>% group_by(Genusorder) %>% 
  summarise(mean = mean(Abundance))
BAL.phy.df %>% group_by(Genusorder) %>% 
  summarise(sd = sd(Abundance))
BAL.phy.df %>% group_by(genusorder) %>% summarise(mean = mean(Abundance))

std.error <- function(x) sd(x)/sqrt(length(x))
BAL.phy.df %>% group_by(genusorder) %>% 
  summarise(sd = sd(Abundance))
BAL.phy.df %>% group_by(genusorder) %>% 
  summarise(sem = std.error(Abundance))

# png("BAL_genus_relative_abundance.png", width = 20, height = 10, units = "in", res = 300)
BAL.phy.df %>% 
  ggplot(aes(x =name, y = Abundance, fill = Genusorder)) +
  geom_bar(stat = "identity") +
  geom_col(position = position_stack(reverse = TRUE)) + 
  labs(x = "",
       y = "Relative Abundance",
       title = "") + #Genus Relative Abundance
  scale_fill_manual(values = c(cols, "gray57")) + 
  facet_grid(Country ~ Visit, scales = "free") +
  labs(fill = "Genus") + 
  theme_classic() + 
  theme(axis.text.x = element_text(size = 7, angle = 90, vjust = 0.5, hjust = 1),
        axis.text.y = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.key.height = unit(.5, 'cm'),
        legend.key.width = unit(.5, 'cm'),
        strip.text = element_text(size = 14),
        # strip.background = element_blank(),
        axis.title.y = element_text(size = 16),
        plot.title = element_text(size = 16)) 
# dev.off()

## 1). we plot BAL samples 
## ii. phylum level relative abundance 
BAL.phy_phylum <- tax_glom(BAL.phy, "Phylum", NArm = TRUE)
# Get top 10 genera
orderPhylum = names(sort(taxa_sums(BAL.phy_phylum), decreasing=TRUE))
top10_phylum <- orderPhylum [1:10]
tax[match(top10_phylum, rownames(tax)),"Phylum"]
# Phylum                        
# otu002 "Firmicutes"                  
# otu005 "Proteobacteria"              
# otu018 "Bacteroidota"                
# otu026 "Actinobacteriota"            
# otu031 "Fusobacteriota" 
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
phyla_summary = BAL.phy.df1 %>% group_by(Group, phylumorder) %>% 
  summarise(mean = mean(Abundance))
phyla_summary1 = BAL.phy.df1 %>% group_by(Group, phylumorder) %>% 
  summarise(sd = sd(Abundance))
BAL.phy.df1 %>% group_by(phylumorder) %>% summarise(mean = mean(Abundance))

std.error <- function(x) sd(x)/sqrt(length(x))
BAL.phy.df1 %>% group_by(phylumorder) %>% 
  summarise(sd = sd(Abundance))
BAL.phy.df1 %>% group_by(phylumorder) %>% 
  summarise(sem = std.error(Abundance))

# png("BAL_phylum_relative_abundance.png", width = 20, height = 10, units = "in", res = 300)
BAL.phy.df1 %>% 
  ggplot(aes(x =name, y = Abundance, fill = phylumorder)) +
  geom_bar(stat = "identity") +
  geom_col(position = position_stack(reverse = TRUE)) + 
  labs(x = "",
       y = "Relative Abundance",
       title = "Phylum Relative Abundance") +
  scale_fill_manual(values = c(cols[c(1:7,10:12)], "gray57")) + 
  facet_grid(Country ~ Visit, scales = "free") +
  labs(fill = "Phylum") + 
  theme_classic() + 
  theme(axis.text.x = element_text(size = 7, angle = 90, vjust = 0.5, hjust = 1),
        axis.text.y = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.key.height = unit(.5, 'cm'),
        legend.key.width = unit(.5, 'cm'),
        strip.text = element_text(size = 12)) 
# dev.off()


## 2). we plot prewash samples 
cols_genus = data.frame(Genus = levels(factor(BAL.phy.df$Genusorder)), Color = c(cols, "gray57"))
cols_phylum = data.frame(Genus = levels(factor(BAL.phy.df1$phylumorder)), Color = c(cols[c(1:7,10:12)], "gray57"))
save(prewash_ID, cols_genus, cols_phylum, file = "color.RData")

## 3). we plot qPCR








# 
# # load("../clinical.RData")
# # colnames(CT1)[1] = "Study_ID"
# # df = CT1 %>% inner_join(meta, by = "Study_ID")%>%
# #   dplyr::select(-ends_with(".y"))
# # keep0 = df$Study_ID
# # BAL.phy0 = subset_samples(BAL.phy, Study_ID %in% keep0)
# # saveRDS(BAL.phy0, "./matched_CT_BAL_T2.rds")
# 
# # temp = paste0(df$name, "BA2")
# # temp = meta[meta$Study_ID %in% temp,]
# # colnames(df)[6] = "Subject_ID" 
# # df= df %>% inner_join(temp, "Subject_ID")%>%
# #   dplyr::select(-ends_with(".y"))
# # idx = which(grepl(".x", colnames(df), fixed))
# # colnames(df)[idx] = str_remove(colnames(df)[idx], ".x")
# # remove = colnames(df) %in% c("Study_ID" )
# # keep = df$Subject_ID  #39
# # # matched sample ID at T1 and T2 that are having matching CT data 
# # BAL.phy = subset_samples(BAL.phy, Subject_ID %in% keep)
# # saveRDS(BAL.phy, "./matched_CT_BAL_2timepoints.rds")
# 
# BAL.phy = readRDS("matched_CT_BAL_2timepoints.rds")
# BAL.phy.T1 = subset_samples(BAL.phy, Visit == "BA2")
# BAL.phy.T2 = subset_samples(BAL.phy, Visit == "BA7")
# otu.T1 = otu_table(BAL.phy.T1)
# otu.T2 = otu_table(BAL.phy.T2)
# 
# ## match CT data and the BAL data 
# ## only keep the matched data 
# load("../clinical.RData")
# load("../pft.RData")
# meta = sample_data(BAL.phy)
# colnames(meta)[1] = "Study_ID"
# colnames(CT1)[1] = "Study_ID"
# co_occ = CT1 %>% left_join(meta, by = "Study_ID")
# 
# colnames(shannon.BAL)[1] = "Study_ID"
# colnames(demo)[which(colnames(demo) == "SampleID")] = "Study_ID"
# 
# ######################## Part 0: demo ######################## 
# # source("./demo.R")
# 
# 
# ######################## Part 1: visualization of relative abundance ########################
# ## i. genus level 
# BAL.phy_genus <- tax_glom(BAL.phy, "Genus", NArm = TRUE)
# # Get top 20 genera
# prewash.orderTaxa = names(sort(taxa_sums(BAL.phy_genus), decreasing=TRUE))
# 
# # apply(BAL.phy_genus_relabun@otu_table, 2, summary)
# # sum(rowSums(BAL.phy@otu_table)/sum(rowSums(BAL.phy@otu_table)) > 0.02)
# 
# top20_genera <- prewash.orderTaxa[1:20]
# # Transform Taxa counts to relative abundance
# BAL.phy_genus_relabun <- transform_sample_counts(BAL.phy_genus, function(OTU) OTU/sum(OTU) * 100)
# 
# # Extract the top 20 taxa 
# BAL.phy_genus_top20 <- prune_taxa(top20_genera, BAL.phy_genus_relabun)
# 
# set.seed(111)
# qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
# col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
# cols = sample(col_vector, 20)
# cols = cols[-3]
# cols = c("#E41A1C", cols[c(1:2)], cols[c(3:19)])
# cols[8] = "#E6AB02"
# cols[9] = "#A6CEE3"
# cols[19] = "#620001"
# cols = cols[c(8, 2:7, 1, 9:20)]
# 
# orderID = sort(meta$Study_ID)[c(1:18,73:78,19:72)]
# 
# #### showing the top 20 taxa separately, and the rest in gray 
# BAL.phy.df = psmelt(BAL.phy_genus_relabun)  %>% 
#   mutate(name = fct_relevel(Sample, orderID)) 
# BAL.phy.df$Genus[which(BAL.phy.df$OTU %in% top20_genera == FALSE)] = "Low Abundance"
# top20_genera_name = BAL.phy.df$Genus[match(top20_genera, BAL.phy.df$OTU)]
# BAL.phy.df = BAL.phy.df %>% mutate(Genusorder = fct_relevel(Genus, c(top20_genera_name, "Low Abundance")))
# BAL.phy.df$Country = factor(BAL.phy.df$Country, levels = c("USA", "AUS"))
# BAL.phy.df$Group = paste0(BAL.phy.df$Country, " ", BAL.phy.df$Visit)
# BAL.phy.df$Group = factor(BAL.phy.df$Group, levels = c("USA BA2","AUS BA2",  
#                                                        "USA BA7", "AUS BA7"))
# BAL.phy_genus_relabun_BA2 = subset_samples(BAL.phy_genus_relabun, Visit == "BA2")
# BAL.phy_genus_relabun_BA7 = subset_samples(BAL.phy_genus_relabun, Visit == "BA7")
# BAL.phy_genus_relabun_BA2_OTU = as.matrix(BAL.phy_genus_relabun_BA2@otu_table)
# 
# 
# # try to keep the genus that are present in more than 15% of the samples 
# otu = otu_table(BAL.phy_genus_relabun)
# otu.df = as.data.frame(otu)
# otu.df = as.data.frame(t(otu.df))
# otu.df$Study_ID = rownames(otu.df)
# rownames(otu.df) = NULL
# 
# dim(otu) #352 78
# otu = t(otu)
# otu = as.data.frame(otu)
# zero.df = as.data.frame(lapply(otu, function(x){ sum(x==0)}))
# zero.df = as.data.frame(t(zero.df))
# nonzero.df = data.frame(nonzero = dim(otu)[1]-zero.df$V1)
# rownames(nonzero.df) = rownames(zero.df)
# rownames(nonzero.df) = tax_table(BAL.phy)[match(rownames(nonzero.df), 
#                          rownames(otu_table(BAL.phy))), "Genus"]
# summary(nonzero.df$nonzero)
# hist(nonzero.df$nonzero)
# # FOR ANCOM: keep the genera that were prevalent in more than 10% of the samples 
# keep = rownames(nonzero.df)[which(nonzero.df$nonzero > dim(otu)[1] * 0.1)]
# keep 
# 
# ## clustering on BA2 samples 
# rownames(BAL.phy_genus_relabun_BA2_OTU) = 
#   BAL.phy_genus_relabun@tax_table[match(rownames(BAL.phy_genus_relabun_BA2_OTU), 
#                                         rownames(BAL.phy_genus_relabun@tax_table)), "Genus"]
# tempname = substring(colnames(BAL.phy_genus_relabun_BA2_OTU),1,1)
# tempcountry = ifelse(tempname %in% c("M", "P"), "AUS", "USA")
# my_sample_col = data.frame(Country = tempcountry)
# rownames(my_sample_col) = colnames(BAL.phy_genus_relabun_BA2_OTU)
# 
# 
# 
# 
# 
# # cluster the samples based on the microbiome at BA2
# # png("../Results/cluster_BA2_3clusters.png", width = 10, height = 20,
#     # units = "in", res = 300)
# # heatmap = pheatmap(BAL.phy_genus_relabun_BA2_OTU, cutree_cols = 3,
# #                    annotation_col = my_sample_col,
# #                    fontsize = 4)
# # dev.off()
# 
# ## clustering on BA2+BA7 samples 
# BAL.phy_genus_relabun_OTU = as.matrix(BAL.phy_genus_relabun@otu_table)
# 
# rownames(BAL.phy_genus_relabun_BA2_OTU) = 
#   BAL.phy_genus_relabun@tax_table[match(rownames(BAL.phy_genus_relabun_BA2_OTU), 
#                                         rownames(BAL.phy_genus_relabun@tax_table)), "Genus"]
# tempname = substring(colnames(BAL.phy_genus_relabun_BA2_OTU),1,1)
# tempcountry = ifelse(tempname %in% c("M", "P"), "AUS", "USA")
# my_sample_col = data.frame(Country = tempcountry)
# rownames(my_sample_col) = colnames(BAL.phy_genus_relabun_BA2_OTU)
# heatmap = pheatmap(BAL.phy_genus_relabun_BA2_OTU, cutree_cols = 5,
#                    annotation_col = my_sample_col,
#                    fontsize = 4)
# 
# 
# # cluster = cutree(heatmap$tree_col, k = 5)
# # BAL.phy_genus_relabun_BA2_OTU = rbind(cluster, BAL.phy_genus_relabun_BA2_OTU)
# 
# # png("./Results/genus_relative_abundance.png", width = 20, height = 10, units = "in", res = 300)
# BAL.phy.df %>% 
#   ggplot(aes(x =name, y = Abundance, fill = Genusorder)) +
#   geom_bar(stat = "identity") +
#   geom_col(position = position_stack(reverse = TRUE)) + 
#   labs(x = "",
#        y = "Relative Abundance",
#        title = "Genus Relative Abundance") +
#   scale_fill_manual(values = c(cols, "gray57")) + 
#   facet_grid(Country ~ Visit, scales = "free") +
#   labs(fill = "Genus") + 
#   theme_classic() + 
#   theme(axis.text.x = element_text(size = 7, angle = 90, vjust = 0.5, hjust = 1),
#     axis.text.y = element_text(size = 12),
#     legend.text = element_text(size = 14),
#     strip.text = element_text(size = 14),
#     # strip.background = element_blank(),
#     axis.title.y = element_text(size = 16),
#     plot.title = element_text(size = 16)) 
# # dev.off()
# 
# 
# ## ii. phylum level relative abundance 
# BAL.phy_phylum <- tax_glom(BAL.phy, "Phylum", NArm = TRUE)
# # Get top 10 genera
# orderPhylum = names(sort(taxa_sums(BAL.phy_phylum), decreasing=TRUE))
# top10_phylum <- orderPhylum [1:10]
# # Transform Taxa counts to relative abundance
# BAL.phy_phylum_relabun <- transform_sample_counts(BAL.phy_phylum, function(OTU) OTU/sum(OTU) * 100)
# 
# # Extract the top 10 taxa 
# BAL.phy_phylum_top10 <- prune_taxa(top10_phylum, BAL.phy_phylum_relabun)
# 
# #### showing the top 10 taxa separately, and the rest in gray 
# BAL.phy.df1 = psmelt(BAL.phy_phylum_relabun)  %>% 
#   mutate(name = fct_relevel(Sample, orderID)) 
# BAL.phy.df1$Phylum[which(BAL.phy.df1$OTU %in% top10_phylum == FALSE)] = "Low Abundance"
# top10_phylum_name = BAL.phy.df1$Phylum[match(top10_phylum, BAL.phy.df1$OTU)]
# BAL.phy.df1 = BAL.phy.df1 %>% mutate(phylumorder = fct_relevel(Phylum, c(top10_phylum_name, "Low Abundance")))
# BAL.phy.df1$Country = factor(BAL.phy.df1$Country, levels = c("USA", "AUS"))
# BAL.phy.df1$Group = paste0(BAL.phy.df1$Country, " ", BAL.phy.df1$Visit)
# BAL.phy.df1$Group = factor(BAL.phy.df1$Group, levels = c("USA BA2","AUS BA2",  
#                                                        "USA BA7", "AUS BA7"))
# 
# # png("./Results/phylum_relative_abundance.png", width = 20, height = 10, units = "in", res = 300)
# BAL.phy.df1 %>% 
#   ggplot(aes(x =name, y = Abundance, fill = phylumorder)) +
#   geom_bar(stat = "identity") +
#   geom_col(position = position_stack(reverse = TRUE)) + 
#   labs(x = "",
#        y = "Relative Abundance",
#        title = "Phylum Relative Abundance") +
#   scale_fill_manual(values = c(cols[c(1:7,10:12)], "gray57")) + 
#   facet_grid(Country ~ Visit, scales = "free") +
#   labs(fill = "Phylum") + 
#   theme_classic() + 
#   theme(axis.text.x = element_text(size = 7, angle = 90, vjust = 0.5, hjust = 1),
#         axis.text.y = element_text(size = 12),
#         legend.text = element_text(size = 10),
#         strip.text = element_text(size = 12)) 
# # dev.off()
# 
# 
# ## iii. species level 
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
#   mutate(name = fct_relevel(Sample, orderID))
# BAL.phy.df2$Species = paste0(BAL.phy.df2$Genus, " ", BAL.phy.df2$Species)
# BAL.phy.df2$Species[which(BAL.phy.df2$OTU %in% top20_species == FALSE)] = "Low Abundance"
# top20_species_name = BAL.phy.df2$Species[match(top20_species, BAL.phy.df2$OTU)]
# BAL.phy.df2 = BAL.phy.df2 %>% mutate(Speciesorder = fct_relevel(Species, c(top20_species_name, "Low Abundance")))
# BAL.phy.df2$Country = factor(BAL.phy.df2$Country, levels = c("USA", "AUS"))
# BAL.phy.df2$Group = paste0(BAL.phy.df2$Country, " ", BAL.phy.df2$Visit)
# BAL.phy.df2$Group = factor(BAL.phy.df2$Group, levels = c("USA BA2","AUS BA2",  
#                                                        "USA BA7", "AUS BA7"))
# 
# # png("./Results/species_relative_abundance.png", width = 20, height = 10, units = "in", res = 300)
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
# 
# 
# 

# 
# ### investigate whether Staphylococcus is DA; not present 
# # staphy = BAL.phy_genus_abs.df %>% filter(Genus == "Staphylococcus")
# # kruskal.test(Abundance ~ Group, data = staphy)
# 
# Haemophilus = BAL.phy_genus_abs.df %>% filter(Genus == "Haemophilus")
# ggplot(Haemophilus, aes(x =name, y = Abundance, fill = Country)) +
#   geom_bar(stat = "identity") +
#   geom_col(position = position_stack(reverse = TRUE)) + 
#   labs(x = "",
#        y = expression("copies/g sample"),
#        title = "") +
#   scale_fill_manual(values = c(cols, "gray57")) +
#   facet_grid(~ Visit, scales = "free") +
#   labs(fill = "Genus") + 
#   theme_classic() + 
#   theme(axis.text.x = element_text(size = 7, angle = 90, vjust = 0.5, hjust = 1),
#         axis.text.y = element_text(size = 12),
#         legend.text = element_text(size = 14),
#         strip.text = element_text(size = 14),
#         axis.title.y = element_text(size = 16),
#         plot.title = element_text(size = 16)) 
# 
# pvals = rep(NA, nlevels(factor(BAL.phy_genus_abs.df$Genus)))
# for (i in 1:nlevels(factor(BAL.phy_genus_abs.df$Genus))){
#   genus = levels(factor(BAL.phy_genus_abs.df$Genus))[i]
#   df = BAL.phy_genus_abs.df %>% filter(Genus == genus)
#   pvals[i] = kruskal.test(Abundance ~ Group, data = df)$p.value
# }
# names(pvals) = levels(factor(BAL.phy_genus_abs.df$Genus))
# padj = p.adjust(pvals, method = "hochberg", n = nlevels(factor(BAL.phy_genus_abs.df$Genus)))
# which(padj < 0.05)
# 
# ########################alpha diversity:  shannon diversity ######################## 
# shannon.BAL = read.csv("../../Data/shannon_BAL_decontam.csv")
# colnames(shannon.BAL)[1] = "Study_ID"
# shannon.BAL = shannon.BAL %>% inner_join(meta, by = "Study_ID")
# 
# ## shannon by country
# test.aov <- aov(Shannon ~ Visit*Country, shannon.BAL)
# summary(test.aov)
# # Df Sum Sq Mean Sq F value  Pr(>F)   
# # Visit          1   2.18   2.178   4.942 0.02927 * 
# #   Country        1   3.70   3.699   8.393 0.00495 **
# #   Visit:Country  1   0.73   0.729   1.654 0.20238   
# # Residuals     74  32.61   0.441                         
# # ---
# #   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# aov_residuals <- residuals(object = test.aov)
# shapiro.test(x = aov_residuals) # 0.086
# TukeyHSD(test.aov, "Visit:Country")
# # $`Visit:Country`
# # diff         lwr         upr     p adj
# # BA7:USA-BA2:USA  0.04414646 -0.66819611  0.75648904 0.9984481
# # BA2:AUS-BA2:USA -0.68130283 -1.28667761 -0.07592805 0.0211720 *
# # BA7:AUS-BA2:USA -0.21820983 -0.82358461  0.38716496 0.7793874
# # BA2:AUS-BA7:USA -0.72544929 -1.33082407 -0.12007451 0.0123540 *
# # BA7:AUS-BA7:USA -0.26235629 -0.86773107  0.34301849 0.6666990
# # BA7:AUS-BA2:AUS  0.46309300 -0.01180205  0.93798805 0.0586489 .
# 
# # png("./Results/shannon_country.png", height = 10, width = 10, units = "in", res = 300)
# ggplot(shannon.BAL, aes(x = Group, y =  Shannon, fill = Group)) + 
#   geom_boxplot() +
#   ylab("Shannon Diversity") +
#   xlab("") + 
#   theme_classic() +
#   scale_x_discrete(breaks=c("USABA2", "AUSBA2", "USABA7", "AUSBA7"),
#                    labels=c("USA BA2", "AUS BA2", "USA BA7", "AUS BA7"))+
#   theme(legend.title = element_text(size = 16),
#         axis.title = element_text(size = 18, face = "bold"),
#         axis.text.x = element_text(size = 16, face = "bold"),
#         axis.text.y = element_text(size = 16)) +
#   # annotate("text", x=1.5,y=4.2, label="BA2",size =7) +
#   # annotate("text", x=3.5,y=4.2, label="BA7",size =7) +
#   annotate("segment", x=1,xend=2,y=3.4,yend=3.4) +
#   annotate("segment", x=2,xend=3,y=3.6,yend=3.6) +
#   annotate("segment", x=2,xend=4,y=3.8,yend=3.8) +
#   annotate("text", x=1.5,y=3.45, label="*",size=7) +
#   annotate("text", x=2.5,y=3.65, label="*",size=7) +
#   annotate("text", x=3,y=4, label=".",size=10) 
# # dev.off()
# 
# ########################beta diversity ######################## 
# # source("./beta_diversity.R")
# 
# 
# ########################
# # load meta data 
# # load demographic data 
# 
# ######################## Part : CT ########################
# ## CT was only performed at BA7
# # CT was performed closer to BA7 
# CT1 = CT1 %>% inner_join(meta, by = "Study_ID") 
# table(CT1$Country.x)
# # AUS USA 
# # 27  12 
# ### Bronchiectasis
# ### 4/12 in USA; 5/27 in AUS
# ## Mucus Plugging 
# ### 1/15 in USA; 5/50 in AUS
# ## Abnormal 
# CT1 = CT1[,-which(grepl(".y", colnames(CT1), fixed = T))]
# colnames(CT1)[which(grepl(".x", colnames(CT1), fixed = T))] = 
#   str_remove(colnames(CT1)[which(grepl(".x", colnames(CT1), fixed = T))],
#              ".x")
# CT1 = CT1[,-which(colnames(CT1) %in% c("is.neg", "Sample_or_Control",
#                            "Visit.1", "name"))]
# 
# # plot Bronchiectasis, Mucus Plugging and Abnormal separately
# CT1.long = melt(CT1 %>% dplyr::select(c("Study_ID", "Subject_ID", "Country", "Bronchiectasis", "Mucus Plugging", "Abnormal")), 
#                 id.vars = c("Study_ID", "Subject_ID", "Country"), variable.name = "CT")
# CT1.long$Country = factor(CT1.long$Country, levels = c("USA", "AUS"))
# orderID1 = substring(orderID, 1, nchar(orderID)-3)
# orderID1 = orderID1[-which(duplicated(orderID1))]
# CT1.long = CT1.long %>% mutate(name = fct_relevel(`Subject_ID`, orderID1)) 
# # png("./Results/CT.png", width = 20, height = 10, units = "in", res = 300)
# ggplot(CT1.long, aes(x = name, y = value, fill = CT)) +
#   geom_bar(position="stack", stat="identity") +
#   ylab("Diseased Percentage") +
#   xlab("") + 
#   theme_classic() + 
#   geom_vline(xintercept = 13.5, linetype = "dashed") +
#   annotate("text", x = 6, y = 5, label = "USA", size = 6) + 
#   annotate("text", x = 25.5, y = 5, label = "AUS", size = 6) +
#   theme(axis.text.x = element_text(angle = 90),
#         legend.title = element_text(size = 16),
#         axis.title = element_text(size = 18, face = "bold"),
#         axis.text.y = element_text(size = 16)) 
# # dev.off()
# 
# 
# # png("./Results/CT1.png", width = 20, height = 10, units = "in", res = 300)
# ggplot(CT1.long, aes(x = Country, y = value, fill = CT)) + 
#   geom_boxplot() +
#   ylab("Percentage") +
#   xlab("") + 
#   theme_classic() + 
#   theme(axis.text.x = element_text(size = 18),
#         legend.title = element_text(size = 16),
#         axis.title = element_text(size = 18, face = "bold"),
#         axis.text.y = element_text(size = 16)) 
# # dev.off()
# 
# # plot diseased% separately
# CT1.long.disease = melt(CT1 %>% dplyr::select(c("Study_ID", "Subject_ID", "Country", "Disease")), 
#                          id.vars = c("Study_ID", "Subject_ID", "Country"), variable.name = "CT")
# CT1.long.disease$Country = factor(CT1.long.disease$Country, levels = c("USA", "AUS"))
# 
# disease.USA = CT1.long.disease %>% filter(Country == "USA") %>% dplyr::select(value)
# disease.AUS = CT1.long.disease %>% filter(Country == "AUS") %>% dplyr::select(value)
# 
# wilcox.test(disease.USA$value, disease.AUS$value,exact = F)
# # p-value = 0.02737
# 
# ## CT by country
# # png("./Results/CT_disease_country.png", width = 10, height = 10, units = "in", res = 300)
# ggplot(CT1.long.disease, aes(x = Country, y = value, fill = Country)) + 
#   geom_violin(width=.7) +
#   geom_boxplot(width=0.2, color="black", alpha=0.6) +
#   ylab("Diseased Percentage") +
#   xlab("") + 
#   theme_classic() +
#   theme(legend.position = "none",
#         axis.title = element_text(size = 18, face = "bold"),
#         axis.text.x = element_text(size = 16, face = "bold"),
#         axis.text.y = element_text(size = 16)) +
#   annotate("text", x= 2, y = 4, label = "p-value = 0.027", size = 7)
# # dev.off()
# 
# 
# ## Cytokines were measured at both BA2 and BA7 
# cytokine = read_excel("../../../Clinical Data/cytokine_cleaned_WW.xlsx")
# cytokine1 = cytokine %>% mutate("Study_ID" = paste0(cytokine$SampleID,
#                                           cytokine$Visit)) %>% 
#   inner_join(meta, by = "Study_ID")
# 
# cytokine1 = cytokine1[,-which(grepl(".y", colnames(cytokine1), fixed = T))]
# colnames(cytokine1)[which(grepl(".x", colnames(cytokine1), fixed = T))] = 
#   str_remove(colnames(cytokine1)[which(grepl(".x", colnames(cytokine1), fixed = T))],
#              ".x")
# cytokine1 = cytokine1[,-which(colnames(cytokine1) %in% c("is.neg", "Sample_or_Control",
#                                        "Visit.1", "name"))]
# 
# 
# cytokine.long =  melt(cytokine1 %>% dplyr::select(c("SampleID", "Visit", 
#                                                    "Site","Country","IL8c_pgml")),
#                       id.vars = c("SampleID", "Visit", "Site","Country"),
#                       variable.name = "IL8c_pgml")
# cytokine.long$Group = paste0(cytokine.long$Country, cytokine.long$Visit)
# cytokine.long$Group = factor(cytokine.long$Group, levels = c("USABA2", "AUSBA2",
#                                                              "USABA7", "AUSBA7"))
# 
# ggplot(cytokine.long, aes(x = Group, y = value, fill = Country)) + 
#   geom_boxplot() +
#   ylab("IL8 (pg/mL)") +
#   xlab("") + 
#   theme_classic() + 
#   scale_y_continuous(trans='log10') +
#   scale_x_discrete(breaks=c("USABA2", "AUSBA2", "USABA7", "AUSBA7"),
#                   labels=c("USA \nBA2", "AUS \nBA2", "USA \nBA7", "AUS \nBA7")) + 
#   theme(axis.title = element_text(size = 18, face = "bold"),
#        axis.text.x = element_text(size = 16, face = "bold"),
#        axis.text.y = element_text(size = 16))
# 
# cytokine.long = cytokine.long[,-which(colnames(cytokine.long) == "IL8c_pgml")]
# colnames(cytokine.long)[which(colnames(cytokine.long) == "value")] = "IL8"
# 
# test.aov1 <- aov(IL8 ~ Visit*Country, cytokine.long)
# aov_residuals1 <- residuals(object = test.aov1)
# shapiro.test(x = aov_residuals1) 
# # 2.665e-09 violates assumption
# 
# kruskal.test(IL8 ~ Group, data = cytokine.long)
# kruskal.test(log10(IL8) ~ Group, data = cytokine.long) # not significant 
# 
# 
# ## check the correlation between IL8 at T1 and T2 
# cytokine.wide = reshape(cytokine.long, idvar = "SampleID", timevar = "Visit", direction = "wide")
# cytokine.wide = cytokine.wide[, -which(colnames(cytokine.wide) %in% 
#                                    c("Site.BA7","Country.BA7","Group.BA7"))]
# cytokine.wide = cytokine.wide[,c(1,2,3,5,4,6)]
# colnames(cytokine.wide)[2:4] = c("Site","Country", "Group")
# 
# # png("./Results/IL8_cor_T2_T1.png", width = 8, height = 6, units = "in", res=300)
# ggplot(cytokine.wide, aes(x = log2(IL8.BA2), y = log2(IL8.BA7), color = Country)) + 
#   geom_point() + 
#   ylab("log(IL8) T2") +
#   xlab("log(IL8) T1") +
#   stat_smooth(method = "lm", col = "red") +
#   theme_classic() +
#   annotate("text", x = 16, y = 16.5, label = "p = 0.0033", size = 4)
# # dev.off()
# fit_IL8 = lm(log2(IL8.BA7) ~ log2(IL8.BA2), cytokine.wide)
# summary(fit_IL8) # 0.003326 ** 
# 
# ## check the correlation between IL8 and CT 
# colnames(cytokine.wide)[1] = "Subject_ID"
# IL8CT = cytokine.wide %>% inner_join(CT1.long.disease, by = "Subject_ID")
# colnames(IL8CT)[10] = "Disease"
# IL8CT = IL8CT[,-9]
# ggplot(IL8CT, aes(x = log2(IL8.BA2), y = Disease)) +
#   geom_point() +
#   stat_smooth(method = "lm", col = "red")+
#   theme_classic() +
#   annotate("text", x = 15.5, y = 3.5, label = "p = 0.475", size = 4)
# fit1_IL8 = lm(Disease ~ log2(IL8.BA2), IL8CT)
# summary(fit1_IL8)
# 
# ggplot(IL8CT, aes(x = log2(IL8.BA7), y = Disease)) +
#   geom_point() +
#   stat_smooth(method = "lm", col = "red")+
#   theme_classic() +
#   annotate("text", x = 15.5, y = 3.5, label = "p = 0.578", size = 4)
# fit2_IL8 = lm(Disease ~ log2(IL8.BA7), IL8CT)
# summary(fit2_IL8)
# 
# quantile(CT1$Disease, seq(0,1,.01)) # 1.22663551 50%
# 
# IL8CT = IL8CT %>% mutate(Quant = ifelse(IL8CT$Disease < 1.22663551,0,1))
# # test whether CT high and CT low are disproportional in AUS vs. USA 
# # check the relative abundance of the samples in Q1 and Q2 at time point 1
# Q1 = IL8CT$Subject_ID[which(IL8CT$Quant == 0)]
# Q2 = IL8CT$Subject_ID[which(IL8CT$Quant == 1)]
# 
# table(IL8CT$`Country.x`[which(IL8CT$Quant == 0)])  
# # AUS USA 
# # 17   2 
# table(IL8CT$`Country.x`[which(IL8CT$Quant == 1)])  
# # AUS USA 
# # 10  10
# fisher.df = data.frame("AUS" = c(10,17), 
#                        "USA" = c(10, 2),
#                        row.names = c("high", "low"))
# mosaicplot(fisher.df, color = T)
# fisher.test(fisher.df) # 0.01381
# 
# IL8CT$Quant = factor(IL8CT$Quant, levels = c(0,1))
# ggplot(IL8CT, aes(y = log2(IL8.BA2), x = Quant)) +
#   geom_boxplot() + 
#   theme_classic()
# 
# 
# ggplot(IL8CT, aes(y = log2(IL8.BA7), x = Quant)) +
#   geom_boxplot()+ 
#   theme_classic()
# 
# wilcox.test(log2(IL8CT$IL8.BA2[IL8CT$Quant==0]),
#             log2(IL8CT$IL8.BA2[IL8CT$Quant==1]), paired = F)
# 
# wilcox.test(log2(IL8CT$IL8.BA7[IL8CT$Quant==0]),
#             log2(IL8CT$IL8.BA7[IL8CT$Quant==1]), paired = F)
# 
# ######################## Part : FEV ########################
# # source("./pft.R")
# 
# 
# 
# 
# 
# ############################## Part: significant abundant using DESeq2##############################
# ############################## trying AMCOM-BC for DA analysis, see ANCOMBC.R ############################## 
# # setting USA as reference level
# # setting BA2 as reference level
# ###### 1. univariate; without adjusting for any covariates 
# # quantile(rowSums(BAL.phy@otu_table), seq(0,1,0.01))
# # BAL.phy@sam_data$Group = paste0(BAL.phy@sam_data$Country, 
# #                                 BAL.phy@sam_data$Visit)
# # BAL.phy@sam_data$Group = factor(BAL.phy@sam_data$Group, 
# #                                 levels = c("USABA2", "AUSBA2",  
# #                                            "USABA7", "AUSBA7"))
# # BAL.phy@sam_data$Country = factor(BAL.phy@sam_data$Country, 
# #                                   levels = c("USA", "AUS"))
# # BAL.phy@sam_data$Visit = factor(BAL.phy@sam_data$Visit, 
# #                                   levels = c("BA2", "BA7"))
# # 
# # demo$Male = ifelse(demo$Gender == "Male", 1, 0)
# # demo$Male = factor(demo$Male, levels = c(0,1))
# # 
# # ## "White or Caucasian"=0
# # ## "Native Hawaiian or Other Pacific Islander"=1                                               
# # ## "More than one race"=2
# # demo$Race1 = ifelse(demo$Race == "More than one race", 2,
# #                     ifelse(demo$Race == "White or Caucasian (origins or descendent of people in Europe, North Africa, or Middle East)", 0, 1))
# # demo$Race1 = factor(demo$Race1, levels = c(0,1,2))
# # 
# # ## heterogenous = 0
# # ## homozygous = 1
# # ## other = 2
# # demo$homozygous = ifelse(demo$Genotype_F508 == "Heterogenous", 0,
# #                          ifelse(demo$Genotype_F508 == "Homozygous", 1, 2))
# # demo$homozygous = factor(demo$homozygous, levels = c(0,1,2))
# # 
# # demo$BMI = demo$Weight/((demo$Height/100)^2)
# # 
# # sample_data(BAL.phy) = cbind.data.frame(as.data.frame(BAL.phy@sam_data), 
# #                                         demo[match(rownames(BAL.phy@sam_data), demo$SampleID),])
# # 
# # # saveRDS(BAL.phy, "BAL_phy.rds") 
# # 
# # # save(demo, shannon.BAL, CT1, cytokine, file = "clinical.RData")
# # 
# 
# 
