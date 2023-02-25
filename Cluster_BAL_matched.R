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

dirname = dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(dirname)
BAL.phy = readRDS("../../Data/phy_noncontam.rds")
# only keep the samples
BAL.phy = subset_samples(BAL.phy, Sample_type == "BAL")
meta = BAL.phy@sam_data
meta$Visit[which(meta$Visit == "")] = "IASM"
meta$Sequence_ID = rownames(meta)
## Remove IASM samples
BAL.phy = subset_samples(BAL.phy, Visit %in% c("BA2", "BA7"))
meta = BAL.phy@sam_data

BAL.BA2.phy = subset_samples(BAL.phy, Visit == "BA2")
BAL.BA7.phy = subset_samples(BAL.phy, Visit == "BA7")
BA2.subID = BAL.BA2.phy@sam_data$Subject_ID
BA7.subID = BAL.BA7.phy@sam_data$Subject_ID
length(intersect(BA2.subID, BA7.subID)) #45
intersect(BA2.subID, BA7.subID)
sum(BAL.phy@sam_data$Visit == "BA2" & BAL.phy@sam_data$Country == "USA") #19
sum(BAL.phy@sam_data$Visit == "BA2" & BAL.phy@sam_data$Country == "AUS") #32
sum(BAL.phy@sam_data$Visit == "BA7" & BAL.phy@sam_data$Country == "USA") #15
sum(BAL.phy@sam_data$Visit == "BA7" & BAL.phy@sam_data$Country == "AUS") #43
sum(BAL.phy@sam_data$Visit == "BA7") #58
sum(BAL.phy@sam_data$Visit == "BA2") #51

meta1 = sample_data(BAL.phy)
USABA2.ID = meta1$Subject_ID[meta1$Visit == "BA2" & meta1$Country == "USA"] # 19
USABA7.ID = meta1$Subject_ID[meta1$Visit == "BA7" & meta1$Country == "USA"] # 15
AUSBA2.ID = meta1$Subject_ID[meta1$Visit == "BA2" & meta1$Country == "AUS"] # 32
AUSBA7.ID = meta1$Subject_ID[meta1$Visit == "BA7" & meta1$Country == "AUS"] # 43
intersect(USABA2.ID, USABA7.ID) # 15
intersect(AUSBA2.ID, AUSBA7.ID) # 30

commonID.T1T2 = intersect(BA2.subID, BA7.subID) 

BAL.phy = subset_samples(BAL.phy, Subject_ID %in% commonID.T1T2)
meta = sample_data(BAL.phy)
otu = otu_table(BAL.phy)
tax = tax_table(BAL.phy)
colnames(otu) == rownames(meta)
rownames(meta) = meta$Sample_ID
colnames(otu) = meta$Sample_ID
colnames(otu) == rownames(meta)
BAL.phy = phyloseq(otu, tax, meta)

# ## clustering on BA2+BA7 samples 
meta = sample_data(BAL.phy)
sort(meta$Sample_ID)
orderID = sort(meta$Sample_ID)[c(1:22, 83:90, 23:82)]

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

#### showing the top 20 taxa separately, and the rest in gray 
BAL.phy.df = psmelt(BAL.phy_genus_relabun)  %>% 
  mutate(name = fct_relevel(Sample, orderID)) 
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

# try to keep the genus that are present in more than 15% of the samples 
otu = otu_table(BAL.phy_genus_relabun)
otu.df = as.data.frame(otu)
otu.df = as.data.frame(t(otu.df))
otu.df$`Study ID` = rownames(otu.df)
rownames(otu.df) = NULL

dim(otu) #352 90
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

##### clustering on BA2 + BA7 samples 
### reference https://susan.su.domains/papers/Pregnancy/PNAS_Vaginal_Analysis.Rmd
BAL.phy_genus <- tax_glom(BAL.phy, "Genus", NArm = TRUE)
BAL.phy_genus_relabun <- transform_sample_counts(BAL.phy_genus, function(OTU) OTU/sum(OTU) * 100)

braydist = phyloseq::distance(BAL.phy_genus_relabun, method="bray")
ord = ordinate(BAL.phy_genus_relabun, method = "MDS", distance = braydist)
print(plot_scree(ord)) + theme_classic() +
  ggtitle("MDS-bray ordination eigenvalues")
evs <- ord$value$Eigenvalues
print(evs[1:20])
print(tail(evs))

h_sub5 <- hist(evs[6:length(evs)], 100)
plot(h_sub5$mids, h_sub5$count, log="y", type='h', lwd=10, lend=2)

library(cluster)

NDIM = 46
x = ord$vectors[,1:NDIM]  
pamPCoA = function(x, k) {
  list(cluster = pam(x[,1:2], k, cluster.only = TRUE))
}
gs = clusGap(x, FUN = pamPCoA, K.max = 12, B = 50)
plot_clusgap(gs) + scale_x_continuous(breaks=c(seq(0, 12, 2)))
#  We select 3 clusters based on the gap statistics 

library(gridExtra)
library(grid)
## Perform PAM 3-fold clusters
K = 3
x = ord$vectors[,1:NDIM]
clust = as.factor(pam(x, k=K, cluster.only=T))

sample_data(BAL.phy_genus_relabun)$Cluster <- clust
Clusters <- as.character(seq(K))

ClusterColors = brewer.pal(6,"Paired")[c(1,3,5)] 
names(ClusterColors) = c("Cluster1", "Cluster2", "Cluster3")
CSTFillScale = scale_fill_manual(values = c("1"="#A6CEE3","2"="#B2DF8A","3"="#FB9A99"),
                                 labels = c("1", "2", "3"))
plot_ordination(BAL.phy_genus_relabun, ord, color="Cluster") 
plot_ordination(BAL.phy_genus_relabun, ord, axes=c(3,4), color="Cluster") 
braydist = phyloseq::distance(BAL.phy_genus_relabun, method="bray")

plot_ordination(BAL.phy_genus_relabun, 
                ordinate(BAL.phy_genus_relabun, method="MDS", distance=braydist), 
                color="Cluster") + ggtitle("MDS -- bray -- By Cluster")

taxa.order = names(sort(taxa_sums(BAL.phy_genus_relabun)))
tax = tax_table(BAL.phy_genus_relabun)
taxa.order = tax[match(taxa.order, rownames(tax)),"Genus"]
taxa.order1 = rownames(taxa.order)
pshm = prune_taxa(names(sort(taxa_sums(BAL.phy_genus_relabun), T))[1:25], 
                   BAL.phy_genus_relabun)
for(Cluster in Clusters) {
  Cluster = 1
  pshm1 = prune_samples(sample_data(pshm)$Cluster == Cluster, pshm)
  print(plot_heatmap(pshm1, taxa.label="Genus", taxa.order=taxa.order1) + 
          ggtitle(paste("Cluster:", Cluster)))
}

sample.order = rownames(sample_data(pshm)[order(get_variable(pshm, "Cluster"))])
heatmap.df = as.data.frame(otu_table(pshm), data.frame)
heatmap.df = heatmap.df[,sample.order]
meta1 = sample_data(pshm)
meta1 = as.data.frame(as.matrix(meta1))
my_sample_col = data.frame(Country = meta1$Country[match(sample.order, meta1$Sample_ID)],
                           Visit = meta1$Visit[match(sample.order, meta1$Sample_ID)],
                           Cluster = meta1$Cluster[match(sample.order, meta1$Sample_ID)])
rownames(my_sample_col) = rownames(meta1)
# write.csv(my_sample_col, "my_sample_col_matchedBAL.csv")



