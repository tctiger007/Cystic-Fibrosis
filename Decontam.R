## load library
library(phyloseq)
library(ggplot2)
library(decontam)
library(stringr)

## Set work directory
#setwd("C:/Users/bency/Box/CF2/Run287") # for lab computer
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

## load data
OTU <- read.csv("../Data/otu_BAL_decipher.csv", row.names=1)
Tax <- read.csv("../Data/tax_BAL_decipher.csv", row.names=1)
# --- already run 
# meta <- read.csv("../Data/meta_BAL.csv", row.names=1)
# meta$Sample_ID = str_remove(meta$Sample_ID, "c")
# meta$Visit = lapply(str_split(meta$Sample_ID, "BA"), "[",2)
# meta$Visit[is.na(meta$Visit)] = ""
# meta$Visit = paste0("BA", meta$Visit)
# meta$Visit[meta$Visit == "BA"] = ""
# meta$Visit[meta$Visit == "BA7c"] = "BA7"
# meta$Subject_ID = str_remove(meta$Sample_ID, meta$Visit)
# meta$Site = gsub("[0-9]", "", meta$Subject_ID)
# IASM.idx = which(grepl("^[0-9]+$", meta$Sample_ID))
# for(i in 1:length(IASM.idx)){
#   idx = IASM.idx[i]
#   temp = paste0("IASM", meta$Sample_ID[idx])
#   meta$Subject_ID[idx] = temp
# }
# meta$Site[IASM.idx] = "IU"
# meta$Study = rep(NA, dim(meta)[1])
# meta$Study[IASM.idx] = "IASM"
# meta$Study[is.na(meta$Site)] = ""
# meta$Study[-c(IASM.idx, which(is.na(meta$Site)))] = "ECF"
# meta$Country = ifelse(meta$Site %in% c("M", "P"), "AUS", 
#                       ifelse(meta$Site %in% c("IU", "WU"), "USA", ""))
# meta$Subject_ID[is.na(meta$Subject_ID)] = ""
# meta$Site[is.na(meta$Site)] = ""
# 
# write.csv(meta, "../Data/meta_BAL1.csv", row.names = T)
# --- already run 

meta <- read.csv("../Data/meta_BAL1.csv", row.names=1)

## remove samples with 0 concentration
## remove samples with 0 OTUs MYCF00192
List=row.names(meta[meta$Concentration==0,])
List=cbind(List, "MYCF00192")
OTU=OTU[,-which(names(OTU) %in% List)]
meta=meta[-which(rownames(meta) %in% List),]

otunames = rownames(OTU)
# replace sequence names as otu001 otu002, ... 
row.names(OTU) <- paste0("otu",sprintf("%0.3d", seq(1, dim(OTU)[1])))
# otu.df = data.frame(otu_seq = row.names(OTU),
#                     id = paste0("otu",sprintf("%0.3d", seq(1, dim(OTU)[1]))))
# write.csv(otu.df, "./picrust2/otu_id.csv")
row.names(Tax) <- row.names(OTU)

## create phyloseq object
OTU <- otu_table(OTU, taxa_are_rows = T)
TAX <- tax_table(as.matrix(Tax))
META <- sample_data(meta)

phy <- phyloseq(OTU, TAX, META)

## Inspect library sizes
df <- as.data.frame(sample_data(phy)) # Put sample_data into a ggplot-friendly data.frame
df$LibrarySize <- sample_sums(phy)
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))
ggplot(data=df, aes(x=Index, y=LibrarySize, color=Sample_or_Control)) + geom_point()

plot_frequency(phy, taxa_names(phy)[c(1,3)], conc="Concentration") + 
  xlab("DNA Concentration (Qubit HS)")

set.seed(100)
## Identify Contaminants - Prevalence
sample_data(phy)$is.neg <- sample_data(phy)$Sample_or_Control == "Control_Sample"
contamdf.prev <- isContaminant(phy, method="prevalence", neg="is.neg")
table(contamdf.prev$contaminant) 
# FALSE  TRUE 
# 2133   27
rownames(contamdf.prev)[contamdf.prev$contaminant == T]

head(which(contamdf.prev$contaminant))
contamdf.prev05 <- isContaminant(phy, method="prevalence", neg="is.neg", threshold=0.5)
table(contamdf.prev05$contaminant)
# FALSE  TRUE 
# 2045   115 

# Make phyloseq object of presence-absence in negative controls and true samples
ps.pa <- transform_sample_counts(phy, function(abund) 1*(abund>0))
ps.pa.neg <- prune_samples(sample_data(ps.pa)$Sample_or_Control == "Control_Sample", ps.pa)
ps.pa.pos <- prune_samples(sample_data(ps.pa)$Sample_or_Control == "True_Sample", ps.pa)
# Make data.frame of prevalence in positive and negative samples
df.pa <- data.frame(pa.pos=taxa_sums(ps.pa.pos), pa.neg=taxa_sums(ps.pa.neg),
                    contaminant=contamdf.prev$contaminant)

png("publications/decontam.png", width = 8, height = 6, units = "in", res = 300)
ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)") +
  theme_classic() +
  theme(axis.text = element_text(face = "bold", size = 15),
        axis.title = element_text(face = "bold", size = 17))
dev.off()

phy.noncontam <- prune_taxa(!contamdf.prev$contaminant, phy)
phy.noncontam      # 2133 taxa
# saveRDS(phy.noncontam, "../Data/phy_noncontam.rds")

set.seed(1)
phy.noncontam.truesample = subset_samples(phy.noncontam, Sample_or_Control == "True_Sample")
phy.noncontam.truesample = rarefy_even_depth(phy.noncontam.truesample, rngseed = 1)
shannon = estimate_richness(phy.noncontam.truesample, split = TRUE, measures = "Shannon")
rownames(shannon) = meta$Sample_ID[match(rownames(shannon), rownames(meta))]
# write.csv(shannon, "../Data/shannon_BAL_decontam.csv")
simpson = estimate_richness(phy.noncontam.truesample, split = TRUE, measures = "Simpson")
rownames(simpson) = meta$Sample_ID[match(rownames(simpson), rownames(meta))]
# write.csv(simpson, "../Data/simpson_BAL_decontam.csv")
observed = estimate_richness(phy.noncontam.truesample, split = TRUE, measures = "Observed")
rownames(observed) = meta$Sample_ID[match(rownames(observed), rownames(meta))]
# write.csv(observed, "../Data/observed_BAL_decontam.csv")
Chao1 = estimate_richness(phy.noncontam.truesample, split = TRUE, measures = "Chao1")
rownames(Chao1) = meta$Sample_ID[match(rownames(Chao1), rownames(meta))]
# write.csv(Chao1, "../Data/Chao1_BAL_decontam.csv")

library(ape)
random_tree <- rtree(ntaxa(phy.noncontam), rooted = TRUE, tip.label = taxa_names(phy.noncontam))
phy.noncontam_tree <- phyloseq(OTU, TAX, META, random_tree)

### plot phylogenetic tree 

unifrac_unweighted_BAL <-UniFrac(phy.noncontam_tree, weighted = FALSE, normalized = TRUE,
                                  parallel = FALSE, fast = TRUE)

unifrac_weighted_BAL <-UniFrac(phy.noncontam_tree, weighted = TRUE, normalized = TRUE,
                                parallel = FALSE, fast = TRUE)

unifrac_unweighted_BAL <- as.matrix(unifrac_unweighted_BAL)
unifrac_weighted_BAL <- as.matrix(unifrac_weighted_BAL)
# write.csv(unifrac_unweighted_BAL, "../Data/UniFrac_unweighted_BAL.csv")
# write.csv(unifrac_weighted_BAL, "../Data/UniFrac_weighted_BAL.csv")

########## only keep the ECF samples ################# 
phy.noncontam_tree = subset_samples(phy.noncontam_tree, Sample_type == "BAL")
phy.noncontam_tree = subset_samples(phy.noncontam_tree, Study == "ECF")
#109*9 
phy.noncontam_tree@sam_data$Sequence_ID = rownames(phy.noncontam_tree@sam_data)
rownames(phy.noncontam_tree@sam_data) = phy.noncontam_tree@sam_data$Sample_ID

# png("../Results/ECF_tree_Site.png", height = 10, 
#     width = 10, units = "in", res = 300)
# plot_tree(phy.noncontam_tree, color="Site", 
#           label.tips="taxa_names", ladderize="left", plot.margin=0.3)
# dev.off()

# png("../Results/ECF_tree_Country.png", height = 10, 
#     width = 10, units = "in", res = 300)
# plot_tree(phy.noncontam_tree, color="Country", 
#           label.tips="taxa_names", ladderize="left", plot.margin=0.3)
# dev.off()

# png("../Results/ECF_tree_phylum.png", height = 10,
#     width = 10, units = "in", res = 300)
# plot_heatmap(phy.noncontam_tree, taxa.label="Phylum")
# dev.off()


plot_richness(phy.noncontam, measures="Shannon", x="Site")
plot_richness(phy.noncontam, measures="Shannon", x="Country")
# save.image("./decontam_species_v138.RData")

## keep only ECF samples 
phy.noncontam.ECF = subset_samples(phy.noncontam.truesample, 
                                   Study == "ECF")

unifrac_unweighted_ECF <-UniFrac(phy.noncontam_tree, weighted = FALSE, normalized = TRUE,
                                 parallel = FALSE, fast = TRUE)


