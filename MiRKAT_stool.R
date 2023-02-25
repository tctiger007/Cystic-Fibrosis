### within each cluster, try performing the association between pfts 
### and beta diversity using MIRKAT (needs to set up the tree and medication history)
### and alpha diversity 
### select important taxa that are associated with pfts, especially FEV0.5
library(MiRKAT)
library(GUniFrac)
library(vegan)
library(gautils2)
library(Matrix)
library(rstudioapi)
library(ape)
library(readxl)
library(microbiome)
setwd(dirname(getActiveDocumentContext()$path))

# set up tree in phyloseq object
stool.phy = readRDS("stool_phy.rds")
merge = readRDS("merge1.rds")
otu = stool.phy@otu_table
otu = t(otu)
otu_table(stool.phy) = otu
otu@taxa_are_rows

stool.phy_genus = tax_glom(stool.phy, "Genus", NArm = TRUE)
stool.phy_genus_relabun = transform_sample_counts(stool.phy_genus, 
                                                  function(OTU) OTU/sum(OTU) * 100)
stool.tree = rarefy_even_depth(stool.phy_genus_relabun, rngseed = 1)
stool.tree = rtree(ntaxa(stool.tree), rooted = TRUE, 
                 tip.label = taxa_names(stool.tree))

otu.temp = otu_table(stool.phy_genus_relabun, taxa_are_rows = TRUE)
tax.temp = tax_table(stool.phy_genus_relabun)
meta.temp = sample_data(stool.phy_genus_relabun)
# otu.temp.keep = otu.temp[match(stool.tree$tip.label, rownames(otu.temp)),]
tax.temp.keep = tax.temp[rownames(tax.temp) %in% stool.tree$tip.label,]
otu.temp.keep = otu.temp[rownames(otu.temp) %in% stool.tree$tip.label == TRUE,]
# tax.temp.keep = tax.temp[match(stool.tree$tip.label, rownames(tax.temp)),]
# rownames(otu.temp) == rownames(tax.temp)
rownames(otu.temp.keep) == rownames(tax.temp.keep)
colnames(otu.temp.keep) == rownames(meta.temp)
# stool.tree$tip.label == rownames(tax.temp.keep)
otu.temp.keep = t(otu.temp.keep)

stool.phy_genus_relabun_tree = phyloseq(otu.temp.keep, tax.temp.keep,
                                      meta.temp, stool.tree)


PFT = readRDS("../../PFT.rds")
PFT.raw = read_excel("../../../../Clinical Data/VPECF Endpoint Data for James_2019.xlsx",
                     sheet = "PFTs")
PFT.raw$`Participant ID` = str_remove(PFT.raw$`Participant ID`, "1C")
PFT.raw$`Visit ID` = ifelse(PFT.raw$`Visit ID` == "Visit 1A", "BA2", "BA7")
colnames(PFT.raw)[7] = "Visit"
colnames(PFT.raw)[1] = "SampleID"
PFT.raw = PFT.raw %>% dplyr::select(c(SampleID, Visit, 
                                      `FVC % Predicted from Kisling/Jones`,
                                      `FEV0.5 % Predicted from Kisling/Jones`,
                                      `FVC Z-Score Predicted from Kisling/Jones`,
                                      `FEV0.5 Z-Score from Kisling/Jones`))
# get the Kisling/Jones adjusted pfev0.5
PFT = PFT %>% left_join(PFT.raw, by = c("SampleID" = "SampleID",
                                        "Visit" = "Visit"))

#### Pfev0.5 and Pfvc are Lum normalized. 
### Kisling is more complete and we use Kisling adjusted values 
## 
colnames(PFT)[22:25] = c("Pfvc_Kisling", "Pfev0.5_Kisling",
                         "Zfvc_Kisling", "Zfev0.5_Kisling")

### using Lum adjusted PFTs 
sum(!(is.na(PFT$Pfev0.5))) # 59
PFT = PFT[-which((is.na(PFT$Pfev0.5))),]
PFT = PFT %>% mutate(Sample_ID = paste0(PFT$SampleID, PFT$Visit))
PFT = PFT %>% inner_join(merge, by = "Sample_ID")
stool.phy_genus_relabun_tree = subset_samples(stool.phy_genus_relabun_tree,
                                              Sample_ID %in% PFT$Sample_ID_stool)
meta = meta(stool.phy_genus_relabun_tree)
meta$Sequence_ID = rownames(meta)
meta = meta %>% left_join(PFT, by = "Sequence_ID")

meta$Sequence_ID == rownames(stool.phy_genus_relabun_tree@sam_data)
rownames(meta) = meta$Sequence_ID
stool.phy_genus_relabun_tree = phy_substitute_metadata(stool.phy_genus_relabun_tree, 
                                                       new_metadata = meta)

stool.phy_genus_relabun_tree = subset_samples(stool.phy_genus_relabun_tree, 
                                              Sequence_ID %in% 
                                                rownames(meta)[!is.na(meta$Gender.x)])

meta = meta(stool.phy_genus_relabun_tree)
otu = otu_table(stool.phy_genus_relabun_tree)
tree = stool.phy_genus_relabun_tree@phy_tree
unifracs = GUniFrac(otu, tree, 
                    alpha = c(0, 0.5, 1))$unifracs

D.weighted = unifracs[,,"d_1"]
D.unweighted = unifracs[,,"d_UW"]
D.generalized = unifracs[,,"d_0.5"]
D.BC = as.matrix(vegdist(otu, method="bray"))

K.weighted = D2K(D.weighted)
K.unweighted = D2K(D.unweighted)
K.generalized = D2K(D.generalized)
K.BC = D2K(D.BC)

Ks = list(K.weighted = K.weighted, K.unweighted = K.unweighted, 
          K.generalized = K.generalized, K.BC = K.BC)
meta$Sequence_ID == rownames(D.weighted)
PFT$Sequence_ID == meta$Sequence_ID
PFT = PFT[match(meta$Sequence_ID, PFT$Sequence_ID),]
PFT$Sequence_ID == meta$Sequence_ID
PFT$Sequence_ID == rownames(D.weighted)



PFT1 = PFT[,which(colnames(PFT) %in% c("SampleID.x", "DOB",
                                       "Date","Stool_Age", "Pfev0.5.x"))]

# write.csv(PFT1, "PFT_with_stool.csv")
PFT1 = read.csv("PFT_with_stool.csv")
PFT1$Male = as.numeric(meta$Gender.x == "Male")
load("../../clinical.RData")
PFT1$F508 = demo$Genotype_F508[match(PFT1$SampleID.x, demo$Study_ID)]
PFT1$F508 = ifelse(PFT1$F508 == "Homozygous",0, 1)
PFT1$Country = meta$Country.x
PFT1$AUS = ifelse(PFT1$Country == "Aus",1,0)
covar = PFT1[,c(12:15,6:11)]
covar = covar[,-3]
covar = covar[,-3]


Matrix::rankMatrix(covar)

set.seed(123)
MiRKAT(y = PFT$Pfev0.5.x, 
       Ks = Ks, 
       returnKRV = TRUE, returnR2 = TRUE,
       out_type = "C", nperm = 999, 
       method = "davies", omnibus = "cauchy")
# $p_values
# K.weighted  K.unweighted K.generalized          K.BC 
# 0.16411268    0.02348627    0.10434677    0.01933755 
# 
# $omnibus_p
# [1] 0.03655878

# #covar1 combined PE and PE_ad; 
# #covar2 combined PE and PE_ad and removed PX_Pseudomonas since it's only in one sample;
# #covar3 removed PX_Pseudomonas since it's only in one sample;
MiRKAT(y = PFT$Pfev0.5.x,
       X = covar[,1:2], Ks = Ks,
       returnKRV = TRUE, returnR2 = TRUE,
       out_type = "C", nperm = 999,
       method = "davies", omnibus = "cauchy")
# $p_values
# K.weighted  K.unweighted K.generalized          K.BC 
# 0.17969655    0.04043132    0.10578718    0.05672689 
# 
# $omnibus_p
# [1] 0.07013784
MiRKAT(y = PFT$Pfev0.5.x,
       X = covar, Ks = Ks,
       returnKRV = TRUE, returnR2 = TRUE,
       out_type = "C", nperm = 999,
       method = "davies", omnibus = "cauchy")
# $p_values
# K.weighted  K.unweighted K.generalized          K.BC 
# 0.4968673     0.4064731     0.5863714     0.2520555 
# 
# $omnibus_p
# [1] 0.4204212


MiRKAT(y = PFT$Pfvc.x, 
       Ks = Ks, 
       returnKRV = TRUE, returnR2 = TRUE,
       out_type = "C", nperm = 999, 
       method = "davies", omnibus = "cauchy")

# $p_values
# K.weighted  K.unweighted K.generalized          K.BC 
# 0.26541148    0.20257212    0.20707646    0.08156818 
# 
# $omnibus_p
# [1] 0.157895

MiRKAT(y = PFT$Pfvc.x,
       X = covar[,1:2], Ks = Ks,
       returnKRV = TRUE, returnR2 = TRUE,
       out_type = "C", nperm = 999,
       method = "davies", omnibus = "cauchy")
# $p_values
# K.weighted  K.unweighted K.generalized          K.BC 
# 0.3363144     0.3058821     0.2787320     0.1708339 
# 
# $omnibus_p
# [1] 0.2590636

MiRKAT(y = PFT$Pfev0.5.x,
       X = covar, Ks = Ks,
       returnKRV = TRUE, returnR2 = TRUE,
       out_type = "C", nperm = 999,
       method = "davies", omnibus = "cauchy")
# $p_values
# K.weighted  K.unweighted K.generalized          K.BC 
# 0.4968673     0.4064731     0.5863714     0.2520555 
# 
# $omnibus_p
# [1] 0.4204212

