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
library(latex2exp)
library(ggpmisc)

# wesanderson::wes_palette("FantasticFox1")
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

PFT = readRDS("../PFT.rds")
PFT.raw = read_excel("../../../Clinical Data/VPECF Endpoint Data for James_2019.xlsx",
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
PFTBA2 = PFT %>% filter(Visit == "BA2")
PFTBA7 = PFT %>% filter(Visit == "BA7")
PFTBA2.ID = PFTBA2$SampleID  # 28
PFTBA7.ID = PFTBA7$SampleID  #31
length(intersect(PFTBA2.ID, PFTBA7.ID))  # 21
# which(PFT$SampleID %in% names(which(table(PFT$SampleID) == 1)))
# # Absence of PFTs at either BA2 or BA7 visit
# temp = PFT[-which(PFT$SampleID %in% names(which(table(PFT$SampleID) == 1))),] #54
# # Absence of either BA2 or BA7 entry
# temp = temp[-which(is.na(temp$Pfev0.5)),] #48 


commonID.T1T2 = intersect(PFTBA2.ID, intersect(PFTBA7.ID, intersect(BA2.subID, BA7.subID))) #19
commonID.T2 = intersect(PFTBA7.ID, BA7.subID) # 28

BAL.phy = subset_samples(BAL.phy, Subject_ID %in% commonID.T1T2)
meta = sample_data(BAL.phy)
otu = otu_table(BAL.phy)
tax = tax_table(BAL.phy)
colnames(otu) == rownames(meta)
rownames(meta) = meta$Sample_ID
colnames(otu) = meta$Sample_ID
colnames(otu) == rownames(meta)
BAL.phy = phyloseq(otu, tax, meta)

##### clustering on BA2 + BA7 samples 
### reference https://susan.su.domains/papers/Pregnancy/PNAS_Vaginal_Analysis.Rmd
# BAL.phy = readRDS("matched_CT_BAL_keep_IU121/matched_CT_BAL_2timepoints.rds")
BAL.phy_genus = tax_glom(BAL.phy, "Genus", NArm = TRUE)
BAL.phy_genus_relabun = transform_sample_counts(BAL.phy_genus, function(OTU) OTU/sum(OTU) * 100)
# saveRDS(BAL.phy_genus_relabun, "BAL.phy_genus_relabun.rds")
braydist = phyloseq::distance(BAL.phy_genus_relabun, method="bray")
ord = ordinate(BAL.phy_genus_relabun, method = "MDS", distance = braydist)
# png("ord_eigenvalues.png", width = 8, height = 8, units = "in", res = 300)
print(plot_scree(ord)) + theme_classic() +
  ggtitle("MDS-bray ordination eigenvalues")
# dev.off()
evs <- ord$value$Eigenvalues
print(evs[1:20])
print(tail(evs))

h_sub5 <- hist(evs[6:length(evs)], 100)
plot(h_sub5$mids, h_sub5$count, log="y", type='h', lwd=10, lend=2)

NDIM <- 27
x <- ord$vectors[,1:NDIM]  # rows=sample, cols=MDS axes, entries = value
pamPCoA = function(x, k) {
  list(cluster = pam(x[,1:2], k, cluster.only = TRUE))
}
gs = clusGap(x, FUN = pamPCoA, K.max = 12, B = 50)
# png("gap_stat.png", width = 8, height = 8, units = "in", res = 300)
plot_clusgap(gs) + scale_x_continuous(breaks=c(seq(0, 12, 2))) + 
  theme_classic()
# dev.off()
#  We select 3 clusters based on the gap statistics 

## Perform PAM 3-fold clusters:
K <- 3
x <- ord$vectors[,1:NDIM]
clust <- as.factor(pam(x, k=K, cluster.only=T))

sample_data(BAL.phy_genus_relabun)$Cluster <- clust
Clusters <- as.character(seq(K))

ClusterColors = brewer.pal(6,"Paired")[c(1,3,5)] 
names(ClusterColors) = c("Cluster1", "Cluster2", "Cluster3")
CSTFillScale = scale_fill_manual(values = c("1"="#A6CEE3","2"="#B2DF8A","3"="#FB9A99"),
                                 labels = c("1", "2", "3"))
plot_ordination(BAL.phy_genus_relabun, ord, color="Cluster") 
braydist <- phyloseq::distance(BAL.phy_genus_relabun, method="bray")

# png("ord_genus_MDS_bray_3cluster.png", width = 8, height = 8, units = "in", res=300)
plot_ordination(BAL.phy_genus_relabun, 
                ordinate(BAL.phy_genus_relabun, method="MDS", distance=braydist), 
                color="Cluster") + ggtitle("MDS -- bray") +
  geom_point(size = 6) + 
  theme_classic() +
  scale_color_manual(values=c("orange2", "olivedrab4", "mediumorchid")) + 
  theme(axis.title = element_text(size = 20, face = "bold"),
        axis.text = element_text(size = 18),
        plot.title = element_text(size = 18, face = "bold"),
        legend.title = element_text(size = 18, face = "bold"),
        legend.text = element_text(size = 18, face = "bold"))
# dev.off()

taxa.order <- names(sort(taxa_sums(BAL.phy_genus_relabun)))
tax = tax_table(BAL.phy_genus_relabun)
taxa.order = tax[match(taxa.order, rownames(tax)),"Genus"]
taxa.order1 = rownames(taxa.order) # reversed order of the most abund taxa
pshm = prune_taxa(names(sort(taxa_sums(BAL.phy_genus_relabun), T))[1:25], 
                   BAL.phy_genus_relabun)

pshm.C1 = prune_samples(sample_data(pshm)$Cluster == 1, pshm)
pshm.C2 = prune_samples(sample_data(pshm)$Cluster == 2, pshm)
pshm.C3 = prune_samples(sample_data(pshm)$Cluster == 3, pshm)
set.seed(1)
sample.order1 = sort(rownames(sample_data(pshm.C1)))
sample.order2 = sort(rownames(sample_data(pshm.C2)))
sample.order3 = sort(rownames(sample_data(pshm.C3)))

p.c1 = plot_heatmap(pshm.C1, taxa.label="Genus", taxa.order=taxa.order1,
                    sample.order = sample.order1) + 
  ggtitle("Cluster: 1")

p.c2 = plot_heatmap(pshm.C2, taxa.label="Genus", taxa.order=taxa.order1,
                    sample.order = sample.order2) + 
  ggtitle("Cluster: 2")
p.c3 = plot_heatmap(pshm.C3, taxa.label="Genus", taxa.order=taxa.order1,
                    sample.order = sample.order3) + 
  ggtitle("Cluster: 3")

sample.order = c(sample.order1, sample.order2,
                 sample.order3)
pshm = merge_phyloseq(pshm.C1, pshm.C2, pshm.C3)
pshm = ps_reorder(pshm, sample.order)
taxa.order2 = rev(taxa.order1)[1:25]
pshm.tax = tax_table(pshm)
pshm.otu = otu_table(pshm)
pshm.tax = pshm.tax[match(taxa.order2, rownames(pshm.tax)),]
pshm.otu = pshm.otu[match(taxa.order2, rownames(pshm.otu)),]
pshm.meta = sample_data(pshm)
pshm = phyloseq(pshm.tax, pshm.otu, pshm.meta)

# taxa.order.c = names(sort(taxa_sums(pshm.C1), decreasing = T))
heatmap.df = as.data.frame(otu_table(pshm, taxa_are_rows = T))
rownames(heatmap.df) = pshm.tax[match(rownames(heatmap.df), rownames(pshm.tax)), "Genus"]

meta2 = sample_data(pshm)
meta2 = as.data.frame(as.matrix(meta2))
PFT = PFT %>% mutate(Sample_ID = paste0(PFT$SampleID, PFT$Visit))
meta2 = meta2 %>% inner_join(PFT, by = "Sample_ID")
meta2 = meta2[,-grep(".y", colnames(meta2), fixed = T)]
colnames(meta2) = str_remove_all(colnames(meta2), ".x")

# load("../clinical.RData")
my_sample_col = meta2 %>% dplyr::select(c("Country","Visit","Cluster","Pfev0.5", "Pfvc"))
colnames(my_sample_col)[4] = "FEV0.5%" 
colnames(my_sample_col)[5] = "FVC%" 
rownames(my_sample_col) = meta2$Sample_ID

my_phylum_col = tax_table(pshm)
my_phylum_col = my_phylum_col[,"Phylum"]
my_phylum_col = as.data.frame(my_phylum_col)
rownames(my_phylum_col) = tax_table(pshm)[,"Genus"]
ann_colors = list(
  Phylum = c(Actinobacteriota="#8c8c8c", Bacteroidota="#ffc425",
             Campylobacterota="#f37735", Firmicutes="#00aedb",
             Fusobacteriota="#00b159", Myxococcota="#FC717F",
             Proteobacteria="#d11141"),
  Cluster = c(`1`="orange2", `2`="olivedrab4", `3`="mediumorchid"),
  Country = c(USA="#00AD9A", AUS="#E16A86"),
  Visit = c(BA2="gold2", BA7="hotpink4"))
# saveRDS(ann_colors, "../publications/publication_picrust2/ann_colors.rds")
# saveRDS(my_sample_col, "../publications/publication_picrust2/my_sample_col.rds")


dev.off()
# png("./cluster_heatmap.png", width = 12, height = 12,
#     units = "in", res = 300)
pheatmap(heatmap.df,
         color=colorRampPalette(c("#000033", "#66CCFF"))(50),
         cluster_rows=F, cluster_cols=F,
         annotation_row = my_phylum_col,
         annotation_col = my_sample_col,
         border_color = "NA",
         annotation_colors = ann_colors,
         gaps_col = c(15, 29),
         fontsize = 10.5) +
  theme_classic()
# dev.off()

my_sample_col %>% group_by(Cluster) %>% summarise(mean = mean(`FEV0.5%`))
# A tibble: 3 × 2
# Cluster  mean
# <chr>   <dbl>
#   1 1        90.8
# 2 2       101. 
# 3 3       111. 
my_sample_col$sampleID = rownames(my_sample_col)
my_sample_col$sampleID = factor(my_sample_col$sampleID, levels = sample.order)
mean.pFev = my_sample_col %>% group_by(Cluster) %>% 
  summarise(mean = mean(`FEV0.5%`, na.rm = T))
mean.pFev.c1 = mean.pFev$mean[1]
mean.pFev.c2 = mean.pFev$mean[2]
mean.pFev.c3 = mean.pFev$mean[3]

my_sample_col$Country = factor(my_sample_col$Country)
# saveRDS(my_sample_col, "cluster.rds")
# ggplot(my_sample_col, aes(x=sampleID, y=zFev, color = Country)) +
#   geom_point() +
#   geom_segment(aes(x=1,xend=length(sample.order1),
#                    y=mean.zFev.c1,yend=mean.zFev.c1)) +
#   geom_segment(aes(x=length(sample.order1)+1,xend=length(sample.order1)+length(sample.order2),
#                    y=mean.zFev.c2,yend=mean.zFev.c2)) +
#   geom_segment(aes(x=length(sample.order1)+length(sample.order2)+1,
#                    xend=length(sample.order),
#                    y=mean.zFev.c3,yend=mean.zFev.c3)) +
#   theme_classic()
#   
my_sample_col$Country = factor(my_sample_col$Country,
                               levels = c("USA", "AUS"))
# png("pFEV_3cluster.png", width = 8, height = 8, units = "in", res=300)
ggplot(my_sample_col, aes(x=Cluster, y=`FEV0.5%`, fill = Cluster)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(shape = Country), size = 3) +
  scale_fill_manual(values=c("orange2", "olivedrab4", "mediumorchid")) +
  scale_shape_manual(values=c(16, 17))+
  annotate("segment", x = 1, xend = 3, y = 160, yend = 160) +
  annotate("text",label="*", x = 2, y = 162, size = 10) +
  theme_classic() +
  # ylab(bquote(bold(paste(FEV[0.5],"%")))) + 
  ylab(TeX(r'(FEV$_{0.5}$%)', 
           bold=TRUE)) +
  geom_hline(yintercept = 80, linetype = "dashed") +
  geom_hline(yintercept = 120, linetype = "dashed") +
 theme(axis.title = element_text(size = 22, face = "bold"),
      axis.text = element_text(size = 20),
      plot.title = element_text(size = 20, face = "bold"),
      legend.title = element_text(size = 20, face = "bold"),
      legend.text = element_text(size = 20, face = "bold"))
# ,
      # legend.position = "none")
# dev.off()


# png("pFVC_3cluster.png", width = 8, height = 8, units = "in", res=300)
ggplot(my_sample_col, aes(x=Cluster, y=`FVC%`, fill = Cluster)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(shape = Country), size = 3) +
  scale_fill_manual(values=c("orange2", "olivedrab4", "mediumorchid")) +
  scale_shape_manual(values=c(16, 17))+
  theme_classic() +
  geom_hline(yintercept = 80, linetype = "dashed") +
  geom_hline(yintercept = 120, linetype = "dashed") +
  theme(axis.title = element_text(size = 22, face = "bold"),
        axis.text = element_text(size = 20),
        plot.title = element_text(size = 20, face = "bold"),
        legend.title = element_text(size = 20, face = "bold"),
        legend.text = element_text(size = 20, face = "bold"))#,
        # legend.position = "none")
# dev.off()

aov.fev = aov(my_sample_col$`FEV0.5%`~ my_sample_col$Cluster)
summary(aov.fev) # p = 0.0213 * 

aov_residuals0 <- residuals(object = aov.fev)
shapiro.test(x = aov_residuals0) # 0.4222
TukeyHSD(aov.fev)$`my_sample_col$Cluster`
# diff        lwr      upr     p adj
# 2-1 0.8177954 -0.3691376 2.004728 0.2246613
# 3-1 1.5912405  0.2445259 2.937955 0.0175465
# 3-2 0.7734452 -0.5911866 2.138077 0.3584833


aov.fvc = aov(my_sample_col$`FVC%`~ my_sample_col$Cluster)
summary(aov.fvc) # p = 0.389 
aov_residuals1 <- residuals(object = aov.fvc)
shapiro.test(x = aov_residuals1) # 0.046
TukeyHSD(aov.fvc)$`my_sample_col$Cluster`

### test whether samples from different countries fall into unique clusters
table(my_sample_col$Country[my_sample_col$Cluster==1]) #15
# AUS USA 
# 12   3 
table(my_sample_col$Country[my_sample_col$Cluster==2]) #14
# AUS USA 
# 5   9 
table(my_sample_col$Country[my_sample_col$Cluster==3]) #9
# AUS USA 
# 1   8 
save(my_sample_col, file = "my_sample_col.RData")
cluster.num = as.table(rbind(c(3,9,8),
                             c(12,5,1)))
dimnames(cluster.num) = list(Country = c("USA", "AUS"),
                             Cluster = c("One","Two","Three"))
  
cluster.num
chisq.test(cluster.num, correct=F, simulate.p.value = T)
# X-squared = 13.171, df = NA, p-value = 0.0009995
# X-squared = 11.915, df = NA, p-value = 0.002999

library(rstatix)
chisq_test(cluster.num, correct=F, simulate.p.value = T)
pairwise_prop_test(cluster.num, 
                   p.adjust.method = "BH")
# group1 group2       p  p.adj p.adj.signif
# * <chr>  <chr>    <dbl>  <dbl> <chr>       
#   1 One    Two    0.0411  0.0617 ns          
#   2  One    Three  0.00429 0.0129 *           
#   3 Two    Three  0.409   0.409  ns 



sample.order.transit = unlist(lapply(str_split(sample.order, "BA"), "[[",1))
sample.order.transit = unique(sample.order.transit)
meta2$Visit = as.factor(meta2$Visit)
# meta2$Subject_ID = factor(meta2$Subject_ID,
#                          levels = c(sample.order.transit))
## study the transition of clusters for each subject 
# png("transition.png", width = 6, height = 8,
#     units = "in", res = 300)
ggplot(meta2, aes(x=Visit, y = Subject_ID, color = Cluster)) +
  geom_point(size = 5)+
  xlab("")+
  ylab("")+
  scale_color_manual(values=c("orange2", "olivedrab4", "mediumorchid")) + 
  theme_classic() +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 16))
# dev.off()


# png("transition_alluvial.png", width = 6, height = 8,
#     units = "in", res = 300)
ggplot(meta2,
       aes(x = Visit, stratum = Cluster, alluvium = Subject_ID,
           fill = Cluster, label = Cluster)) +
  scale_fill_manual(values=c("orange2", "olivedrab4", "mediumorchid")) + 
  geom_flow(stat = "alluvium", lode.guidance = "frontback",
            color = "darkgray") +
  geom_stratum() +
  theme_classic() +
  ggtitle("Cluster transition between two visits") +
  theme(axis.title = element_text(size = 20, face = "bold"),
        axis.text = element_text(size = 18),
        plot.title = element_text(size = 18, face = "bold"),
        legend.title = element_text(size = 18, face = "bold"),
        legend.text = element_text(size = 18, face = "bold"))
# dev.off()

# 
# ### try Markov train
pshm <- prune_taxa(names(sort(taxa_sums(BAL.phy_genus_relabun), T))[1:25],
                   BAL.phy_genus_relabun)
# View(sample_data(pshm))
meta = sample_data(pshm)
meta = as.data.frame(do.call(cbind,meta@.Data))
colnames(meta) = colnames(sample_data(pshm))
meta.BA2 = meta %>% dplyr::filter(Visit == "BA2")
meta.BA7 = meta %>% dplyr::filter(Visit == "BA7")
meta.wide = meta.BA2 %>% inner_join(meta.BA7, by = "Subject_ID")
meta.wide$Cluster.x = as.numeric(meta.wide$Cluster.x)
meta.wide$Cluster.y = as.numeric(meta.wide$Cluster.y)

ttab <- table(meta.wide$Cluster.x, meta.wide$Cluster.y) # prevstate=row, curstate=col
trans <- matrix(ttab, 3)
trans <- trans/rowSums(trans)  # Normalize row sums to 1
CSTtrans <- trans
CSTs <- c("Cluster1","Cluster2","Cluster3")
colnames(CSTtrans) <- CSTs
rownames(CSTtrans) <- CSTs
t_persist <- -1/log(diag(CSTtrans))
CSTtrans 
#            Cluster1  Cluster2  Cluster3
# Cluster1 0.3333333 0.5555556 0.1111111
# Cluster2 0.1666667 0.3333333 0.5000000
# Cluster3 0.5000000 0.2500000 0.2500000

meta.USA = meta %>% filter(Country == "USA")
meta.USA.BA2 = meta.USA %>% dplyr::filter(Visit == "BA2")
meta.USA.BA7 = meta.USA %>% dplyr::filter(Visit == "BA7")
meta.USA.wide = meta.USA.BA2 %>% inner_join(meta.USA.BA7, by = "Subject_ID")
meta.USA.wide$Cluster.x = as.numeric(meta.USA.wide$Cluster.x)
meta.USA.wide$Cluster.y = as.numeric(meta.USA.wide$Cluster.y)

ttab.USA <- table(meta.USA.wide$Cluster.x, meta.USA.wide$Cluster.y) # prevstate=row, curstate=col
trans.USA <- matrix(ttab.USA, 3)
trans.USA <- trans.USA/rowSums(trans.USA)  # Normalize row sums to 1
CSTtrans.USA <- trans.USA
CSTs <- c("Cluster1","Cluster2","Cluster3")
colnames(CSTtrans.USA) <- CSTs
rownames(CSTtrans.USA) <- CSTs
t_persist.USA <- -1/log(diag(CSTtrans.USA))
CSTtrans.USA # Paper
#             Cluster1 Cluster2 Cluster3
# Cluster1      0.0     1.00     0.00
# Cluster2      0.0     0.40     0.60
# Cluster3      0.5     0.25     0.25
t_persist.USA # Paper
# Cluster1  Cluster2  Cluster3 
# 0.0000000 1.0913567 0.7213475 
# summary(meta.USA.wide$GDColl-samdf[meta.USA.wide$PrevID,"GDColl"]) # Paper


meta.AUS = meta %>% filter(Country == "AUS")
meta.AUS.BA2 = meta.AUS %>% dplyr::filter(Visit == "BA2")
meta.AUS.BA7 = meta.AUS %>% dplyr::filter(Visit == "BA7")
meta.AUS.wide = meta.AUS.BA2 %>% inner_join(meta.AUS.BA7, by = "Subject_ID")
meta.AUS.wide$Cluster.x = as.numeric(meta.AUS.wide$Cluster.x)
meta.AUS.wide$Cluster.y = as.numeric(meta.AUS.wide$Cluster.y)

ttab.AUS <- table(meta.AUS.wide$Cluster.x, meta.AUS.wide$Cluster.y) # prevstate=row, curstate=col
trans.AUS <- matrix(ttab.AUS, ncol = 3, byrow = F)
trans.AUS <- rbind(trans.AUS, c(0,0,0))
trans.AUS <- trans.AUS/rowSums(trans.AUS)  # Normalize row sums to 1
CSTtrans.AUS <- trans.AUS
CSTs.AUS <- c("Cluster1","Cluster2","Cluster3")
colnames(CSTtrans.AUS) <- CSTs
rownames(CSTtrans.AUS) <- CSTs
t_persist.AUS <- -1/log(diag(CSTtrans.AUS))
CSTtrans.AUS # Paper
# Cluster1 Cluster2 Cluster3
# Cluster1    0.375      0.5    0.125
# Cluster2    1.000      0.0    0.000
# Cluster3      NaN      NaN      NaN
t_persist.AUS
# Cluster1  Cluster2  Cluster3
# 1.019545 0.000000      NaN 


### within each cluster, try performing the association between pfts 
### and beta diversity using MIRKAT (needs to set up the tree and medication history)
### and alpha diversity 
### select important taxa that are associated with pfts, especially FEV0.5
library(MiRKAT)
library(GUniFrac)
library(vegan)
library(gautils2)
library(Matrix)
# set up tree in phyloseq object
BAL.tree = rarefy_even_depth(BAL.phy_genus_relabun, rngseed = 1)
BAL.tree = rtree(ntaxa(BAL.tree), rooted = TRUE, 
                 tip.label = taxa_names(BAL.tree))
otu.temp = otu_table(BAL.phy_genus_relabun, taxa_are_rows = TRUE)
tax.temp = tax_table(BAL.phy_genus_relabun)
meta.temp = sample_data(BAL.phy_genus_relabun)
# otu.temp.keep = otu.temp[match(BAL.tree$tip.label, rownames(otu.temp)),]
tax.temp.keep = tax.temp[rownames(tax.temp) %in% BAL.tree$tip.label,]
otu.temp.keep = otu.temp[rownames(otu.temp) %in% BAL.tree$tip.label == TRUE,]
# tax.temp.keep = tax.temp[match(BAL.tree$tip.label, rownames(tax.temp)),]
# rownames(otu.temp) == rownames(tax.temp)
rownames(otu.temp.keep) == rownames(tax.temp.keep)
colnames(otu.temp.keep) == rownames(meta.temp)
# BAL.tree$tip.label == rownames(tax.temp.keep)
otu.temp.keep = t(otu.temp.keep)
BAL.phy_genus_relabun_tree = phyloseq(otu.temp.keep, tax.temp.keep, 
                                      meta.temp, BAL.tree)
 
unifracs = GUniFrac(otu.temp.keep, BAL.tree, 
                    alpha = c(0, 0.5, 1))$unifracs
D.weighted = unifracs[,,"d_1"]
D.unweighted = unifracs[,,"d_UW"]
D.generalized = unifracs[,,"d_0.5"]
D.BC = as.matrix(vegdist(otu.temp.keep, method="bray"))

D.weighted.c1 = D.weighted[rownames(D.weighted) %in% sample.order1, colnames(D.weighted) %in% sample.order1]
D.weighted.c2 = D.weighted[rownames(D.weighted) %in% sample.order2, colnames(D.weighted) %in% sample.order2]
D.weighted.c3 = D.weighted[rownames(D.weighted) %in% sample.order3, colnames(D.weighted) %in% sample.order3]
D.unweighted.c1 = D.unweighted[rownames(D.unweighted) %in% sample.order1, colnames(D.unweighted) %in% sample.order1]
D.unweighted.c2 = D.unweighted[rownames(D.unweighted) %in% sample.order2, colnames(D.unweighted) %in% sample.order2]
D.unweighted.c3 = D.unweighted[rownames(D.unweighted) %in% sample.order3, colnames(D.unweighted) %in% sample.order3]
D.generalized.c1 = D.generalized[rownames(D.generalized) %in% sample.order1, colnames(D.generalized) %in% sample.order1]
D.generalized.c2 = D.generalized[rownames(D.generalized) %in% sample.order2, colnames(D.generalized) %in% sample.order2]
D.generalized.c3 = D.generalized[rownames(D.generalized) %in% sample.order3, colnames(D.generalized) %in% sample.order3]
D.BC.c1 = D.BC[rownames(D.BC) %in% sample.order1, colnames(D.BC) %in% sample.order1]
D.BC.c2 = D.BC[rownames(D.BC) %in% sample.order2, colnames(D.BC) %in% sample.order2]
D.BC.c3 = D.BC[rownames(D.BC) %in% sample.order3, colnames(D.BC) %in% sample.order3]

K.weighted = D2K(D.weighted)
K.unweighted = D2K(D.unweighted)
K.generalized = D2K(D.generalized)
K.BC = D2K(D.BC)

K.weighted.c1 = D2K(D.weighted.c1)
K.unweighted.c1 = D2K(D.unweighted.c1)
K.generalized.c1 = D2K(D.generalized.c1)
K.BC.c1 = D2K(D.BC.c1)
K.weighted.c2 = D2K(D.weighted.c2)
K.unweighted.c2 = D2K(D.unweighted.c2)
K.generalized.c2 = D2K(D.generalized.c2)
K.BC.c2 = D2K(D.BC.c2)
K.weighted.c3 = D2K(D.weighted.c3)
K.unweighted.c3 = D2K(D.unweighted.c3)
K.generalized.c3 = D2K(D.generalized.c3)
K.BC.c3 = D2K(D.BC.c3)

Ks = list(K.weighted = K.weighted, K.unweighted = K.unweighted, 
          K.generalized = K.generalized, K.BC = K.BC)
Ks.c1 = list(K.weighted = K.weighted.c1, K.unweighted = K.unweighted.c1, K.BC = K.BC.c1)
Ks.c2 = list(K.weighted = K.weighted.c2, K.unweighted = K.unweighted.c2, K.BC = K.BC.c2)
Ks.c3 = list(K.weighted = K.weighted.c3, K.unweighted = K.unweighted.c3, K.BC = K.BC.c3)

load("../clinical.RData")
colnames(demo)[which(colnames(demo) == "SampleID")] = "Sample_ID"
meta.temp1 = meta.temp %>% data.frame() %>% inner_join(demo, by = "Sample_ID")
rownames(meta.temp1) = meta.temp1$Sample_ID
meta.temp = meta.temp1 
meta.temp = meta.temp[,-grep(".y", colnames(meta.temp), fixed = T)]
colnames(meta.temp)[grep(".x", colnames(meta.temp), fixed = T)] = 
  stringr::str_remove(colnames(meta.temp)[grep(".x", colnames(meta.temp), fixed = T)], ".x")

meta.temp$Group = paste0(meta.temp$Country, meta.temp$Visit)
meta.temp$Group = factor(meta.temp$Group, levels = c("USABA2", "USABA7",
                                                     "AUSBA2", "AUSBA7"))

# set up medication and add it into the meta data 
Male = as.numeric(meta.temp$Gender == "Male")
F508 = ifelse(meta.temp$Genotype_F508 == "Homozygous",0,
              ifelse(meta.temp$Genotype_F508 == "Heterogenous",1,2))
Group = as.numeric(meta.temp$Group)

load("../onlyT2/ABX_accumulate.RData") # 58*3
ABX.freq.accumulate.T1 = ABX.freq.accumulate.T1[ABX.freq.accumulate.T1$Subject %in% commonID.T1T2,]
ABX.freq.accumulate.T2 = ABX.freq.accumulate.T2[ABX.freq.accumulate.T2$Subject %in% commonID.T1T2,]
colnames(meta)[which(colnames(meta) == "Subject_ID")] = "Subject"
meta.1 = as.data.frame(meta)
ABX.freq.accumulate.T1.wide = ABX.freq.accumulate.T1 %>% spread(SubType,n)
ABX.freq.accumulate.T2.wide = ABX.freq.accumulate.T2 %>% spread(SubType,n)
for(j in 1:2){
  df.name = paste0("ABX.freq.accumulate.T",j,".wide")
  # eval(parse(df.name))
  df = get(df.name)
  na.idx = which(is.na(df), arr.ind = T)
  for (i in 1:dim(na.idx)[1]){
    row = na.idx[i,1]
    col = na.idx[i,2]
    df[row, col] = 0
  }
  colnames(df)[2:7] = paste0(colnames(df)[2:7], ".T", j)
  assign(df.name, df)
}

ABX.freq.accumulate.T1.wide$Sample_ID = paste0(ABX.freq.accumulate.T1.wide$Subject,"BA2")
ABX.freq.accumulate.T2.wide$Sample_ID = paste0(ABX.freq.accumulate.T2.wide$Subject,"BA7")
colnames(ABX.freq.accumulate.T1.wide) = str_remove_all(colnames(ABX.freq.accumulate.T1.wide), ".T1")
colnames(ABX.freq.accumulate.T2.wide) = str_remove_all(colnames(ABX.freq.accumulate.T2.wide), ".T2")

# rankMatrix(ABX.freq.accumulate.T1.wide[,-c(1,9)]) ## 4; not full rank;
# rankMatrix(ABX.freq.accumulate.T2.wide[,-c(1,9)]) ## 7; full rank

meta.temp1 = meta.temp %>% data.frame() %>% 
  left_join(ABX.freq.accumulate.T1.wide, by = "Sample_ID") 

meta.temp1[match(ABX.freq.accumulate.T2.wide$Sample_ID, meta.temp1$Sample_ID),
           28:34] = ABX.freq.accumulate.T2.wide[,1:7]
meta.temp = meta.temp1

medication = meta.temp[,c(29:34)]
rankMatrix(medication)
dim(medication)

colnames(meta.temp)[grep("-", colnames(meta.temp))] = 
  str_replace(colnames(meta.temp)[grep("-", colnames(meta.temp))], "-", "_")

meta.1 = meta.temp
Male = as.numeric(meta.1$Gender == "Male")
F508 = ifelse(meta.1$Genotype_F508 == "Homozygous",0,
              ifelse(meta.1$Genotype_F508 == "Heterogenous",1,2))
Group = as.numeric(meta.1$Group)
medication = meta.1[,29:35]
medication[is.na(medication$Other),dim(medication)[2]] = 0
rownames(meta.1) = meta.1$Sample_ID

BAL.phy_genus_relabun_tree = phy_substitute_metadata(BAL.phy_genus_relabun_tree, meta.1)

covar = cbind(Male, F508, Group, medication)
## PX_Staphylococcus is nested in Group, we remove it 
covar = covar[,-which(colnames(covar) == "PX_Staphylococcus")]
rownames(covar) = rownames(meta.1)
PFT = PFT[match(rownames(covar), PFT$Sample_ID),]
PFT$Sample_ID == rownames(covar)

set.seed(123)
# FEV0.5 is significantly associated with beta diversity 
MiRKAT(y = PFT$Pfev0.5, 
       Ks = Ks, 
       returnKRV = TRUE, returnR2 = TRUE,
       out_type = "C", nperm = 999, 
       method = "davies", omnibus = "cauchy")
covar1 = covar
covar1$PE_total = covar$PE + covar$PE_ad
covar1 = covar1[,-c(4,5)]
covar2 = covar1[,-5]
covar3 = covar[,-7]
#covar1 combined PE and PE_ad; 
#covar2 combined PE and PE_ad and removed PX_Pseudomonas since it's only in one sample;
#covar3 removed PX_Pseudomonas since it's only in one sample;
# MiRKAT(y = PFT$Pfev0.5, 
#        X = covar, Ks = Ks, 
#        returnKRV = TRUE, returnR2 = TRUE,
#        out_type = "C", nperm = 999, 
#        method = "davies", omnibus = "cauchy")
MiRKAT(y = PFT$Pfev0.5, 
       X = covar1, Ks = Ks, 
       returnKRV = TRUE, returnR2 = TRUE,
       out_type = "C", nperm = 999, 
       method = "davies", omnibus = "cauchy")
MiRKAT(y = PFT$Pfev0.5, 
       X = covar2, Ks = Ks, 
       returnKRV = TRUE, returnR2 = TRUE,
       out_type = "C", nperm = 999, 
       method = "davies", omnibus = "cauchy")
# MiRKAT(y = PFT$Pfev0.5, 
#        X = covar3, Ks = Ks, 
#        returnKRV = TRUE, returnR2 = TRUE,
#        out_type = "C", nperm = 999, 
#        method = "davies", omnibus = "cauchy")

MiRKAT(y = PFT$Pfvc, 
        Ks = Ks, 
       returnKRV = TRUE, returnR2 = TRUE,
       out_type = "C", nperm = 999, 
       method = "davies", omnibus = "cauchy")

MiRKAT(y = PFT$Pfvc, 
       X = covar1, Ks = Ks, 
       returnKRV = TRUE, returnR2 = TRUE,
       out_type = "C", nperm = 999, 
       method = "davies", omnibus = "cauchy")

MiRKAT(y = PFT$Pfvc, 
       X = covar2, Ks = Ks, 
       returnKRV = TRUE, returnR2 = TRUE,
       out_type = "C", nperm = 999, 
       method = "davies", omnibus = "cauchy")



shannon.BAL = read.csv("../../Data/Shannon_BAL_decontam.csv")
colnames(shannon.BAL)[1] = "Sample_ID" 
shannon.c1 = shannon.BAL %>% filter(Sample_ID %in% sample.order1)
shannon.c2 = shannon.BAL %>% filter(Sample_ID %in% sample.order2)
shannon.c3 = shannon.BAL %>% filter(Sample_ID %in% sample.order3)
shannon.clusters = rbind.data.frame(shannon.c1, shannon.c2, shannon.c3)
shannon.clusters$Cluster = c(rep(1,dim(shannon.c1)[1]),
                           rep(2,dim(shannon.c2)[1]),
                           rep(3,dim(shannon.c3)[1]))
dim(shannon.clusters)
dim(PFT)
PFT = PFT %>% inner_join(shannon.clusters, by = "Sample_ID")
summary(lm(Shannon~Pfev0.5, PFT))
# write.csv(PFT, "../../Code/publications/stool_BAL/PFT.csv")
summary(lm(Shannon~Pfev0.5, PFT %>% filter(Visit == "BA2")))
summary(lm(Shannon~Pfev0.5, PFT %>% filter(Visit == "BA7")))
summary(lm(Shannon~Pfev0.5, PFT %>% filter(Cluster == 1)))
summary(lm(Shannon~Pfev0.5, PFT %>% filter(Cluster == 2)))
summary(lm(Shannon~Pfev0.5, PFT %>% filter(Cluster == 3))) 
# only in Cluster 3: shannon diversity is sig. associated with FEV0.5.
# the higher the diversity is, the worse the FEV is??? 
# very small sample size 


ggplot(PFT, aes(x = Pfev0.5, y = Shannon, group = Visit)) + 
  geom_point() + 
  geom_smooth(method = lm) +
  theme_classic() 

ggplot(PFT, aes(x = Pfev0.5, y = Shannon, group = Cluster)) + 
  geom_point() + 
  geom_smooth(method = lm) +
  theme_classic() 
## Cluster 3 has the highest shannon diversity and good pFev0.5 
## pFev0.5 is inversely correlated with diversity 


View(PFT)
summary(lm(PFT$Pfev0.5~PFT$Shannon + PFT$Country))

# https://biodatamining.biomedcentral.com/articles/10.1186/s13040-018-0173-9#Sec8
FEV05 = PFT %>% dplyr::select("Sample_ID", "Pfev0.5")
FVC = PFT %>% dplyr::select("Sample_ID", "Pfvc")
OTU = otu_table(BAL.phy_genus_relabun) %>% as.data.frame()
dim(tax)
rownames(OTU) = tax[, "Genus"]
FEV05$Sample_ID == colnames(OTU); FVC$Sample_ID == colnames(OTU)
OTU = OTU[-which(rowSums(OTU) == 0),]
OTU.mean = apply(OTU, 1, mean)
names(OTU.mean) = rownames(OTU)
plot(OTU.mean)
hist(OTU.mean)
sort(OTU.mean, decreasing = T)[1:15]
# Streptococcus                  Haemophilus 
# 32.255436                     6.933735 
# Gemella                  Ramlibacter 
# 5.011736                     4.859394 
# Neisseria               Bradyrhizobium 
# 4.330947                     4.094520 
# Veillonella                 Prevotella_7 
# 3.137104                     2.857502 
# Clostridium sensu stricto 10                       Rothia 
# 2.631579                     2.484340 
# Alloprevotella                Acinetobacter 
# 2.396708                     2.297224 
# Cloacibacterium               Actinobacillus 
# 1.797201                     1.710633 
# Corynebacterium 
# 1.675542 
# number of samples that have 0% for a given taxa 
OTU.zero = as.data.frame(apply(OTU, 1, function(x)sum(x == 0)))
OTU.zero$percentage = 1-OTU.zero$`apply(OTU, 1, function(x) sum(x == 0))`/39
colnames(OTU.zero)[1] = "numZero"
sum(OTU.zero$percentage > 0.1)
OTU.per10 = rownames(OTU.zero)[which(OTU.zero$percentage > .1)]

## calculate spearman correlations between each pair of bacterial genera and pfts
OTU = t(OTU)  # sample as row, taxa as col
FEV05$Sample_ID==rownames(OTU)
FVC$Sample_ID==rownames(OTU)
cor = matrix(NA, nrow = dim(OTU)[2], ncol = 16)
apply(cytokine, 2, function(x)sum(is.na(x))) # too many NA's in TNFa; let's remove it
cytokine = cytokine[,-which(colnames(cytokine) == "TNFac_pgml")]
cytokine = cytokine %>% mutate(Sample_ID = paste0(cytokine$SampleID,
                                       cytokine$Visit))
cytokine = cytokine[match(rownames(OTU), cytokine$Sample_ID),]
colnames(cytokine)
cytokine[,3:8] = as.data.frame(lapply(cytokine[,3:8],function(x)as.numeric(x)))
apply(cytokine, 2, function(x)sum(is.na(x))) # too many NA's in TNFa; let's remove it

for(i in 1:dim(OTU)[2]){
  fev = FEV05$Pfev0.5
  fvc = FVC$Pfvc
  Nec = cytokine$Nec_ugml; IL8 = cytokine$IL1Bc_ngml
  IL6 = cytokine$IL6c_ngml; IL1B = cytokine$IL1Bc_ngml
  cellCount = cytokine$TotalCellCount_cellsml; neutrophil = cytokine$P_neutrophils
  cor[i,1] = cor.test(fev, OTU[,i], method="pearson")$estimate
  cor[i,2] = cor.test(fev, OTU[,i], method="pearson")$p.value
  cor[i,3] = cor.test(fvc, OTU[,i], method="pearson")$estimate
  cor[i,4] = cor.test(fvc, OTU[,i], method="pearson")$p.value
  cor[i,5] = cor.test(Nec, OTU[,i], method="spearman", exact = F)$estimate
  cor[i,6] = cor.test(Nec, OTU[,i], method="spearman", exact = F)$p.value
  cor[i,7] = cor.test(IL8, OTU[,i], method="spearman", exact = F)$estimate
  cor[i,8] = cor.test(IL8, OTU[,i], method="spearman", exact = F)$p.value
  cor[i,9] = cor.test(IL6, OTU[,i], method="spearman", exact = F)$estimate
  cor[i,10] = cor.test(IL6, OTU[,i], method="spearman", exact = F)$p.value
  cor[i,11] = cor.test(IL1B, OTU[,i], method="spearman", exact = F)$estimate
  cor[i,12] = cor.test(IL1B, OTU[,i], method="spearman", exact = F)$p.value
  cor[i,13] = cor.test(cellCount, OTU[,i], method="spearman", exact = F, na.rm = T)$estimate
  cor[i,14] = cor.test(cellCount, OTU[,i], method="spearman", exact = F, na.rm = T)$p.value
  cor[i,15] = cor.test(neutrophil, OTU[,i], method="spearman", exact = F, na.rm = T)$estimate
  cor[i,16] = cor.test(neutrophil, OTU[,i], method="spearman", exact = F, na.rm = T)$p.value
}
colnames(cor) = c("FEV_rho", "FEV_p-value","FVC_rho", "FVC_p-value",
                  "Nec_rho", "Nec_p-value", "IL8_rho", "IL8_p-value", 
                  "IL6_rho", "IL6_p-value", "IL1B_rho", "IL1B_p-value", 
                  "cellCount_rho", "cellCount_p-value", 
                  "neutrophil_rho", "neutrophil_p-value")
rownames(cor) = colnames(OTU)
cor = as.data.frame(cor)
cor$`FEV_p-value-adj` = p.adjust(cor$`FEV_p-value`, method = "fdr")
cor$`FVC_p-value-adj` = p.adjust(cor$`FVC_p-value`, method = "fdr")
cor$`Nec_p-value-adj` = p.adjust(cor$`Nec_p-value`, method = "fdr")
cor$`IL8_p-value-adj` = p.adjust(cor$`IL8_p-value`, method = "fdr")
cor$`IL6_p-value-adj` = p.adjust(cor$`IL6_p-value`, method = "fdr")
cor$`IL1B_p-value-adj` = p.adjust(cor$`IL1B_p-value`, method = "fdr")
cor$`cellCount_p-value-adj` = p.adjust(cor$`cellCount_p-value`, method = "fdr")
cor$`neutrophil_p-value-adj` = p.adjust(cor$`neutrophil_p-value`, method = "fdr")
cor$numSamples = 39 - OTU.zero$numZero
# retain the taxa that are in more than 10% samples 
cor.per10 = cor[which(rownames(cor) %in% OTU.per10),]
# heatmap(cor.per10)
dim(cor.per10)
cor.df = cor.per10[,-25]
cor.df = cor.df[,1:16]
cor.df.sig = cor.df
for(i in seq(2,16,length.out=8)){
  pval = cor.df.sig[,i]
  cor.df.sig[which(pval>=0.1), i-1] = NA
  cor.df.sig[which(pval>=0.1), i] = NA
}

dim(cor.df.sig)  #66 16
# View(cor.df.sig)
cor.df.sig = cor.df.sig[-which(apply(cor.df.sig, 1, function(x)sum(is.na(x))) == 16),]
cor.df.sig1 = cor.df.sig
colnames(cor.df.sig1) = NA
temp1 = data.frame()
for (j in seq(1,15,length.out=8)){
  temp1 = rbind.data.frame(temp1, cor.df.sig1[,c(j,j+1)])
}
cor.df.sig.long = temp1
colnames(cor.df.sig.long) = c("rho", "p-value")
rownames(cor.df.sig.long) = NULL
cor.df.sig.long$Genus = rep(rownames(cor.df.sig), 8)
dim(cor.df.sig.long) #cor.df.sig: 44 16; cor.df.sig.long: 352   3
cor.df.sig.long$Traits = c(rep("FEV",dim(cor.df.sig)[1]), rep("FVC",dim(cor.df.sig)[1]),
                           rep("Nec",dim(cor.df.sig)[1]), rep("IL8",dim(cor.df.sig)[1]),
                           rep("IL6",dim(cor.df.sig)[1]), rep("IL1B",dim(cor.df.sig)[1]),
                           rep("cellCount",dim(cor.df.sig)[1]), rep("neutrophil",dim(cor.df.sig)[1]))

cor.df.sig.long$Genus = factor(cor.df.sig.long$Genus, 
                               levels = rev(rownames(cor.df.sig)))
## no sig for IL8 and IL1B 

# cor.df.sig.long = cor.df.sig.long[-which(cor.df.sig.long$Traits %in% 
#                                     c("IL8", "IL1B")),]
cor.df.sig.long$Traits = factor(cor.df.sig.long$Traits, 
                                levels = c("FEV", "FVC", "Nec", "IL8", "IL6", "IL1B",      
                                  "neutrophil", "cellCount"))
cor.df.sig.long$abs.rho = abs(cor.df.sig.long$rho)
# png("cor.png", width = 8, height = 10, unit = "in", res = 300)
ggplot(cor.df.sig.long, aes(x=Traits, y=Genus, fill=rho)) + 
  # geom_tile() +
  geom_point(aes(size = abs.rho, color = rho)) +
  ylab("") + 
  xlab("") + 
  scale_color_gradient2(low = "#075AFF",
                       mid = "#FFFFCC",
                       high = "#FF0000") +
  theme_classic() + 
  scale_x_discrete("", labels = c(TeX(r'(FEV$_{0.5}$%)'), "FVC%", "Neutrophil elastase", "IL-8",
                                  "IL-6", TeX(r'(IL-1$\beta$)'), "Neutrophil%", "Cell count"), 
                   breaks = c("FEV", "FVC", "Nec", "IL8", "IL6", "IL1B",
                                         "neutrophil", "cellCount")) + 
  guides(size = FALSE, fill = FALSE) + 
  theme(axis.title = element_text(size = 22, face = "bold"),
      axis.text.x = element_text(size = 14, angle = 90, vjust = 0.5),
      axis.text.y = element_text(size = 14),
      plot.title = element_text(size = 20, face = "bold"),
      legend.title = element_blank(),
      legend.text = element_text(size = 20, face = "bold"))
# dev.off()

# glmm + lasso
## need to normalize the data first
## we use clr transformation after adding 1 to each taxa 
## check publications/variable selection/variable selection.R

# library(glmmLasso)
# clr = microbiome::transform(BAL.phy, 'clr')

# save(OTU, cytokine, PFT, covar, 
#      file = "../publications/variable selection/variable_selection.RData")
# 
# OTU = as.data.frame(OTU)
# dim(OTU)
# OTU.adj = OTU + 1
# z_OTU.adj = log(OTU.adj)
# clrx_OTU <- apply(z_OTU.adj, 2, function(x) x - rowMeans(z_OTU.adj))
# 
# OTU = clrx_OTU
# OTU$Sample_ID = rownames(OTU)
# clinic = OTU %>% inner_join(cytokine, by = "Sample_ID") %>%
#   inner_join(PFT, by = "Sample_ID")
# rownames(covar) == OTU$Sample_ID
# clinic = cbind(clinic, covar[,4:9])
# clinic = clinic[,-grep(".y", colnames(clinic), fixed = T)]
# colnames(clinic)[grep(".x", colnames(clinic), fixed = T)] = str_remove(colnames(clinic)[grep(".x", colnames(clinic), fixed = T)],
#                                                                        ".x")
# 
# summary(lm(data=clinic, `Pfev0.5`~Rothia)) # 0.0413 *
# cor(clinic$Rothia, clinic$`Pfev0.5`) # -0.3326298
# cor.test(clinic$Rothia, clinic$`Pfev0.5`) 
# 
# # png("rothia_fev.png", width = 8, height = 8, units = "in", res = 300)
# ggplot(data = clinic, aes(x = Rothia, y = `Pfev0.5`)) + 
#   geom_point() + 
#   stat_smooth(method = "lm", size = 2, 
#               color = wes_palette("FantasticFox1")[3]) + 
#   # stat_poly_eq(parse=T, aes(label = ..eq.label..), formula=y~x) + 
#   ylab(TeX(r'(FEV$_{0.5}$%)', bold=TRUE)) +
#   xlab("Rothia%") + 
#   theme_classic() + 
#   theme(axis.title = element_text(size = 22, face = "bold"),
#         axis.text = element_text(size = 20),
#         plot.title = element_text(size = 20, face = "bold"),
#         legend.title = element_text(size = 20, face = "bold"),
#         legend.text = element_text(size = 20, face = "bold")) + 
#   annotate("text", x = 15, y = 150, size = 8,
#            label = "R = -0.333 \n p-val = 0.041")
# # dev.off()
# 
# 
# summary(lm(data=clinic, `Pfev0.5`~Streptococcus)) # 0.00206 **
# cor(clinic$Streptococcus, clinic$`Pfev0.5`) # -0.484,
# cor.test(clinic$Streptococcus, clinic$`Pfev0.5`) # 0.002062
# 
# # png("strep_fev.png", width = 8, height = 8, units = "in", res = 300)
# ggplot(data = clinic, aes(x = Streptococcus, y = `Pfev0.5`)) + 
#   geom_point() + 
#   stat_smooth(method = "lm", size = 2, 
#               color = wes_palette("FantasticFox1")[1]) + 
#   # stat_poly_eq(parse=T, aes(label = ..eq.label..), formula=y~x) + 
#   ylab(TeX(r'(FEV$_{0.5}$%)', bold=TRUE)) +
#   xlab("Streptococcus%") + 
#   theme_classic() + 
#   theme(axis.title = element_text(size = 22, face = "bold"),
#         axis.text = element_text(size = 20),
#         plot.title = element_text(size = 20, face = "bold"),
#         legend.title = element_text(size = 20, face = "bold"),
#         legend.text = element_text(size = 20, face = "bold")) + 
#   annotate("text", x = 80, y = 150, size = 8,
#            label = "R = -0.484 \n p-val = 0.00206")
# # dev.off()
# 
# summary(lm(data = clinic, `Pfev0.5`~Streptococcus))
# # Streptococcus  -0.29671    0.08933  -3.321  0.00206 ** 
# 
# summary(lm(data = clinic, `Pfvc`~Streptococcus))
# cor(clinic$Streptococcus, clinic$`Pfvc`) # -0.2371121
# cor.test(clinic$Streptococcus, clinic$`Pfvc`) # 0.1518
# 
# # png("strep_fvc.png", width = 8, height = 8, units = "in", res = 300)
# ggplot(data = clinic, aes(x = Streptococcus, y = `Pfvc`)) + 
#   geom_point() + 
#   stat_smooth(method = "lm", size = 2, 
#               color = wes_palette("FantasticFox1")[3]) + 
#   # stat_poly_eq(parse=T, aes(label = ..eq.label..), formula=y~x) + 
#   ylab('FVC%') +
#   xlab("Streptococcus%") + 
#   theme_classic() + 
#   theme(axis.title = element_text(size = 22, face = "bold"),
#         axis.text = element_text(size = 20),
#         plot.title = element_text(size = 20, face = "bold"),
#         legend.title = element_text(size = 20, face = "bold"),
#         legend.text = element_text(size = 20, face = "bold")) + 
#   annotate("text", x = 80, y = 148, size = 8,
#            label = "R = -0.237 \n p-val = 0.152")
# # dev.off()
# 
# 
# 
# # 
# # #### select taxa that are associated with fev0.5
# # rownames(cor.df.sig)[which(!is.na(cor.df.sig$FEV_rho))]
# # clinic$SampleID = factor(clinic$SampleID)
# # bics = rep(NA, 201)
# # ## use the taxa that are significantly associated with FEV0.5 with p-value < 0.05
# # for (k in seq(0,200,1)){
# #   lm = glmmLasso(fix=Pfev0.5~Streptococcus + Gemella + Granulicatella + 
# #                    Rothia + `Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium` + 
# #                    Hephaestia + Novosphingobium +#Actinomyces +  
# #                       Cupriavidus + Bergeyella + #Sphingomonas + 
# #                    Enhydrobacter + PE + 
# #                       PE_ad + AZITH + ERAD + Other, 
# #                   rnd=list(SampleID=~1), data = clinic, 
# #                   lambda = k, 
# #                   family = poisson(link = "log"), 
# #                   switch.NR=T, final.re=T, control = list()) #Sphingomonas +Actinomyces were from spearman
# #   bics[k+1] = lm$bic
# # }
# # names(bics) = seq(0,200,1)
# # which.min(bics) # lambda = 164 ; 58 
# # lm1= glmmLasso(fix=Pfev0.5~Streptococcus + Gemella + Granulicatella + 
# #                  Rothia + `Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium` + 
# #                  Hephaestia + Novosphingobium +#Actinomyces +  
# #                  Cupriavidus + Bergeyella + #Sphingomonas + 
# #                  Enhydrobacter + PE + 
# #                  PE_ad + AZITH + ERAD + Other,  
# #                rnd=list(SampleID=~1), data = clinic, 
# #                lambda = 58, 
# #                family = poisson(link = "log"), 
# #                switch.NR=T, final.re=T, control = list())
# # temp1 = lm1$coefficients
# # temp1
# # # 
# # #  (Intercept) 
# # # 4.660768138 
# # # Streptococcus 
# # # -0.001353555 
# # # Gemella 
# # # -0.001039727 
# # # Granulicatella 
# # # -0.017224803 
# # # Rothia 
# # # -0.005249771 
# # # `Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium` 
# # # -0.025457070 
# # # Hephaestia 
# # # -0.154900441 
# # # Novosphingobium 
# # # 0.009379517 
# # # Cupriavidus 
# # # -0.124976689 
# # # Bergeyella 
# # # 0.037458722 
# # # Enhydrobacter 
# # # 0.031283536 
# # # AZITH 
# # # 0.007494605 
# # 
# # lm1= glmmLasso(fix=Pfev0.5~Streptococcus + Gemella + Granulicatella + 
# #                  Rothia + `Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium` + 
# #                  Hephaestia + Novosphingobium +#Actinomyces +  
# #                  Cupriavidus + Bergeyella + #Sphingomonas + 
# #                  Enhydrobacter + PE + 
# #                  PE_ad + AZITH + ERAD + Other,  
# #                rnd=list(SampleID=~1), data = clinic, 
# #                lambda = 58, 
# #                family = poisson(link = "log"), 
# #                switch.NR=T, final.re=T, control = list())
# # 
# # library(lme4)
# # glmer(Pfev0.5~ Streptococcus + Gemella + Granulicatella +
# #         Rothia + `Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium` +
# #         Hephaestia + Novosphingobium + Cupriavidus + Bergeyella +
# #         Enhydrobacter + AZITH + (1|SampleID), 
# #       family=poisson(link = "log"), data=clinic)
# 
# #### add country as a covariate
# clinic$Group = paste0(clinic$Country, clinic$Visit)
# clinic$Group = factor(clinic$Group, levels = c("USABA2", "USABA7",
#                                                      "AUSBA2", "AUSBA7"))
# clinic$Country = factor(clinic$Country, levels = c("USA", "AUS"))
# clinic$SampleID = factor(clinic$SampleID)
# cor.df.sig.long %>% dplyr::filter(Traits == "FEV") %>% 
#   dplyr::filter(`p-value` < 0.1) %>% 
#   dplyr::select(Genus)
# # Streptococcus, Gemella, Granulicatella, Rothia, Hephaestia, 
# # Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium, Novosphingobium,
# # Cupriavidus, Bergeyella, Enhydrobacter
# 
# bics1 = rep(NA, 201)
# ## use the taxa that are significantly associated with FEV0.5 with p-value < 0.1
# ## lambda is 0 to 201
# ##
# for (k in seq(0,200,1)){
#   lm = glmmwLasso(fix=Pfev0.5~Streptococcus + Gemella + Granulicatella + Rothia +         
#                    Hephaestia + `Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium` + 
#                    Novosphingobium + Cupriavidus + 
#                    Bergeyella + Enhydrobacter + PE + 
#                    PE_ad + ERAD + Other + as.factor(Country), 
#                  rnd=list(SampleID=~1), data = clinic, 
#                  lambda = k, 
#                  family = poisson(link = "log"), 
#                  switch.NR=T, final.re=T, control = list())
#   bics1[k+1] = lm$bic
# }
# names(bics1) = seq(0,200,1) ### names of the bics results 0 to 200
# which.min(bics1) # lambda = 43
# lm2 = glmmLasso(fix=Pfev0.5~Streptococcus + Gemella + Granulicatella + Rothia +         
#                   Hephaestia + `Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium` + 
#                   Novosphingobium + Cupriavidus + 
#                   Bergeyella + Enhydrobacter + PE + 
#                   PE_ad + ERAD + Other + as.factor(Country), 
#                rnd=list(SampleID=~1), data = clinic, 
#                lambda = 43, 
#                family = poisson(link = "log"), 
#                switch.NR=T, final.re=T, control = list())
# # lm2$coefficients
# # 4.682436423 
# # Streptococcus 
# # -0.001108606 
# # Granulicatella 
# # -0.016966109 
# # Rothia 
# # -0.004900300 
# # Hephaestia 
# # -0.152904734 
# # `Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium` 
# # -0.019684970 
# # Novosphingobium 
# # 0.009283743 
# # Cupriavidus 
# # -0.123936706 
# # Bergeyella 
# # 0.034176002 
# # as.factor(Country)AUS 
# # -0.065853643 
# # 
# # ##2. use the taxa that are significantly associated with FEV0.5 with p-value < 0.1
# # bics2 = rep(NA, 201)
# # for (k in seq(0,200,1)){
# #   lm = glmmLasso(fix=Pfev0.5~Streptococcus + Granulicatella + Rothia +         
# #                    Hephaestia + Actinomyces + Novosphingobium +
# #                    Cupriavidus + Bergeyella + Sphingomonas + 
# #                    Prevotella_7 + `Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium` +
# #                    Lactobacillus + Klebsiella + Sphingomonas + Lactococcus + 
# #                    PE + 
# #                    PE_ad + AZITH + ERAD + Other, 
# #                  rnd=list(SampleID=~1), data = clinic, 
# #                  lambda = k, 
# #                  family = poisson(link = "log"), 
# #                  switch.NR=T, final.re=T, control = list())
# #   bics2[k+1] = lm$bic
# # }
# # names(bics2) = seq(0,200,1)
# # which.min(bics2)
# # bics2[164] # lambda = 163
# # lm3 = glmmLasso(fix=Pfev0.5~Streptococcus + Granulicatella + Rothia +         
# #                  Hephaestia + Actinomyces + Novosphingobium +
# #                  Cupriavidus + Bergeyella + Sphingomonas + 
# #                  Prevotella_7 + `Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium` +
# #                  Lactobacillus + Klebsiella + Sphingomonas + Lactococcus + 
# #                  PE + 
# #                  PE_ad + AZITH + ERAD + Other, 
# #                rnd=list(SampleID=~1), data = clinic, 
# #                lambda = 163, 
# #                family = poisson(link = "log"), 
# #                switch.NR=T, final.re=T, control = list())
# # temp3 = lm3$coefficients
# # 
# # # (Intercept) 
# # # 4.664781524 
# # # Streptococcus 
# # # -0.001504412 
# # # Granulicatella 
# # # -0.016960967 
# # # Rothia 
# # # -0.005371521 
# # # Hephaestia 
# # # -0.159784283 
# # # Novosphingobium 
# # # 0.009171818 
# # # Cupriavidus 
# # # -0.110230905 
# # # Bergeyella 
# # # 0.034429351 
# # # `Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium` 
# # # -0.017633875 
# # 
# # 
# # #################################################################
# # #FVC
# # rownames(cor.df.sig)[which(!is.na(cor.df.sig$FVC_rho))]
# # 
