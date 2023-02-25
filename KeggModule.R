library(phyloseq)
library(ALDEx2)
library(rstudioapi)
library(data.table) 
library(grid)
library(KEGGREST)
library(readxl)
library(readr)
library(tibble)
library(microbiome)

setwd(dirname(getActiveDocumentContext()$path))
ARG_module = read_excel("ARG_module.xlsx")

listDatabases()
mod = keggList("module")
head(mod)

modules = ARG_module$Module
KEGGs = list()
for(i in 1:length(modules)){
  mod = modules[i]
  KEGGs[[i]] = keggGet(mod)[[1]]$DEFINITION
}
names(KEGGs) = modules

KEGGs1 <- as.data.frame(do.call(rbind, KEGGs))
# write.csv(KEGGs1, "KEGGs1.csv")
KO = read_tsv("picrust2_out/KO_metagenome_out/pred_metagenome_unstrat.tsv")
KO = column_to_rownames(KO, var = "function")

KEGGs1 <- read.csv("KEGGs1.csv")
rownames(KEGGs1) = KEGGs1$X
KEGGs1 = KEGGs1[,-1]
colnames(KEGGs1) = c("KEGG1","KEGG2","KEGG3","KEGG4","KEGG5","KEGG6",
                     "KEGG7","KEGG8","KEGG9","KEGG10","KEGG11")
KEGGs1.t = t(KEGGs1)
mod.sum = matrix(NA, nrow = 25, ncol = dim(KO)[2])
colnames(mod.sum) = colnames(KO)

for (i in 1:dim(KEGGs1.t)[2]){
  KEGG = KEGGs1.t[,i]
  KEGG = KEGG[KEGG != ""]
  mod.temp = matrix(0, nrow = 1, ncol = dim(KO)[2])
  colnames(mod.temp) = colnames(KO)
  for (j in 1:length(KEGG)){
    if(is.na(match(KEGG[j], rownames(KO)))){next}
    mod.temp = mod.temp + KO[match(KEGG[j], rownames(KO)),]
  }
  mod.temp = t(mod.temp)
  mod.sum[i,] = mod.temp[,1]
}
rownames(mod.sum) = modules
# saveRDS(mod.sum, "ARG_module.rds")

### DE analysis
ARG_module1 = readRDS("ARG_module.rds")

BAL.phy = readRDS("../../BAL_phy.rds")
IDs = c("IU120", "IU121", "IU123", "IU124", "IU125", "IU127", "IU128", "IU129", "IU130",
        "IU131", "IU132", "WU119", "WU121", "WU123", "WU124", "M119",  "M120",  "M121",
        "M122",  "M124",  "M126",  "M129", "M130",  "M133",  "M134",  "M135",  "M136",
        "M137",  "M140",  "M141",  "M142", "M143",  "M145",  "M146",  "M148",  "M149",
        "M151",  "M152",  "M153",  "M154", "M155",  "M160",  "M161",  "P216",  "P240")
BAL.phy = subset_samples(BAL.phy, Subject_ID %in% IDs)
meta = meta(BAL.phy)
genus = aggregate_taxa(BAL.phy, "Genus", verbose = FALSE)
tax = tax_table(genus)


rela  = transform_sample_counts(genus, function(x) x / sum(x) )
otu = otu_table(rela)

mean(otu[which(rownames(otu) == "Pseudomonas"),])
mean(otu[which(rownames(otu) == "Staphylococcus"),])

meta$Sequence_ID == colnames(ARG_module1)
colnames(ARG_module1) = meta$Sample_ID

T1 = subset_samples(BAL.phy, Visit == "BA2" )
T2 = subset_samples(BAL.phy, Visit == "BA7" )

ARG_module1 = round(ARG_module1)

# Subset-BA2
ARG_module.T1 = ARG_module1[, sample_names(T1)]

# Subset-BA7
ARG_module.T2 = ARG_module1[, sample_names(T2)]


#### T1
# KO
set.seed(12345)
system.time({
  aldex2_KO1 = aldex(ARG_module.T1, sample_data(T1)$Country, mc.samples = 500, test = "t",
                     effect = TRUE, denom = "iqlr", verbose = TRUE)
})


head(aldex2_KO1, 10)


#### T2
# KO
set.seed(12345)
system.time({
  aldex2_KO2 = aldex(ARG_module.T2, sample_data(T2)$Country, mc.samples = 500, test = "t",
                     effect = TRUE, denom = "iqlr", verbose = TRUE)
})

hist(aldex2_KO1$effect, breaks = 20, xlab = "effect size", col = "yellow", main = "KO")
hist(aldex2_KO2$effect, breaks = 20, xlab = "effect size", col = "yellow", main = "KO")

par(mfrow = c(3,2))

plot(aldex2_KO1$effect, aldex2_KO1$we.ep, log = "y", cex = 0.4, col = "blue", pch = 19,
     xlab = "Effect size", ylab = "P value", main = "(KO) Effect size plot")
points(aldex2_KO1$effect, aldex2_KO1$wi.eBH, cex = 0.5, col = "red", pch = 19)
abline(h = 0.05, lty = 2, col = "grey")
legend("bottom", legend = c("P value", "BH-adjusted"), pch = 19, col = c("blue", "red"))

plot(aldex2_KO1$diff.btw, aldex2_KO1$we.ep, log = "y", cex = 0.4, col = "blue", pch = 19,
     xlab = "Difference", ylab = "P value", main = "(KO) Volcano plot")
points(aldex2_KO1$diff.btw, aldex2_KO1$wi.eBH, cex = 0.5, col = "red", pch = 19)
abline(h = 0.05, lty = 2, col = "grey")

# reference https://usfomicshub.github.io/Workshops/Microbiome_Workshop_Materials/microbiome_workshop_demos/day3/Ranalysis/PartII/

library(tidyverse)
library(ALDEx2)
library(readr)
library(dplyr)
library(phyloseq)
library(microbiome)

dfk=ARG_module1
kotu=otu_table(dfk,taxa_are_rows=TRUE)
mdata = meta
md=as.data.frame(mdata)
rownames(md)=mdata$Sample_ID
md=md[,-which(colnames(md)=="Sample_ID")]
md=sample_data(md)

kegg=phyloseq(kotu,md)
#kegg is the kegg table
#pw is the pathway table (much higher level)
kegg_clr=microbiome::transform(kegg,"clr")

#rda is one of many ordination methods
ord_kegg=phyloseq::ordinate(kegg_clr,"RDA")
#get the top eigenvalues of the first few PC axes

sapply(ord_kegg$CA$eig[1:8], function(x) x / sum(ord_kegg$CA$eig))
# PC1        PC2        PC3        PC4        PC5        PC6        PC7 
# 0.23979809 0.19697412 0.10678005 0.07123763 0.05537256 0.05132079 0.04539295 
# PC8 
# 0.04003278 

p_kg=plot_ordination(kegg,ord_kegg,color="Country", shape = "Visit") #,shape="Genotype")
p_kg1=p_kg+geom_polygon(aes(fill=Country))
p_kg1

#aldex2_EC1 KO1 PW1 EC2 KO2 PW2
aldex2_KO1_sig = aldex2_KO1[aldex2_KO1$wi.eBH< 0.1, ]
aldex2_KO1_sig$path=rownames(aldex2_KO1_sig)
aldex2_KO1_sig2=arrange(aldex2_KO1_sig,desc(abs(effect)))

## no DE KO2
aldex2_KO2_sig = aldex2_KO2[aldex2_KO2$wi.eBH< 0.1, ]


text_left <- textGrob("Increase in USA", gp=gpar(fontsize=18))
text_right <- textGrob("Increase in AUS", gp=gpar(fontsize=18))
text_center = textGrob("Effect", gp=gpar(fontsize=18))


aldex2_KO1_sig_copy = aldex2_KO1_sig
df <- aldex2_KO1_sig_copy[order(aldex2_KO1_sig_copy$effect, decreasing = T),]
df$path = factor(df$path, df$path)

## plot
wesanderson::wes_palette(name = "FantasticFox1", n = 5)[5]
df = df %>% mutate(color = ifelse(effect >0, 
                                  wesanderson::wes_palette(name = "FantasticFox1", n = 5)[5], 
                                  wesanderson::wes_palette(name = "FantasticFox1", n = 5)[3]))

# png("./ARG_module.png", width = 12, height = 8.5, units = "in", res = 300)
ggplot(df, aes(y = path, x = effect)) +
  geom_point(aes(size = -log(wi.eBH)), color = df$color) + 
  # xlab("Effect") +
  ylab("") + xlab("") +
  scale_size_continuous("Adjusted p-value", 
                        breaks = c(3.5, 4),
                        labels = c("0.03", "0.02")) + 
  scale_color_manual("Effect", values = c("#B40F20", "#46ACC8")) +
  # breaks = c("#B40F20", "#46ACC8"),
  # labels = c("Increase in AUS",
  #            "Increase in USA")) +
  theme_classic() + 
  xlim(c(-1.2,1.2))+
  theme(legend.position = "bottom",
        legend.title = element_text(size = 17),
        legend.text = element_text(size = 16),
        plot.margin = unit(c(1.5,1,2,1), "lines"),
        axis.title.x = element_text(vjust=0),
        axis.title = element_text(size = 17, face = "bold"),
        axis.text = element_text(size = 17,face="bold")) +
  annotation_custom(text_left,xmin=-1.2,xmax=-.4,ymin=-2.7,ymax=-2.7) +
  annotation_custom(text_right,xmin=.4,xmax=1.2,ymin=-2.7,ymax=-2.7)  +
  annotation_custom(text_center,xmin=-.1,xmax=.1,ymin=-2.7,ymax=-2.7)  +
  geom_segment(aes(x =-0.2, xend = -1.1, y = -2, yend = -2),
               arrow = arrow(length = unit(0.2, "cm")), color = "black") +
  geom_segment(aes(x =0.2, xend = 1.1, y = -2, yend = -2),
               arrow = arrow(length = unit(0.2, "cm")), color = "black") +
  coord_cartesian(ylim = c(0,3),  clip = 'off')
# dev.off()


module_converter = read_excel("ARG_module.xlsx")
module_converter$Module == rownames(ARG_module1)
USABA2 = meta(T1) %>% filter(Country == "USA") %>% dplyr::select(Sample_ID)
USABA7 = meta(T2) %>% filter(Country == "USA") %>% dplyr::select(Sample_ID)
AUSBA2 = meta(T1) %>% filter(Country == "AUS") %>% dplyr::select(Sample_ID)
AUSBA7 = meta(T2) %>% filter(Country == "AUS") %>% dplyr::select(Sample_ID)
sampleID = c(sort(USABA2$Sample_ID), 
             sort(AUSBA2$Sample_ID), 
             sort(USABA7$Sample_ID), 
             sort(AUSBA7$Sample_ID))
ARG_module_Rela = ARG_module1/rowSums(ARG_module1)
ARG_module_Rela = ARG_module_Rela[,sampleID]


library(pheatmap)       
orderID <- rownames(aldex2_KO1)[order(aldex2_KO1$wi.eBH, decreasing = F)]
ARG_module_Rela = ARG_module_Rela[orderID,]
rownames(ARG_module_Rela) = paste0(rownames(ARG_module_Rela), ": ",
                                   module_converter$Name[match(rownames(ARG_module_Rela), module_converter$Module)])
rownames(ARG_module_Rela)[c(1:3)] = paste0(rownames(ARG_module_Rela)[c(1:3)], c("(*)","(\U00B7)","(\U00B7)"))

# rownames(ARG_module_Rela) = paste0(module_converter$Module, ": ",
                                   # module_converter$Name)
my_sample_col = data.frame(Group = c(rep("USABA2", dim(USABA2)[1]),
                                       rep("AUSBA2", dim(AUSBA2)[1]),
                                       rep("USABA7", dim(USABA7)[1]),
                                       rep("AUSBA7", dim(AUSBA7)[1])))

ann_colors = list(
  Group = c(USABA2 = "#DBEDD5", AUSBA2 = "#E0989D",
            USABA7 = "#F8D9BA", AUSBA7 = "#A8C2DA"))
rownames(my_sample_col) = sampleID
newnames <- lapply(
  rownames(ARG_module_Rela)[1:3],
  function(x) bquote(bold(.(x))))
newnames = c(newnames, as.list(rownames(ARG_module_Rela)[4:25]))


dev.off()
png("../Figures/keggmodule.png",width = 16,height = 6.5,units = "in",res=300)
pheatmap(ARG_module_Rela, 
         cluster_rows = F,
         cluster_cols = F,
         border_color = NA,
         color=colorRampPalette(c("#000033", "#66CCFF"))(50),
         annotation_col = my_sample_col,
         annotation_colors = ann_colors,
         fontsize_col = 7,
         fontsize_row = 11,
         labels_row = as.expression(newnames))
dev.off()

                            