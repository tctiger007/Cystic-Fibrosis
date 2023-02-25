library(phyloseq)
library(ALDEx2)
library(rstudioapi)
library(data.table) 
library(grid)

setwd(dirname(getActiveDocumentContext()$path))
picrust2 = "."
list.files(picrust2, recursive = TRUE)

p2_KO = paste0(picrust2, "/KO_L3.csv")
mapfile = "/Users/wangfei/opt/anaconda3/envs/picrust2/lib/python3.8/site-packages/picrust2/default_files/description_mapfiles"
# list.files(mapfile, recursive = TRUE)
# mapfile_EC = paste0(mapfile, "/ec_level4_info.tsv.gz")
mapfile_KO = paste0(mapfile, "/ko_info.tsv.gz")
# mapfile_PW = paste0(mapfile, "/metacyc_pathways_info.txt.gz")

mapKO = as.data.frame(fread(mapfile_KO, header = FALSE, sep = "\t"))
colnames(mapKO) = c("function","description")

p2KO = as.data.frame(fread(p2_KO))
colnames(p2KO)[1] = "function"

rownames(p2KO) = p2KO$"function"
p2KO = as.matrix(p2KO[,-1])
p2KO = round(p2KO)


BAL.phy = readRDS("../../BAL_phy.rds")
IDs = c("IU120", "IU121", "IU123", "IU124", "IU125", "IU127", "IU128", "IU129", "IU130",
        "IU131", "IU132", "WU119", "WU121", "WU123", "WU124", "M119",  "M120",  "M121",
        "M122",  "M124",  "M126",  "M129", "M130",  "M133",  "M134",  "M135",  "M136",
        "M137",  "M140",  "M141",  "M142", "M143",  "M145",  "M146",  "M148",  "M149",
        "M151",  "M152",  "M153",  "M154", "M155",  "M160",  "M161",  "P216",  "P240")
BAL.phy = subset_samples(BAL.phy, Subject_ID %in% IDs)
meta = meta(BAL.phy)
# write.csv(meta, "meta.csv")
meta$Sequence_ID == colnames(p2KO)
colnames(p2KO) = meta$Sample_ID

T1 = subset_samples(BAL.phy, Visit == "BA2" )
T2 = subset_samples(BAL.phy, Visit == "BA7" )


# Subset-BA2
p2KO1 = p2KO[, sample_names(T1)]

# Subset-BA7
p2KO2 = p2KO[, sample_names(T2)]


#### T1
# KO
set.seed(12345)
system.time({
  aldex2_KO1 = aldex(p2KO1, sample_data(T1)$Country, mc.samples = 500, test = "t",
                     effect = TRUE, denom = "iqlr", verbose = TRUE)
})


head(aldex2_KO1, 10)


#### T2
# KO
set.seed(12345)
system.time({
  aldex2_KO2 = aldex(p2KO2, sample_data(T2)$Country, mc.samples = 500, test = "t",
                     effect = TRUE, denom = "iqlr", verbose = TRUE)
})

hist(aldex2_KO1$effect, breaks = 20, xlab = "effect size", col = "yellow", main = "KO")
hist(aldex2_KO2$effect, breaks = 20, xlab = "effect size", col = "yellow", main = "KO")

par(mfrow = c(3,2))

plot(aldex2_KO1$effect, aldex2_KO1$we.ep, log = "y", cex = 0.4, col = "blue", pch = 19,
     xlab = "Effect size", ylab = "P value", main = "(KO) Effect size plot")
points(aldex2_KO1$effect, aldex2_KO1$we.eBH, cex = 0.5, col = "red", pch = 19)
abline(h = 0.05, lty = 2, col = "grey")
legend("bottom", legend = c("P value", "BH-adjusted"), pch = 19, col = c("blue", "red"))

plot(aldex2_KO1$diff.btw, aldex2_KO1$we.ep, log = "y", cex = 0.4, col = "blue", pch = 19,
     xlab = "Difference", ylab = "P value", main = "(KO) Volcano plot")
points(aldex2_KO1$diff.btw, aldex2_KO1$we.eBH, cex = 0.5, col = "red", pch = 19)
abline(h = 0.05, lty = 2, col = "grey")

# https://usfomicshub.github.io/Workshops/Microbiome_Workshop_Materials/microbiome_workshop_demos/day3/Ranalysis/PartII/

library(tidyverse)
library(ALDEx2)
library(readr)
library(dplyr)
#wide file?

dfk=read_tsv("picrust2_out/KO_metagenome_out/pred_metagenome_unstrat.tsv.gz")
colnames(dfk) == c("function", meta$Sequence_ID)
colnames(dfk) = c("function", meta$Sample_ID)
dfk=as.data.frame(dfk)
rownames(dfk)=dfk$`function`
dfk=dfk[,-which(colnames(dfk)=="function")]


library(phyloseq)
library(microbiome)
mdata = read.table("metadata.txt", header = T)
kotu=otu_table(dfk,taxa_are_rows=TRUE)
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
# PC1        PC2        PC3        PC4        PC5        PC6        PC7      PC8
# 0.24750628 0.07272967 0.05248926 0.04416388 0.04069374 0.03516856 0.03065760 0.02685558 


p_kg=plot_ordination(kegg,ord_kegg,color="Country", shape = "Visit") #,shape="Genotype")
p_kg1=p_kg+geom_polygon(aes(fill=Country))
p_kg1


#aldex2_EC1 KO1 PW1 EC2 KO2 PW2
aldex2_KO1_sig = aldex2_KO1[aldex2_KO1$we.eBH< 0.05, ]
aldex2_KO1_sig$path=rownames(aldex2_KO1_sig)
aldex2_KO1_sig2=arrange(aldex2_KO1_sig,desc(abs(effect)))
# aldex2_KO1_sig2$path = mapKO$description[match(rownames(aldex2_KO1_sig2), mapKO$`function`)]
aldex2_KO1_sig3=aldex2_KO1_sig2[1:20,] 
print(nrow(aldex2_KO1_sig)/nrow(aldex2_KO1))  #0.09489051

# write.csv(aldex2_KO1_sig, "KO1_L3_sig.csv")
## no DE KO2
aldex2_KO2_sig = aldex2_KO2[aldex2_KO2$we.eBH< 0.05, ]


text_left <- textGrob("Increase in USA", gp=gpar(fontsize=18))
text_right <- textGrob("Increase in AUS", gp=gpar(fontsize=18))
text_center = textGrob("Effect", gp=gpar(fontsize=18))


aldex2_KO1_sig_copy = aldex2_KO1_sig
df <- aldex2_KO1_sig_copy[order(aldex2_KO1_sig_copy$effect, decreasing = T),]
df$path = factor(df$path, df$path)

# relevel = sort(aldex2_KO1_sig_copy$effect, decreasing = T)
## plot
wesanderson::wes_palette(name = "FantasticFox1", n = 5)[5]
df = df %>% mutate(color = ifelse(effect >0, 
                                  wesanderson::wes_palette(name = "FantasticFox1", n = 5)[5], 
                                  wesanderson::wes_palette(name = "FantasticFox1", n = 5)[3]))
png("./KO_L3_T1.png", width = 12, height = 8, units = "in", res = 300)
ggplot(df, aes(y = path, x = effect)) +
  geom_point(aes(size = -log(we.eBH)), color = df$color) + 
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
  coord_cartesian(ylim = c(0,26),  clip = 'off')
dev.off()












