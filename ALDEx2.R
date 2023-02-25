library(phyloseq)
library(ALDEx2)
library(rstudioapi)
library(data.table) 
library(grid)
library(microbiome)

setwd(dirname(getActiveDocumentContext()$path))
picrust2 = "picrust2_out"
list.files(picrust2, recursive = TRUE)

p2_EC = paste0(picrust2, "/EC_metagenome_out/pred_metagenome_unstrat.tsv.gz")
p2_KO = paste0(picrust2, "/KO_metagenome_out/pred_metagenome_unstrat.tsv.gz")
p2_PW = paste0(picrust2, "/pathways_out/path_abun_unstrat.tsv.gz")

# mapfile = "/Users/wangfei/opt/anaconda3/envs/picrust2/lib/python3.8/site-packages/picrust2/default_files/description_mapfiles"
mapfile = "~/Documents/picrust2-2.5.0/picrust2/default_files/description_mapfiles"

list.files(mapfile, recursive = TRUE)
mapfile_EC = paste0(mapfile, "/ec_level4_info.tsv.gz")
mapfile_KO = paste0(mapfile, "/ko_info.tsv.gz")
mapfile_PW = paste0(mapfile, "/metacyc_pathways_info.txt.gz")

mapEC = as.data.frame(fread(mapfile_EC, header = FALSE))
colnames(mapEC) = c("function","description")
mapKO = as.data.frame(fread(mapfile_KO, header = FALSE, sep = "\t"))
colnames(mapKO) = c("function","description")
mapPW = as.data.frame(fread(mapfile_PW, header = FALSE))
colnames(mapPW) = c("pathway","description")

p2EC = as.data.frame(fread(p2_EC))
rownames(p2EC) = p2EC$"function"
p2EC = as.matrix(p2EC[,-1])
p2EC = round(p2EC)

p2KO = as.data.frame(fread(p2_KO))
rownames(p2KO) = p2KO$"function"
p2KO = as.matrix(p2KO[,-1])
p2KO = round(p2KO)

p2PW = as.data.frame(fread(p2_PW))
rownames(p2PW) = p2PW$"pathway"
p2PW = as.matrix(p2PW[,-1])
p2PW = round(p2PW)


BAL.phy = readRDS("../Data/BAL_phy.rds")
IDs = c("IU120", "IU121", "IU123", "IU124", "IU125", "IU127", "IU128", "IU129", "IU130",
        "IU131", "IU132", "WU119", "WU121", "WU123", "WU124", "M119",  "M120",  "M121",
        "M122",  "M124",  "M126",  "M129", "M130",  "M133",  "M134",  "M135",  "M136",
        "M137",  "M140",  "M141",  "M142", "M143",  "M145",  "M146",  "M148",  "M149",
        "M151",  "M152",  "M153",  "M154", "M155",  "M160",  "M161",  "P216",  "P240")
BAL.phy = subset_samples(BAL.phy, Subject_ID %in% IDs)
meta = meta(BAL.phy)
meta$Sequence_ID == colnames(p2EC); meta$Sequence_ID == colnames(p2KO); meta$Sequence_ID == colnames(p2PW)
colnames(p2EC) = meta$Sample_ID
colnames(p2KO) = meta$Sample_ID
colnames(p2PW) = meta$Sample_ID

T1 = subset_samples(BAL.phy, Visit == "BA2" )
T2 = subset_samples(BAL.phy, Visit == "BA7" )


# Subset-BA2
p2EC1 = p2EC[, sample_names(T1)]
p2KO1 = p2KO[, sample_names(T1)]
p2PW1 = p2PW[, sample_names(T1)]

# Subset-BA7
p2EC2 = p2EC[, sample_names(T2)]
p2KO2 = p2KO[, sample_names(T2)]
p2PW2 = p2PW[, sample_names(T2)]


#### T1
# EC
set.seed(12345)
system.time({
  aldex2_EC1 = aldex(p2EC1, sample_data(T1)$Country, mc.samples = 500, test = "t",
                     effect = TRUE, denom = "iqlr", verbose = TRUE)
})


# KO
set.seed(12345)
system.time({
  aldex2_KO1 = aldex(p2KO1, sample_data(T1)$Country, mc.samples = 500, test = "t",
                     effect = TRUE, denom = "iqlr", verbose = TRUE)
})

# Pathway
set.seed(12345)
system.time({
  aldex2_PW1 = aldex(p2PW1, sample_data(T1)$Country, mc.samples = 500, test = "t",
                     effect = TRUE, denom = "iqlr", verbose = TRUE)
})

head(aldex2_EC1, 10)



#### T2
# EC
set.seed(12345)
system.time({
  aldex2_EC2 = aldex(p2EC2, sample_data(T2)$Country, mc.samples = 500, test = "t",
                     effect = TRUE, denom = "iqlr", verbose = TRUE)
})


# KO
set.seed(12345)
system.time({
  aldex2_KO2 = aldex(p2KO2, sample_data(T2)$Country, mc.samples = 500, test = "t",
                     effect = TRUE, denom = "iqlr", verbose = TRUE)
})

# Pathway
set.seed(12345)
system.time({
  aldex2_PW2 = aldex(p2PW2, sample_data(T2)$Country, mc.samples = 500, test = "t",
                     effect = TRUE, denom = "iqlr", verbose = TRUE)
})


par(mfrow = c(2,2))
hist(aldex2_EC1$effect, breaks = 20, xlab = "effect size", col = "yellow", main = "EC")
hist(aldex2_KO1$effect, breaks = 20, xlab = "effect size", col = "yellow", main = "KO")
hist(aldex2_PW1$effect, breaks = 20, xlab = "effect size", col = "yellow", main = "Pathway")


hist(aldex2_EC2$effect, breaks = 20, xlab = "effect size", col = "yellow", main = "EC")
hist(aldex2_KO2$effect, breaks = 20, xlab = "effect size", col = "yellow", main = "KO")
hist(aldex2_PW2$effect, breaks = 20, xlab = "effect size", col = "yellow", main = "Pathway")

par(mfrow = c(3,2))
plot(aldex2_EC1$effect, aldex2_EC1$we.ep, log = "y", cex = 0.4, col = "blue", pch = 19,
     xlab = "Effect size", ylab = "P value", main = "(EC) Effect size plot")
points(aldex2_EC1$effect, aldex2_EC1$we.eBH, cex = 0.5, col = "red", pch = 19)
abline(h = 0.05, lty = 2, col = "grey")
legend("bottom", legend = c("P value", "BH-adjusted"), pch = 19, col = c("blue", "red"))

plot(aldex2_EC1$diff.btw, aldex2_EC1$we.ep, log = "y", cex = 0.4, col = "blue", pch = 19,
     xlab = "Difference", ylab = "P value", main = "(EC) Volcano plot")
points(aldex2_EC1$diff.btw, aldex2_EC1$we.eBH, cex = 0.5, col = "red", pch = 19)
abline(h = 0.05, lty = 2, col = "grey")

plot(aldex2_KO1$effect, aldex2_KO1$we.ep, log = "y", cex = 0.4, col = "blue", pch = 19,
     xlab = "Effect size", ylab = "P value", main = "(KO) Effect size plot")
points(aldex2_KO1$effect, aldex2_KO1$we.eBH, cex = 0.5, col = "red", pch = 19)
abline(h = 0.05, lty = 2, col = "grey")
legend("bottom", legend = c("P value", "BH-adjusted"), pch = 19, col = c("blue", "red"))

plot(aldex2_KO1$diff.btw, aldex2_KO1$we.ep, log = "y", cex = 0.4, col = "blue", pch = 19,
     xlab = "Difference", ylab = "P value", main = "(KO) Volcano plot")
points(aldex2_KO1$diff.btw, aldex2_KO1$we.eBH, cex = 0.5, col = "red", pch = 19)
abline(h = 0.05, lty = 2, col = "grey")

plot(aldex2_PW1$effect, aldex2_PW1$we.ep, log = "y", cex = 0.4, col = "blue", pch = 19,
     xlab = "Effect size", ylab = "P value", main = "(PW) Effect size plot")
points(aldex2_PW1$effect, aldex2_PW1$we.eBH, cex = 0.5, col = "red", pch = 19)
abline(h = 0.05, lty = 2, col = "grey")
legend("bottom", legend = c("P value", "BH-adjusted"), pch = 19, col = c("blue", "red"))

plot(aldex2_PW1$diff.btw, aldex2_PW1$we.ep, log = "y", cex = 0.4, col = "blue", pch = 19,
     xlab = "Difference", ylab = "P value", main = "(PW) Volcano plot")
points(aldex2_PW1$diff.btw, aldex2_PW1$we.eBH, cex = 0.5, col = "red", pch = 19)
abline(h = 0.05, lty = 2, col = "grey")

# reference https://usfomicshub.github.io/Workshops/Microbiome_Workshop_Materials/microbiome_workshop_demos/day3/Ranalysis/PartII/

library(tidyverse)
library(ALDEx2)
library(readr)
library(dplyr)
#wide file?
#this is taken from the prep_input code
df=read_tsv("picrust2_out/pathways_out/path_abun_unstrat.tsv.gz")
colnames(df)[-1] == meta$Sequence_ID
colnames(df) = c("pathway", meta$Sample_ID)
dfus=as.data.frame(df)
# rownames(dfs)=paste(df$pathway,df$sequence,sep="|")
rownames(dfus) = df$pathway
dfus=dfus[,-which(colnames(dfus) == "pathway")]


#nothing in the pathway, look in the kegg?
dfk=read_tsv("picrust2_out/KO_metagenome_out/pred_metagenome_unstrat.tsv.gz")
colnames(dfk) == c("function", meta$Sequence_ID)
colnames(dfk) = c("function", meta$Sample_ID)
dfk=as.data.frame(dfk)
rownames(dfk)=dfk$`function`
dfk=dfk[,-which(colnames(dfk)=="function")]
rownames(dfout)=dfout$pathway
dfout=dfout[,-which(colnames(dfout) == "pathway")]
# df_tmp= df[, -which(colnames(df) == "sequence")]
# df_abun=aggregate(.~pathway,data=df_tmp,FUN=sum)
# rownames(df_abun)=df_abun$pathway
# df_abun=df_abun[,-which(colnames(df_abun)=="pathway")]
#convert to relative abundance
# df_rabun=data.frame(sweep(df_abun,2,colSums(df_abun),"/"),check.names=F)
# asin transformation
# asinT=function(x) asin(sqrt(x))
# df_rabun2=df_rabun %>% mutate_all(asinT)
#finally, round them so they
# df_rrabun=round(df_rabun)

library(phyloseq)
library(microbiome)
mdata = read.table("metadata.txt", header = T)
kotu=otu_table(dfk,taxa_are_rows=TRUE)
md=as.data.frame(mdata)
rownames(md)=mdata$Sample_ID
md=md[,-which(colnames(md)=="Sample_ID")]
md=sample_data(md)

kegg=phyloseq(kotu,md)
pwo=otu_table(dfout,taxa_are_rows=TRUE)
pw=phyloseq(pwo,md)
#kegg is the kegg table
#pw is the pathway table (much higher level)
pw_clr=microbiome::transform(pw,"clr")
kegg_clr=microbiome::transform(kegg,"clr")

#rda is one of many ordination methods
ord_pw=phyloseq::ordinate(pw_clr,"RDA")
ord_kegg=phyloseq::ordinate(kegg_clr,"RDA")
#get the top eigenvalues of the first few PC axes
head(ord_pw$CA$eig)

sapply(ord_pw$CA$eig[1:8], function(x) x / sum(ord_pw$CA$eig))
# PC1        PC2        PC3        PC4        PC5        PC6        PC7 
# 0.23995713 0.10488448 0.04792208 0.04648600 0.04191258 0.03236779 0.02957236 
# PC8 
# 0.02831186 
sapply(ord_kegg$CA$eig[1:8], function(x) x / sum(ord_kegg$CA$eig))
# PC1        PC2        PC3        PC4        PC5        PC6        PC7 
# 0.24750628 0.07272967 0.05248926 0.04416388 0.04069374 0.03516856 0.03065760 
# PC8 
# 0.02685558 

# 3-4 axes for pathway and kegg
p_pw=plot_ordination(pw,ord_pw,color="Country", shape = "Visit") #,shape="Genotype")
p_pw
p_pw1=p_pw+geom_polygon(aes(fill=Country))
p_pw1

# p_kg=plot_ordination(kegg,ord_kegg,color="Country", shape = "Visit") #,shape="Genotype")
# p_kg1=p_kg+geom_polygon(aes(fill=Country))
# p_kg1

# conds=factor(mdata$Group, levels = c("USABA2", "AUSBA2",
#                                      "USABA7", "AUSBA7"))
# #kegg
# x.kf=aldex(round(dfk),conds,effects=T,denom="all")
# 
# #unstratified pathway
# x.pwu=aldex(round(dfus),conds,effects=TRUE,denom="all")
# 
# #stratified pathway
# x.pws=aldex(round(dfs),conds,effects=TRUE)


#aldex2_EC1 KO1 PW1 EC2 KO2 PW2
aldex2_EC1_sig = aldex2_EC1[aldex2_EC1$we.eBH< 0.05, ]
aldex2_EC1_sig$path=rownames(aldex2_EC1_sig)
aldex2_EC1_sig2=arrange(aldex2_EC1_sig,desc(abs(effect)))
aldex2_EC1_sig2$path = mapEC$description[match(rownames(aldex2_EC1_sig2), mapEC$`function`)]
aldex2_EC1_sig3=aldex2_EC1_sig2[1:25,]
print(nrow(aldex2_EC1_sig)/nrow(aldex2_EC1))  #0.06811731

aldex2_KO1_sig = aldex2_KO1[aldex2_KO1$we.eBH< 0.05, ]
aldex2_KO1_sig$path=rownames(aldex2_KO1_sig)
aldex2_KO1_sig2=arrange(aldex2_KO1_sig,desc(abs(effect)))
aldex2_KO1_sig2$path = mapKO$description[match(rownames(aldex2_KO1_sig2), mapKO$`function`)]
aldex2_KO1_sig3=aldex2_KO1_sig2[1:25,] 
print(nrow(aldex2_KO1_sig)/nrow(aldex2_KO1))  #0.08619409

aldex2_PW1_sig = aldex2_PW1[aldex2_PW1$we.eBH< 0.05, ]
aldex2_PW1_sig$path=rownames(aldex2_PW1_sig)
aldex2_PW1_sig2=arrange(aldex2_PW1_sig,desc(abs(effect)))
aldex2_PW1_sig2$path = mapPW$description[match(rownames(aldex2_PW1_sig2), mapPW$pathway)]
aldex2_PW1_sig3=aldex2_PW1_sig2[1:25,] 
print(nrow(aldex2_PW1_sig)/nrow(aldex2_PW1))  

## no DE EC2
aldex2_EC2_sig = aldex2_EC2[aldex2_EC2$we.eBH< 0.05, ]

## no DE KO2
aldex2_KO2_sig = aldex2_KO2[aldex2_KO2$we.eBH< 0.05, ]

## no DE PW2
aldex2_PW2_sig = aldex2_PW2[aldex2_PW2$we.eBH< 0.05, ]


text_left <- textGrob("Increase in USA", gp=gpar(fontsize=18))
text_right <- textGrob("Increase in AUS", gp=gpar(fontsize=18))
text_center <- textGrob("Effect", gp=gpar(fontsize=18))


aldex2_EC1_sig3 <- aldex2_EC1_sig3[order(aldex2_EC1_sig3$effect, decreasing = T),]
aldex2_EC1_sig3$path = factor(aldex2_EC1_sig3$path, aldex2_EC1_sig3$path)
aldex2_EC1_sig3 = aldex2_EC1_sig3 %>% mutate(color = ifelse(effect >0, "#B40F20","#46ACC8"))
    
aldex2_KO1_sig3 <- aldex2_KO1_sig3[order(aldex2_KO1_sig3$effect, decreasing = T),]
aldex2_KO1_sig3$path = factor(aldex2_KO1_sig3$path, aldex2_KO1_sig3$path)
aldex2_KO1_sig3 = aldex2_KO1_sig3 %>% mutate(color = ifelse(effect >0, "#B40F20","#46ACC8"))

aldex2_PW1_sig3 <- aldex2_PW1_sig3[order(aldex2_PW1_sig3$effect, decreasing = T),]
aldex2_PW1_sig3$path = factor(aldex2_PW1_sig3$path, aldex2_PW1_sig3$path)
aldex2_PW1_sig3 = aldex2_PW1_sig3 %>% mutate(color = ifelse(effect >0, "#B40F20","#46ACC8"))

png("./EC1.png", width = 15, height = 8.5, units = "in", res = 300)
ggplot(aldex2_EC1_sig3, aes(y = path, x = effect)) +
  geom_point(aes(size = -log(we.eBH)), color = aldex2_EC1_sig3$color) + 
  xlab("") +
  ylab("") +
  scale_size_continuous("Adjusted p-value", 
                        breaks = c(3.5, 4, 4.6),
                        labels = c("0.03", "0.02", "0.01")) + 
  scale_color_manual("Effect", values = c("#B40F20", "#46ACC8"))+
                     # breaks = c("#B40F20", "#46ACC8"),
                     # labels = c("Increase in AUS", 
                     #            "Increase in USA")) +
  theme_classic() + 
  xlim(c(-1.2,1.2)) +
  theme(legend.position = "bottom",
        legend.title = element_text(size = 17),
        legend.text = element_text(size = 14.5),
        plot.margin = unit(c(1.5,1,2,1), "lines"),
        axis.title.x = element_text(vjust=0),
        axis.title = element_text(size = 17, face = "bold"),
        axis.text = element_text(size = 17, face = "bold")) +
  annotation_custom(text_left,xmin=-1.2,xmax=-.4,ymin=-2.7,ymax=-2.7) +
  annotation_custom(text_right,xmin=.4,xmax=1.2,ymin=-2.7,ymax=-2.7)  +
  # annotation_custom(text_center,xmin=-.1,xmax=.1,ymin=-2.7,ymax=-2.7)  +
  geom_segment(aes(x =-0.2, xend = -1.1, y = -2, yend = -2),
               arrow = arrow(length = unit(0.2, "cm")), color = "black") +
  geom_segment(aes(x =0.2, xend = 1.1, y = -2, yend = -2),
               arrow = arrow(length = unit(0.2, "cm")), color = "black") +
  coord_cartesian(ylim = c(0,25),  clip = 'off')

dev.off()

# write.csv(aldex2_KO1_sig3, "aldex2_KO1_sig3_temp.csv")
aldex2_KO1_sig3 = read.csv("aldex2_KO1_sig3_temp.csv")
aldex2_KO1_sig3 <- aldex2_KO1_sig3[order(aldex2_KO1_sig3$effect, decreasing = T),]
aldex2_KO1_sig3$path = factor(aldex2_KO1_sig3$path, aldex2_KO1_sig3$path)
aldex2_KO1_sig3 = aldex2_KO1_sig3 %>% mutate(color = ifelse(effect >0, "#B40F20","#46ACC8"))
text_left <- textGrob("Increase in USA", gp=gpar(fontsize=16))
text_right <- textGrob("Increase in AUS", gp=gpar(fontsize=16))
text_center = textGrob("Center", gp=gpar(fontsize=16))

png("./KO1.png", width = 15, height = 8.5, units = "in", res = 300)
ggplot(aldex2_KO1_sig3, aes(y = path, x = effect)) +
  geom_point(aes(size = -log(we.eBH)), color = aldex2_KO1_sig3$color) + 
  xlab("") +
  ylab("") +
  scale_size_continuous("Adjusted p-value", 
                        breaks = c(3.9, 4.6, 5.3),
                        labels = c("0.02", "0.01", "0.005")) + 
  scale_color_manual("Effect", values = c("#B40F20", "#46ACC8")) +
  # breaks = c("#B40F20", "#46ACC8"),
  # labels = c("Increase in AUS",
  #            "Increase in USA")) +
  theme_classic() + 
  xlim(c(-1.2,1.2)) +
  theme(legend.position = "bottom",
        legend.title = element_text(size = 17),
        legend.text = element_text(size = 14.5),
        plot.margin = unit(c(1.5,1,2,1), "lines"),
        axis.title.x = element_text(vjust=0),
        axis.title = element_text(size = 17, face = "bold"),
        axis.text = element_text(size = 17, face = "bold")) +
  annotation_custom(text_left,xmin=-1.2,xmax=-.4,ymin=-2.7,ymax=-2.7) +
  annotation_custom(text_right,xmin=.4,xmax=1.2,ymin=-2.7,ymax=-2.7)  +
  # annotation_custom(text_center,xmin=-.1,xmax=.1,ymin=-2.7,ymax=-2.7)  +
  geom_segment(aes(x =-0.2, xend = -1.1, y = -2, yend = -2),
               arrow = arrow(length = unit(0.2, "cm")), color = "black") +
  geom_segment(aes(x =0.2, xend = 1.1, y = -2, yend = -2),
               arrow = arrow(length = unit(0.2, "cm")), color = "black") +
  coord_cartesian(ylim = c(0,25),  clip = 'off')
dev.off()


aldex2_PW1_sig2 <- aldex2_PW1_sig2[order(aldex2_PW1_sig2$effect, decreasing = T),]
aldex2_PW1_sig2$path = factor(aldex2_PW1_sig2$path, aldex2_PW1_sig2$path)
aldex2_PW1_sig2 = aldex2_PW1_sig2 %>% mutate(color = ifelse(effect >0, "#B40F20","#46ACC8"))
text_left <- textGrob("Increase in USA", gp=gpar(fontsize=18))
text_right <- textGrob("Increase in AUS", gp=gpar(fontsize=18))

# png("./PW1.png", width = 12.5, height = 7, units = "in", res = 300)
ggplot(aldex2_PW1_sig2,aes(y=path,x=effect))+
  geom_point(aes(size = -log(we.eBH)), color = aldex2_PW1_sig2$color) + 
  xlab("") +
  ylab("") +
  scale_size_continuous("Adjusted p-value",
                        breaks = c(3.2, 3.5, 3.9),
                        labels = c("0.04", "0.03", "0.02")) +
  scale_color_manual("Effect", values = c("#B40F20", "#46ACC8")) +
  # breaks = c("#B40F20", "#46ACC8"),
  # labels = c("Increase in AUS",
  #            "Increase in USA")) +
  xlim(c(-1.2,1.2)) +
  theme_classic() + 
  theme(legend.position = "bottom",
        legend.title = element_text(size = 17),
        legend.text = element_text(size = 14.5),
        plot.margin = unit(c(1.5,1,2,1), "lines"),
        axis.title.x = element_text(vjust=0),
        axis.title = element_text(size = 17, face = "bold"),
        axis.text = element_text(size = 17, face = "bold")) +
  annotation_custom(text_left,xmin=-1.2,xmax=-.4,ymin=-2.7,ymax=-2.7) +
  annotation_custom(text_right,xmin=.4,xmax=1.2,ymin=-2.7,ymax=-2.7)  +
  # annotation_custom(text_center,xmin=-.1,xmax=.1,ymin=-2.7,ymax=-2.7)  +
  geom_segment(aes(x =-0.2, xend = -1.1, y = -2, yend = -2),
               arrow = arrow(length = unit(0.2, "cm")), color = "black") +
  geom_segment(aes(x =0.2, xend = 1.1, y = -2, yend = -2),
               arrow = arrow(length = unit(0.2, "cm")), color = "black") +
  coord_cartesian(ylim = c(0,17),  clip = 'off')
# dev.off()

