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
library(latex2exp)
library(wesanderson)


dirname = dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(dirname)
phy = readRDS("../../Data/phy_noncontam.rds")
meta = meta(phy)
# meta$Visit = ifelse(meta$Visit == "BA2", "T1", 
#                     ifelse(meta$Visit == "BA7", "T2", ""))
BAL.BA2 = meta %>% dplyr::filter(Sample_type == "BAL" & Visit == "BA2") #51
BAL.BA7 = meta %>% dplyr::filter(Sample_type == "BAL" & Visit == "BA7") #58
BA2_BA7_ID = intersect(BAL.BA2$Subject_ID, BAL.BA7$Subject_ID) 
# 45; we keep these 45 matched BAL samples. 
# Now we find the matched prewash samples 
prewash = meta %>% dplyr::filter(Sample_type == "Prewash")
prewash.BA2 = prewash %>% dplyr::filter(Sample_ID %in% paste0(BA2_BA7_ID, "BA2"))
prewash.BA7 = prewash %>% dplyr::filter(Sample_ID %in% paste0(BA2_BA7_ID, "BA7"))
prewash_ID = c(prewash.BA2$Sample_ID, prewash.BA7$Sample_ID) #51

# BAL phyloseq object that contains matched BA2 and BA7 samples 
BAL.phy = subset_samples(phy, Subject_ID %in% BA2_BA7_ID)
BAL.phy = subset_samples(BAL.phy, Sample_type == "BAL")
BAL.meta = meta(BAL.phy)
prewash.phy = subset_samples(phy, Sample_type == "Prewash")
prewash.phy = subset_samples(prewash.phy, Sample_ID %in% prewash_ID)
orderID = sort(BAL.meta$Sample_ID)[c(1:22,83:90,23:82)]
prewash.orderID = sort(prewash.phy@sam_data$Sample_ID)[c(1:14,50:51,15:49)]

##############################Part:  qPCR ##############################
qPCR = read_excel("../../Data/qPCR Results March 2020.xlsx",
                  sheet = "qPCR_Results")
qPCR$`Average copies/g sample Adjust for dilution`[which(qPCR$`Average copies/g sample` == "Undetermined")] = 0

# positive controls are stool samples
qPCR$`Specimen Type` = factor(qPCR$`Specimen Type`,
                              levels = c("BAL", "Prewash", "Negative Control",
                                         "Stool Positive Control"))
wes_palette("FantasticFox1")[c(3,5)] 
# "#46ACC8" "#B40F20"

png("Figures/qPCR.png", width = 6, height = 6, units = "in", res = 300)
ggplot(qPCR %>% filter(`Specimen Type` %in% c("BAL", "Prewash")), 
       aes(x = `Specimen Type`,
                 y = `Average copies/g sample Adjust for dilution`,
                 fill = `Specimen Type`)) +
  geom_boxplot() +
  ylab("Average copies/g sample") +
  xlab("") +
  scale_fill_manual(values = wes_palette("FantasticFox1")[c(4,3)]) +
  scale_x_discrete(breaks = c("BAL", "Prewash"), #"Negative Control", "Stool Positive Control"),
                  labels = c("BAL", "Prewash")) +  #"Neg. Control", "Pos. Control")) +
  scale_y_continuous(trans = log10_trans(),
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x))) +
  theme_classic() + 
  theme(axis.text.x = element_text(size = 16, angle = 45, vjust = 1, hjust = 1, face = "bold"),
        axis.text.y = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 16, face = "bold")) +
  annotate("text", x = 2, y = 10^8, label = TeX("p = 7.45 $\\times\\10^{-11}$"))
dev.off()
################################################

qPCR.BAL = qPCR %>% dplyr::filter(`Specimen Type` == "BAL")
qPCR.prewash = qPCR %>% dplyr::filter(`Specimen Type` == "Prewash")
qPCR.df = qPCR.BAL %>% left_join(qPCR.prewash, by = "Sample Name") 

wilcox.test(qPCR.df$`Average copies/g sample Adjust for dilution.x`,
            qPCR.df$`Average copies/g sample Adjust for dilution.y`, 
            paired = TRUE,
            alternative = "two.sided")


######################## qPCR for BAL and prewash samples ######################## 
BA2_BA7_ID = c(paste0(BA2_BA7_ID, "BA2"), paste0(BA2_BA7_ID, "BA7"))
colnames(qPCR)[c(11, 14:15)] = c("Average copies/g sample without adjust",
                                 "Average copies/g sample", "SD")
qPCR.BAL = qPCR %>% dplyr::filter(`Specimen Type` == "BAL") %>% 
  dplyr::filter(`Sample Name` %in% BA2_BA7_ID)
qPCR.prewash = qPCR %>% dplyr::filter(`Specimen Type` == "Prewash") 
qPCR.prewash = qPCR.prewash %>% dplyr::filter(`Sample Name` %in% prewash_ID)
prewash.orderID = sort(qPCR.prewash$`Sample Name`)[c(1:14, 50:51, 15:49)]
qPCR.BAL = qPCR.BAL %>% 
  mutate(name = fct_relevel(`Sample Name`, orderID))
qPCR.prewash = qPCR.prewash %>% 
  mutate(name = fct_relevel(`Sample Name`, prewash.orderID))
qPCR = rbind.data.frame(qPCR.BAL, qPCR.prewash)
qPCR$name1 = ifelse(qPCR$`Specimen Type` == "BAL", qPCR$`Sample Name`, 
                    paste0(qPCR$`Sample Name`, "c"))

png("Figures/qPCR_individual.png", width = 12, height = 8, units = "in", res = 300)
ggplot(qPCR, 
       aes(x =name1, y = `Average copies/g sample`, fill = `Specimen Type`)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin=`Average copies/g sample`-SD, 
                    ymax=`Average copies/g sample`+SD), 
                width=.2, position=position_dodge(.9)) +
  geom_col(position = position_stack(reverse = TRUE)) +
  labs(x = "",
       y = expression("Average copies/g sample"),
       title = "") +
  scale_fill_manual(values = wes_palette("FantasticFox1")[c(4,3)]) +
  # facet_grid(.~`Specimen Type`) + 
  theme_classic() +
  theme(axis.text.x = element_text(size = 6, angle = 90, vjust = 0.5, hjust = 0),
        axis.text.y = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.key.height = unit(.5, 'cm'),
        legend.key.width = unit(.5, 'cm'),
        strip.text = element_text(size = 14),
        axis.title.y = element_text(size = 16),
        plot.title = element_text(size = 16)) +
  scale_y_continuous(trans = "log10", limits = c(1, 2*10^10),
                     breaks = c(10^2,10^4,10^6,10^8,10^10))
dev.off()



############################## Part: qPCR * relative abundance ##############################
BAL.phy_genus <- tax_glom(phy, "Genus", NArm = TRUE)
BAL.phy_genus_relabun <- transform_sample_counts(BAL.phy_genus, function(OTU) OTU/sum(OTU) * 100)
BAL.phy_genus_relabun_ECF = subset_samples(phy, Sample_ID %in% BA2_BA7_ID)
meta = BAL.phy_genus_relabun_ECF@sam_data
meta$sampleNames = paste0(meta$Sample_ID,
                          ifelse(meta$Sample_type == "BAL", "", "c"))
sample_data(BAL.phy_genus_relabun_ECF) = meta

qPCR_ECF = qPCR[match(meta$sampleNames, qPCR$name1),]
otu_table_abs = matrix()
qPCR_copies = qPCR_ECF$`Average copies/g sample`
# qPCR * relative abundance
for (i in 1:dim(qPCR_ECF)[1]){
  temp = BAL.phy_genus_relabun_ECF@otu_table[,i] * qPCR_copies[i]
  otu_table_abs = cbind.data.frame(otu_table_abs, temp)
}
otu_table_abs = otu_table_abs[,-1]
otu_table_abs = otu_table(otu_table_abs, taxa_are_rows = T)
BAL.phy_genus_abs = phyloseq(otu_table_abs,
                             BAL.phy_genus_relabun_ECF@tax_table,
                             BAL.phy_genus_relabun_ECF@sam_data)

orderTaxa1 = names(sort(taxa_sums(BAL.phy_genus_abs), decreasing=TRUE))
top20_genera1 <- orderTaxa1[1:20]

BAL.phy_genus_abs.df = psmelt(BAL.phy_genus_abs)  %>%
  mutate(name = fct_relevel(sampleNames, sort(meta$sampleNames)))
BAL.phy_genus_abs.df$Genus1 = BAL.phy_genus_abs.df$Genus
  
BAL.phy_genus_abs.df$Genus1[which(! BAL.phy_genus_abs.df$OTU %in% top20_genera1)] = "Low Abundance"

top20_genera_name1 = BAL.phy_genus_abs.df$Genus[match(top20_genera1, BAL.phy_genus_abs.df$OTU)]
BAL.phy_genus_abs.df = BAL.phy_genus_abs.df %>%
  mutate(Genusorder = fct_relevel(Genus1, c(top20_genera_name1, "Low Abundance")))
BAL.phy_genus_abs.df$Country = factor(BAL.phy_genus_abs.df$Country, levels = c("USA", "AUS"))
BAL.phy_genus_abs.df$Group = paste0(BAL.phy_genus_abs.df$Country, " ", BAL.phy_genus_abs.df$Visit)
BAL.phy_genus_abs.df$Group = factor(BAL.phy_genus_abs.df$Group, levels = c("USA BA2","AUS BA2",
                                                       "USA BA7", "AUS BA7"))
# BAL.phy_genus_abs.df$Abundance = BAL.phy_genus_abs.df$Abundance/(10^6)
## divide the abundance with 10^6

# png("./Results/genus_abundance_timesqPCR.png", width = 20, height = 10, units = "in", res = 300)
ggplot(BAL.phy_genus_abs.df, aes(x = sampleNames, y = Abundance, fill = Genusorder)) +
  geom_bar(stat = "identity") +
  geom_col(position = position_stack(reverse = TRUE)) +
  labs(x = "",
       y = expression("copies/g sample"),
       title = "") +
  scale_fill_manual(values = c(cols, "gray57")) +
  facet_grid(Country ~ Visit, scales = "free") +
  labs(fill = "Genus") +
  theme_classic() +
  theme(axis.text.x = element_text(size = 7, angle = 90, vjust = 0.5, hjust = 1),
        axis.text.y = element_text(size = 12),
        legend.text = element_text(size = 14),
        strip.text = element_text(size = 14),
        # strip.background = element_blank(),
        axis.title.y = element_text(size = 16),
        plot.title = element_text(size = 16)) +
  scale_y_continuous(trans = "log10", limits = c(1, 2*10^10),
                     breaks = c(10^2,10^4,10^6,10^8,10^10))
# dev.off()




# # BAL samples 
# ggplot(qPCR.BAL, 
#        aes(x =name, y = `Average copies/g sample`)) +
#   geom_bar(stat = "identity") +
#   geom_errorbar(aes(ymin=`Average copies/g sample`-SD, 
#                     ymax=`Average copies/g sample`+SD), 
#                 width=.2, position=position_dodge(.9)) +
#   geom_col(position = position_stack(reverse = TRUE)) +
#   labs(x = "",
#        y = expression("copies/g sample"),
#        title = "") +
#   theme_classic() +
#   theme(axis.text.x = element_text(size = 7, angle = 90, vjust = 0.5, hjust = 1),
#         axis.text.y = element_text(size = 12),
#         legend.text = element_text(size = 14),
#         strip.text = element_text(size = 14),
#         axis.title.y = element_text(size = 16),
#         plot.title = element_text(size = 16)) +
#   scale_y_continuous(trans = "log10", limits = c(1, 2*10^10),
#                      breaks = c(10^2,10^4,10^6,10^8,10^10))
# 
# # prewash
# ggplot(qPCR.prewash, 
#        aes(x =name, y = `Average copies/g sample`)) +
#   geom_bar(stat = "identity") +
#   geom_errorbar(aes(ymin=`Average copies/g sample`-SD, 
#                     ymax=`Average copies/g sample`+SD), 
#                 width=.2, position=position_dodge(.9)) +
#   geom_col(position = position_stack(reverse = TRUE)) +
#   labs(x = "",
#        y = expression("copies/g sample"),
#        title = "") +
#   theme_classic() +
#   theme(axis.text.x = element_text(size = 7, angle = 90, vjust = 0.5, hjust = 1),
#         axis.text.y = element_text(size = 12),
#         legend.text = element_text(size = 14),
#         strip.text = element_text(size = 14),
#         axis.title.y = element_text(size = 16),
#         plot.title = element_text(size = 16)) +
#   scale_y_continuous(trans = "log10", limits = c(1, 2*10^10),
#                      breaks = c(10^2,10^4,10^6,10^8,10^10))
# 

