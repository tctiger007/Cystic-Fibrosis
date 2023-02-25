library(dplyr)
library(ggplot2)
library(readxl)
library(stringr)
library(phyloseq)
library(DESeq2)
library(MASS)
library(phyloseq)
library(reshape2)
library(forcats)

dirname = dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(dirname)
# BAL.phy = readRDS("../../Data/phy_noncontam.rds")
# # only keep the samples
# BAL.phy = subset_samples(BAL.phy, Sample_type == "BAL")
# meta = BAL.phy@sam_data
# meta$Visit[which(meta$Visit == "")] = "IASM"
# meta$Sequence_ID = rownames(meta)
# ## Remove IASM samples
# BAL.phy = subset_samples(BAL.phy, Visit %in% c("BA2", "BA7"))
# meta = BAL.phy@sam_data
# 
# BAL.BA2.phy = subset_samples(BAL.phy, Visit == "BA2")
# BAL.BA7.phy = subset_samples(BAL.phy, Visit == "BA7")
# BA2.subID = BAL.BA2.phy@sam_data$Subject_ID
# BA7.subID = BAL.BA7.phy@sam_data$Subject_ID
# length(intersect(BA2.subID, BA7.subID)) #45

PFT = readRDS("../PFT.rds")
load("my_sample_col.RData")
sample.order1 = c("IU120BA7", "IU121BA7", "IU131BA2", "M134BA2",  "M135BA2",  "M136BA2", 
                  "M137BA2", "M140BA2",  "M140BA7",  "M141BA2",  "M141BA7",  "M151BA2", 
                  "M151BA7", "M153BA2",  "M155BA7" )
sample.order2 = c("IU123BA2", "IU123BA7", "IU124BA2", "IU125BA2", "IU128BA2", "IU129BA2",
                  "IU129BA7", "IU130BA7", "IU131BA7", "M134BA7",  "M135BA7",  "M136BA7" ,
                  "M137BA7",  "M155BA2")
sample.order3 = c("IU120BA2", "IU121BA2", "IU124BA7", "IU125BA7", "IU127BA2", "IU127BA7",
                  "IU128BA7", "IU130BA2", "M153BA7" )
### alpha diversity 
alpha_diversity = function(alpha = c("Shannon", "Chao1", "Observed","Simpson")){
  # alpha = "Shannon"
  alpha.BAL = read.csv(paste0("../../Data/", alpha, "_BAL_decontam.csv"))
  colnames(alpha.BAL)[1] = "Sample_ID" 
  alpha.c1 = alpha.BAL %>% filter(Sample_ID %in% sample.order1)
  alpha.c2 = alpha.BAL %>% filter(Sample_ID %in% sample.order2)
  alpha.c3 = alpha.BAL %>% filter(Sample_ID %in% sample.order3)
  alpha.clusters = rbind.data.frame(alpha.c1, alpha.c2, alpha.c3)
  alpha.clusters$Cluster = c(rep(1,dim(alpha.c1)[1]),
                             rep(2,dim(alpha.c2)[1]),
                             rep(3,dim(alpha.c3)[1]))
  alpha.clusters$Cluster = factor(alpha.clusters$Cluster)
  aov.alpha = aov(get(alpha, alpha.clusters)~ alpha.clusters$Cluster)
  pval = summary(aov.alpha)[[1]][1,5]
  aov_residuals <- residuals(object = aov.alpha)
  print(paste0("shapiro test p value: ", 
               shapiro.test(x = aov_residuals)[[2]])) # 0.8091
  print(TukeyHSD(aov.alpha))
  country = substring(alpha.clusters$Sample_ID,1,1)
  alpha.clusters$Country = ifelse(country == "M", "AUS", "USA")
  alpha.clusters$Country = factor(alpha.clusters$Country, 
                                  levels = c("USA", "AUS"))
  p = ggplot(alpha.clusters, aes(x = Cluster, y = eval(parse(text = alpha)), fill = Cluster) ) +
    geom_boxplot() +
    geom_jitter(aes(shape = Country), size = 3) +
    scale_fill_manual(values=c("orange2", "olivedrab4", "mediumorchid")) + 
    scale_shape_manual(values=c(16, 17))+
    ylab(alpha) + 
    theme_classic() + 
  theme(axis.title = element_text(size = 22, face = "bold"),
        axis.text = element_text(size = 20),
        plot.title = element_text(size = 20, face = "bold"),
        legend.title = element_text(size = 20, face = "bold"),
        legend.text = element_text(size = 20, face = "bold"))#,
        # legend.position = "none")
  return (p)
}

# png("Shannon_3cluster.png", width = 8, height = 8, units = "in", res=300)
alpha_diversity("Shannon") + 
  annotate("segment", x = 1, xend = 2, y = 3.3, yend = 3.3) +
  annotate("segment", x = 1, xend = 3, y = 3.6, yend = 3.6) +
  annotate("text",label="**", x = 1.5, y = 3.4, size = 10) +
  annotate("text",label="**", x = 2, y = 3.7, size = 10)
# dev.off()
# diff        lwr      upr     p adj
# 2-1 0.7588332  0.2082569 1.309410 0.0050717
# 3-1 0.9276497  0.3029563 1.552343 0.0024927
# 3-2 0.1688165 -0.4641880 0.801821 0.7920959

# png("Chao1_3cluster.png", width = 8, height = 8, units = "in", res=300)
alpha_diversity("Chao1")
# dev.off()

# png("Simpson_3cluster.png", width = 8, height = 8, units = "in", res=300)
alpha_diversity("Simpson")+ 
  annotate("segment", x = 1, xend = 2, y = 1.1, yend = 1.1) +
  annotate("segment", x = 1, xend = 3, y = 1.25, yend = 1.25) +
  annotate("text",label="**", x = 1.5, y = 1.15, size = 10) +
  annotate("text",label="**", x = 2, y = 1.3, size = 10) 
# dev.off()
# diff         lwr       upr     p adj
# 2-1 0.23716164  0.07632912 0.3979942 0.0026736
# 3-1 0.27268760  0.09020425 0.4551709 0.0023399
# 3-2 0.03552595 -0.14938522 0.2204371 0.8856788

# png("Observed_3cluster.png", width = 8, height = 8, units = "in", res=300)
alpha_diversity("Observed")+ 
  annotate("segment", x = 1, xend = 2, y = 60, yend = 60) +
  annotate("segment", x = 1, xend = 3, y = 66, yend = 66) +
  annotate("text",label="*", x = 1.5, y = 62, size = 10) +
  annotate("text",label=".", x = 2, y = 72, size = 14)
# dev.off()

Shannon.BAL = read.csv("../../Data/Shannon_BAL_decontam.csv")
Chao.BAL = read.csv("../../Data/Chao1_BAL_decontam.csv")
Observed.BAL = read.csv("../../Data/Observed_BAL_decontam.csv")
Simpson.BAL = read.csv("../../Data/Simpson_BAL_decontam.csv")
colnames(Shannon.BAL)[1] = "Sample_ID"
colnames(Chao.BAL)[1] = "Sample_ID"
colnames(Observed.BAL)[1] = "Sample_ID"
colnames(Simpson.BAL)[1] = "Sample_ID"

commonID.T1T2 = c("IU120", "IU121", "IU123", "IU124", "IU125", "IU127", "IU128", "IU129",
                  "IU130", "IU131", "M134",  "M135",  "M136",  "M137",  "M140",  "M141", 
                  "M151",  "M153",  "M155")
PFT$Sample_ID = paste0(PFT$SampleID, PFT$Visit)
PFT = PFT %>% dplyr::inner_join(Shannon.BAL, by = "Sample_ID") %>% 
  dplyr::inner_join(Chao.BAL, by = "Sample_ID") %>% 
  dplyr::inner_join(Observed.BAL, by = "Sample_ID") %>% 
  dplyr::inner_join(Simpson.BAL, by = "Sample_ID")
PFT = PFT[which(PFT$SampleID %in% commonID.T1T2),]


ggplot(PFT, aes(x = Simpson, y = Shannon)) + 
  geom_point() +
  stat_smooth(method = "lm")

ggplot(PFT, aes(x = Chao1, y = Shannon)) + 
  geom_point() +
  stat_smooth(method = "lm")

ggplot(PFT, aes(x = Observed, y = Shannon)) + 
  geom_point() +
  stat_smooth(method = "lm")

ggplot(PFT, aes(x = Shannon, y = Pfev0.5))+ 
  geom_point() +
  stat_smooth(method = "lm")
summary(lm(data = PFT, Pfev0.5~Shannon))
summary(lm(data = PFT, Pfev0.5~Chao1))
summary(lm(data = PFT, Pfev0.5~Simpson))
summary(lm(data = PFT, Pfev0.5~Observed))
# no sig association between fev0.5 and alpha diversity 

PFT$Cluster = ifelse(PFT$Sample_ID %in% sample.order1, 1,
                     ifelse(PFT$Sample_ID %in% sample.order2, 2, 3))
PFT$Cluster = as.factor(PFT$Cluster)
PFT$Group = paste0(PFT$Country, PFT$Visit)
PFT$Group = factor(PFT$Group, levels = c("USABA2", "USABA7",
                                         "AUSBA2", "AUSBA7"))
logit0 = glm(Cluster~Zfev0.5,
             data = PFT, family = "binomial")
logit1 = glm(Cluster~Shannon,
            data = PFT, family = "binomial")
summary(logit0)
summary(logit1)

logit2 = glm(Cluster~Shannon + Group ,
            data = PFT, family = "binomial")
summary(logit2)


