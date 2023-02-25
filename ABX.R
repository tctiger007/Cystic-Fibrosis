library(ggplot2)
library(readxl)
library(dplyr)
library(lme4)
library(sjstats)
library(lmerTest)
library(glmnet)

commonID.T1T2 = c("IU120","IU121", "IU123", "IU124", "IU125", 
                  "IU127", "IU128", "IU129", "IU130", "IU131",
                  "M134", "M135","M136", "M137", "M140", 
                  "M141", "M151", "M153", "M155")
########### ABX history 
## 1. plot abx frequency 
ABX.freq = read_excel("../Data/ABX_freq.xlsx")
ABX.freq$SubType = factor(ABX.freq$SubType,
                          levels = c("PE", "PE-ad", "AZITH",  
                                     "PX-Staphylococcus", "PX-Pseudomonas",
                                     "ERAD", "Other"))


## find when BA2 and BA7 occurred 
endpoint = read_excel("../Data/Clinical Data/VPECF Endpoint Data for James_2019_cleaned_WW.xlsx",
                      sheet = "PFTs")
endpoint1 = endpoint[which(endpoint$`Participant ID` %in% commonID.T1T2),]
endpoint1 = endpoint1[,c(1,2,6)]
colnames(endpoint1) = c("SampleID","Visit","Age")
endpoint1$Visit = lapply(str_split(endpoint1$Visit," "), "[[", 2)
endpoint1$Visit = ifelse(endpoint1$Visit == "1A", "BA2", "BA7")
endpoint1$Visit = factor(endpoint1$Visit, levels = c("BA2", "BA7"))
endpoint1$`Age(mo)` = endpoint1$Age/30

# png("medication_freq.png", height = 8, width = 8, unit = "in", res = 300)
ggplot() + 
  geom_point(data = ABX.freq, aes(x = Age, y = Subject, color = SubType,
                                  shape = SubType), size = 3) + 
  geom_point(data = endpoint1, aes(x = Age, y=SampleID), size = 2, shape = 3) +
  xlab("Age (Days)") +
  ylab("") + 
  scale_y_discrete(limits=rev) +
  scale_color_manual(name = "Medication Type",
                     labels = c("PE", "PE-ad", "AZITH", "PX-Staphylococcus", "PX-Pseudomonas",
                                "ERAD", "Other"),#, "Visit1", "Visit2"),
                     values = c("#d11141","#d11141", "#00b159", "#00aedb","#00aedb",
                                "#f37735", "#ffc425"))+ #, "black", "black")) +
  scale_shape_manual(name = "Medication Type",
                     labels = c("PE", "PE-ad", "AZITH", "PX-Staphylococcus", "PX-Pseudomonas",
                                "ERAD", "Other"),#, "Visit1", "Visit2"),
                     values = c(19,1,19,19,1,19,19))+#,3,4)) +
  theme_bw() +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.title.x = element_text(face = "bold")) 
# dev.off()

## 2. investigate whether the ABX usage is significantly different between 
## USA and AUS at T1 and T2. 
ABX.freq1 = read_excel("../Data/VPECF ABX SUMMARY.xlsx", sheet = "ABX_TwoTimePoints")
ABX.freq1$Country = rep(c(rep("USA", 10), rep("AUS",9)),2)
ABX.freq1$`PE-sum` = ABX.freq1$PE + ABX.freq1$`PE-ad`
ABX.freq1.BA2 = ABX.freq1 %>% filter(Visit == "BA2")
ABX.freq1.BA7 = ABX.freq1 %>% filter(Visit == "BA7")
kruskal.test(ABX.freq1.BA2$PE~Country, data=ABX.freq1.BA2) # PE: Not sig. 
kruskal.test(ABX.freq1.BA2$`PE-ad`~Country, data=ABX.freq1.BA2) # PE-ad: NA
kruskal.test(ABX.freq1.BA2$`PE-sum`~Country, data=ABX.freq1.BA2) # PE-ad: Not sig
kruskal.test(ABX.freq1.BA2$`PX-Staphylococcus`~Country, data=ABX.freq1.BA2) # anti-staph sig. difference
kruskal.test(ABX.freq1.BA2$ERAD~Country, data=ABX.freq1.BA2) # ERAD: not sig
kruskal.test(ABX.freq1.BA2$Other~Country, data=ABX.freq1.BA2) # Other: not sig

kruskal.test(ABX.freq1.BA7$PE~Country, data=ABX.freq1.BA7) # PE: Not sig. 
kruskal.test(ABX.freq1.BA7$`PE-ad`~Country, data=ABX.freq1.BA7) # PE-ad: Not sig
kruskal.test(ABX.freq1.BA7$`PE-sum`~Country, data=ABX.freq1.BA7) # PE-ad: Not sig
kruskal.test(ABX.freq1.BA7$`PX-Staphylococcus`~Country, data=ABX.freq1.BA7) # anti-staph sig. difference
kruskal.test(ABX.freq1.BA7$ERAD~Country, data=ABX.freq1.BA7) # ERAD: not sig
kruskal.test(ABX.freq1.BA7$Other~Country, data=ABX.freq1.BA7) # Other: sig



## 2. GLM: association between PFT and ABX, microbiome, and etc. 
PFT = readRDS("../Data/PFT.rds")
PFT = PFT %>% mutate(SampleID = paste0(PFT$SampleID, PFT$Visit))
shannon.BAL = read.csv("../Data/shannon_BAL_decontam.csv")
colnames(shannon.BAL)[1] = "SampleID"
cluster = readRDS("cluster.rds")

clinic.df = ABX.freq1 %>% mutate(SampleID = paste0(ABX.freq1$Subject, ABX.freq1$Visit)) %>% 
  inner_join(PFT, by = "SampleID") %>% 
  inner_join(shannon.BAL, by = "SampleID")
clinic.df = clinic.df[,-grep(".y", colnames(clinic.df), fixed = T)]
colnames(clinic.df)[grep(".x", colnames(clinic.df), fixed = T)] = str_remove(colnames(clinic.df)[grep(".x", colnames(clinic.df), fixed = T)], ".x")
clinic.df$Cluster = cluster$Cluster[match(rownames(cluster), clinic.df$SampleID)]
clinic.df$Cluster = factor(clinic.df$Cluster, levels = c(1,2,3))
BAL.phy_genus_relabun = readRDS("BAL.phy_genus_relabun.rds")
OTU = as.data.frame(otu_table(BAL.phy_genus_relabun)@.Data)
rownames(OTU) = tax_table(BAL.phy_genus_relabun)[,"Genus"]
taxa.order = names(sort(rowSums(OTU), decreasing = T))
OTU.top = OTU[match(taxa.order, rownames(OTU)),]
OTU.top = as.data.frame(t(OTU.top))
OTU.top$SampleID = rownames(OTU.top)
clinic.df = clinic.df %>% inner_join(OTU.top, by = "SampleID")
demo = read.csv("../../../Clinical Data/Ben_clin/demographics.csv")
colnames(demo)[12] = "SampleID"
clinic.df = clinic.df %>% inner_join(demo, by = "SampleID")
clinic.df = clinic.df[,-grep(".y", colnames(clinic.df), fixed = T)]
colnames(clinic.df)[grep(".x", colnames(clinic.df), fixed = T)] = str_remove(colnames(clinic.df)[grep(".x", colnames(clinic.df), fixed = T)], ".x")

colnames(clinic.df)
# write.csv(clinic.df, "clinic.csv")

lm = lmer(formula = Pfev0.5 ~ 1 + `PX-Staphylococcus` + Shannon  + (1|Subject),
          data = clinic.df)
summary(lm)
icc(lm)

lm = lmer(formula = Pfev0.5 ~ 1 + Streptococcus  + (1|Subject),
          data = clinic.df)
summary(lm)
icc(lm)

ggplot(clinic.df, aes(x=Streptococcus, y = `Pfev0.5`)) + 
  geom_point() + 
  geom_smooth(method = "lm")

ggplot(clinic.df, aes(x=Gemella, y = `Pfev0.5`)) + 
  geom_point() + 
  geom_smooth(method = "lm")

ggplot(clinic.df, aes(x=Neisseria, y = `Pfev0.5`)) + 
  geom_point() + 
  geom_smooth(method = "lm")

ggplot(clinic.df, aes(x=Staphylococcus, y = `Pfev0.5`)) + 
  geom_point() + 
  geom_smooth(method = "lm")

ggplot(clinic.df, aes(x=Fusobacterium, y = `Pfev0.5`)) + 
  geom_point() + 
  geom_smooth(method = "lm")

ggplot(clinic.df, aes(x=Bradyrhizobium, y = `Pfev0.5`)) + 
  geom_point() + 
  geom_smooth(method = "lm")

lm0 = lm(Pfev0.5~ Streptococcus +  Gemella   +  Neisseria + 
           Bradyrhizobium + `Clostridium sensu stricto 10` + Rothia +   
           Alloprevotella + Cloacibacterium  + Granulicatella  +            
           Moraxella + Staphylococcus + Fusobacterium, data = clinic.df)
summary(lm0)

ggplot(data = clinic.df, 
       aes(x   = `PX-Staphylococcus`,
           y   = Pfev0.5,
           color = Cluster))+
  geom_point(size     = 1, 
             alpha    = .7, 
             position = "jitter")+
  geom_smooth(method   = lm,
              se       = T, 
              size     = 1.5, 
              linetype = 1, 
              alpha    = .7)+
  theme_minimal() +
  scale_color_manual(name = "Cluster", 
                     labels = c("1","2","3"),
                     values = c("orange2","olivedrab4","mediumorchid"))




# i. FEV0.5 is not linearly associated with Shannon diversity
## not linearly associated in three clusters either
lm1 = lm(Pfev0.5~ Shannon, data = clinic.df)
summary(lm1)
temp = clinic.df %>% filter(Cluster == 2)
summary(lm(temp$Pfev0.5~temp$Shannon))
ggplot(clinic.df, aes(x = Shannon, y = Pfev0.5, Color = Cluster)) + 
  geom_point() +
  geom_smooth(method = lm)

# ii. 
lm2 = lm(Pfev0.5~ `PX-Staphylococcus`, data = clinic.df)
summary(lm2)
ggplot(clinic.df, aes( y = Pfev0.5, x = `PX-Staphylococcus`)) +
  geom_point() + 
  geom_smooth(method = lm)
  # geom_boxplot() +
  # facet_grid(~Cluster)

ggplot(clinic.df, aes(y = Pfev0.5, x = SampleID, color = Country)) + geom_point()


colnames(ABX.freq)[1] = "SampleID" 
glm.df = ABX.freq %>% inner_join(PFT, by = "SampleID")


ABX = read_excel("./Data/VPECF ABX SUMMARY.xlsx")
ABX = ABX[-which(is.na(ABX$Subject)),]
ABX$`# Exacerbations oral antibiotics` = lapply(str_split(ABX$`# Exacerbations oral antibiotics (dates)`,
                                                          " "), "[[", 1)
ABX$`# Exacerbations intravenous antibiotics` = lapply(str_split(ABX$`# Exacerbations intravenous antibiotics (dates)`,
                                                                 " "), "[[", 1)

ABX$`# Other antibiotic use` = lapply(str_split(ABX$`# Other antibiotic use (age)`,
                                                " "), "[[", 1)

ABX$`# Exacerbations oral antibiotics`[is.na(ABX$`# Exacerbations oral antibiotics`)] = 0
ABX$`# Exacerbations intravenous antibiotics`[is.na(ABX$`# Exacerbations intravenous antibiotics`)] = 0
ABX$`# Other antibiotic use`[is.na(ABX$`# Other antibiotic use`)] = 0
ABX = ABX %>% mutate(SampleID = str_remove(ABX$Subject, "1C"))

PFT.new = ABX %>% inner_join(PFT, by = "SampleID")





