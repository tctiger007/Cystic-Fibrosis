library(phyloseq)
library(ggplot2)
library(rstudioapi)
library(lubridate)
library(microbiome)
library(dplyr)
library(stringr)

setwd(dirname(getActiveDocumentContext()$path))
stool.phy = readRDS("../../../../stool/stool_phy.rds")
BAL.phy = readRDS("../../BAL_phy.rds")
IDs = c("IU120", "IU121", "IU123", "IU124", "IU125", "IU127", "IU128", "IU129", "IU130",
        "IU131", "IU132", "WU119", "WU121", "WU123", "WU124", "M119",  "M120",  "M121",  
        "M122",  "M124",  "M126",  "M129", "M130",  "M133",  "M134",  "M135",  "M136",  
        "M137",  "M140",  "M141",  "M142", "M143",  "M145",  "M146",  "M148",  "M149",  
        "M151",  "M152",  "M153",  "M154", "M155",  "M160",  "M161",  "P216",  "P240")  
BAL.phy = subset_samples(BAL.phy, Subject_ID %in% IDs)

stool.phy = subset_samples(stool.phy, Subject_ID %in% IDs)
stool.meta = sample_data(stool.phy)

BAL.meta = sample_data(BAL.phy)
stool.meta$DOB = BAL.meta$DOB[match(stool.meta$Subject_ID, BAL.meta$Subject_ID)]
stool.meta$Stool_Age = difftime(mdy(stool.meta$Date), mdy(stool.meta$DOB))

BAL.temp.meta = meta(BAL.meta)
BAL.temp.meta = BAL.temp.meta %>% 
  dplyr::select("Sample_ID", "Subject_ID", "Visit", "Age.at.BAL..day.")

BAL.temp.meta.wide = reshape(BAL.temp.meta, idvar = "Subject_ID", timevar = "Visit", direction = "wide")
rownames(BAL.temp.meta.wide) = BAL.temp.meta.wide$Subject_ID
# write.csv(BAL.temp.meta.wide, "BAL_meta_wide.csv")


stool.phy = subset_samples(stool.phy, Sample_ID != "IU131S5WK")
# remove sample MYCF00111

stool.temp.meta = meta(stool.meta)
stool.temp.meta = stool.temp.meta %>% 
  dplyr::select("Sample_ID", "Subject_ID", "Visit", "Stool_Age")
# write.csv(stool.temp.meta, "stool_meta.csv")
stool.temp.meta = read.csv("stool_meta.csv")
stool.temp.meta.wide = reshape(stool.temp.meta, idvar = "Subject_ID", timevar = "Visit1", direction = "wide")

merge.temp = stool.temp.meta.wide %>% left_join(BAL.temp.meta.wide, by = "Subject_ID")

D1 = merge.temp$Stool_Age.1 - merge.temp$Age.at.BAL..day..BA2
D2 = merge.temp$Stool_Age.2 - merge.temp$Age.at.BAL..day..BA2
D3 = merge.temp$Stool_Age.3 - merge.temp$Age.at.BAL..day..BA2
D4 = merge.temp$Stool_Age.4 - merge.temp$Age.at.BAL..day..BA2
D5 = merge.temp$Stool_Age.5 - merge.temp$Age.at.BAL..day..BA2
D6 = merge.temp$Stool_Age.6 - merge.temp$Age.at.BAL..day..BA2
diff.BA2 = cbind.data.frame(D1, D2, D3, D4, D5, D6)
diff.abs.BA2 = abs(diff.BA2)
which.min = apply(diff.abs.BA2, 1, FUN = which.min)
# [1] 219 141 113 147   0  87  17 106   7   3 300 256 233 226  68   1   1 185   1   1
# [21]   1   1   8  93   2 292   0  29  64
# which.min  [1] 1 1 1 1 1 1 2 1 1 3 1 1 1 1 1 1 2 1 2 1 1 1 1 1 1 1 2 2 2


D1.BA7 = merge.temp$Stool_Age.1 - merge.temp$Age.at.BAL..day..BA7
D2.BA7 = merge.temp$Stool_Age.2 - merge.temp$Age.at.BAL..day..BA7
D3.BA7 = merge.temp$Stool_Age.3 - merge.temp$Age.at.BAL..day..BA7
D4.BA7 = merge.temp$Stool_Age.4 - merge.temp$Age.at.BAL..day..BA7
D5.BA7 = merge.temp$Stool_Age.5 - merge.temp$Age.at.BAL..day..BA7
D6.BA7 = merge.temp$Stool_Age.6 - merge.temp$Age.at.BAL..day..BA7
diff.BA7 = cbind.data.frame(D1.BA7, D2.BA7, D3.BA7, D4.BA7, D5.BA7, D6.BA7)
diff.abs.BA7 = abs(diff.BA7)
which.min.BA7 = apply(diff.abs.BA7, 1, FUN = which.min)
# [1] 2 3 1 1 3 1 3 3 6 5 1 1 2 2 4 3 4 2 5 2 1 4 3 1 2 1 3 3 3

# try to find the samples of stools that were collected around the time of BA2 and BA7 collection
# if the sample was both close to BA2 and BA7, we will re-use the stool sample 
# length(unique(stool.meta$Subject_ID)) #29
# 29 subjects will be matched 

# write.csv(merge.temp, "merge_temp.csv")
# 1 1 1 1 1 1 2 1 1 3 1 1 1 1 1 1 2 1 2 1 1 1 1 1 1 1 2 2 2
# 2 3 1 1 3 1 3 3 6 5 1 1 2 2 4 3 4 2 5 2 1 4 3 1 2 1 3 3 3

which.min == which.min.BA7
# [1] FALSE FALSE  TRUE  TRUE FALSE  TRUE FALSE FALSE FALSE FALSE  TRUE  TRUE FALSE
# [14] FALSE FALSE FALSE FALSE FALSE FALSE FALSE  TRUE FALSE FALSE  TRUE FALSE  TRUE
# [27] FALSE FALSE FALSE


stool.phy.rare = rarefy_even_depth(stool.phy, rngseed=1, 
                                   sample.size=0.9*min(sample_sums(stool.phy)), replace=F)
richness = estimate_richness(stool.phy.rare)
shannon.stool = richness %>% dplyr::select("Shannon")
stool.meta$Sequence_ID = rownames(stool.meta)
shannon.stool$Sequence_ID = rownames(shannon.stool)
stool.meta$Sequence_ID = as.character(stool.meta$Sequence_ID )
shannon.stool = shannon.stool %>% left_join(stool.meta, by = "Sequence_ID")

shannon.BAL = read.csv("../../../Data/Shannon_BAL_decontam.csv")
numbers_only <- function(x) !grepl("\\D", x)
shannon.BAL = shannon.BAL[-which(numbers_only(shannon.BAL$X)),]

shannon.BAL$Visit = paste0("BA", sapply(str_split(shannon.BAL$X, "BA"), "[[",2))
shannon.BAL$Subject_ID = sapply(str_split(shannon.BAL$X, "BA"), "[[",1)
shannon.BA2 = shannon.BAL %>% dplyr::filter(Visit == "BA2")
shannon.BA7 = shannon.BAL %>% dplyr::filter(Visit == "BA7")

matched.BA2 = readxl::read_excel("matched_stool_samples.xlsx", sheet = 1)
matched.BA7 = readxl::read_excel("matched_stool_samples.xlsx", sheet = 2)
stool.BA2 = shannon.stool %>% dplyr::filter(Sequence_ID %in% matched.BA2$Sequence_ID)
stool.BA7 = shannon.stool %>% dplyr::filter(Sequence_ID %in% matched.BA7$Sequence_ID)


stool.BAL.BA2 = stool.BA2 %>% left_join(shannon.BA2, by = "Subject_ID")
stool.BAL.BA7 = stool.BA7 %>% left_join(shannon.BA7, by = "Subject_ID")

# Shannon.x: stool, Shannon.y: BAL


# png("stool_BAL_BA2.png", width = 12, height = 10, units = "in", res = 300)
ggplot(stool.BAL.BA2, aes(x=Shannon.x, y=Shannon.y)) +
  geom_point(size = 3) + 
  ylab("BAL") +
  xlab("Stool") +
  stat_smooth(method = lm, linewidth = 3) +
  ggtitle("BA2") +
  annotate(geom="text", x = 3, y = 3.5, label = "p = 0.00271", size = 10) +
  theme_classic()  + 
  theme(legend.title = element_text(size = 20),
        axis.title = element_text(size = 24, face = "bold"),
        axis.text = element_text(size = 24),
        title = element_text(size = 24, face = "bold")) 
# dev.off()
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)   
# (Intercept)               0.7614     0.3067   2.483  0.01955 * 
#   stool.BAL.BA2$Shannon.x   0.4967     0.1505   3.301  0.00271 **
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.6961 on 27 degrees of freedom
# Multiple R-squared:  0.2876,	Adjusted R-squared:  0.2612 
# F-statistic:  10.9 on 1 and 27 DF,  p-value: 0.002712

# png("stool_BAL_BA7.png", width = 12, height = 10, units = "in", res = 300)
ggplot(stool.BAL.BA7, aes(x=Shannon.x, y=Shannon.y)) +
  geom_point(size = 3) + 
  ylab("BAL") +
  xlab("Stool") +
  stat_smooth(method = lm, linewidth = 3) +
  ggtitle("BA7") +
  annotate(geom="text", x = 3.2, y = 3.2, label = "p = 0.467", size = 10) +
  theme_classic() + 
  theme(legend.title = element_text(size = 20),
        axis.title = element_text(size = 24, face = "bold"),
        axis.text = element_text(size = 24),
        title = element_text(size = 24, face = "bold")) 
# dev.off()
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)               1.8005     0.3806   4.731 8.23e-05 ***
#   stool.BAL.BA7$Shannon.x   0.1231     0.1666   0.739    0.467    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.6775 on 24 degrees of freedom
# Multiple R-squared:  0.02224,	Adjusted R-squared:  -0.0185 
# F-statistic: 0.546 on 1 and 24 DF,  p-value: 0.4671

colnames(stool.BAL.BA2)[which(colnames(stool.BAL.BA2) == "Sample_ID")] = "Sample_ID_stool"
colnames(stool.BAL.BA2)[which(colnames(stool.BAL.BA2) == "X")] = "Sample_ID"
colnames(stool.BAL.BA7)[which(colnames(stool.BAL.BA7) == "Sample_ID")] = "Sample_ID_stool"
colnames(stool.BAL.BA7)[which(colnames(stool.BAL.BA7) == "X")] = "Sample_ID"

### check the association between gut microbiome (shannon) and FEV
PFT = read.csv("PFT.csv")
merge.BA2 = stool.BAL.BA2 %>% left_join(PFT, by="Sample_ID")
merge.BA7 = stool.BAL.BA7 %>% left_join(PFT, by="Sample_ID")
merge = rbind.data.frame(merge.BA2, merge.BA7)


merge1 = merge[!is.na(merge$Pfev0.5),]
unique(merge1$Sequence_ID) #24
merge1$Sequence_ID[which(duplicated(merge1$Sequence_ID))]
merge1 = merge1[-which(duplicated(merge1$Sequence_ID)),]
summary(lm(Pfev0.5~Shannon.x, merge1)) # 0.11
summary(lm(Pfvc~Shannon.x, merge1)) # 0.2563

# saveRDS(merge1, "merge1.rds")
# saveRDS(stool.phy, "stool_phy.rds")
