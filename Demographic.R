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
library(table1)


dirname = dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(dirname)
IDs = c("IU131", "IU130", "IU123", "WU124", "WU119", "WU121",
         "IU128", "WU123", "IU121", "IU127", "IU125", "IU124",
         "IU120", "IU132", "IU129", "M135",  "M151",  "M124", 
         "M122",  "M121",  "M126",  "M136",  "P240",  "M143", 
         "M133",  "M153",  "P216",  "M155",  "M140",  "M161", 
         "M141",  "M149","M142",  "M152",  "M129",  "M137", 
         "M154",  "M119", "M146",  "M134",  "M120",  "M130", 
         "M160",  "M145",  "M148" )
# BAL.phy = readRDS("../../Data/phy_noncontam.rds")
# # only keep the samples
# BAL.phy = subset_samples(BAL.phy, is.neg == FALSE)
# meta = BAL.phy@sam_data
# meta$Visit[which(meta$Visit == "")] = "IASM"
# meta$Sequence_ID = rownames(meta)
# #IASM samples 
# idx = grepl("^[0-9]", meta$Sample_ID)
# meta$Sample_ID[idx] = meta$Subject_ID[idx]
# rownames(meta) = meta$Sample_ID
# colnames(BAL.phy@otu_table) = meta$Sample_ID
# # orderID = c(paste0("IASM", seq(1:20)))
# # orderID = orderID[-which(orderID=="IASM13")]
# # orderID = c(orderID, sort(meta$Sample_ID)[c(20:43, 119:128, 44:118)])
# orderID = sort(meta$Sample_ID)[c(20:43, 119:128, 44:118)]
# BAL.phy = phyloseq(BAL.phy@otu_table, BAL.phy@tax_table, meta)
# ## remove IASM samples 
# ## remove IU121 BA2 BA7 since IU121 received prophylax 
# BAL.phy = subset_samples(BAL.phy, Study == "ECF" & Subject_ID != "IU121")
# orderID = orderID[-which(orderID %in% c("IU121BA2", "IU121BA7"))]
# 
# sum(BAL.phy@sam_data$Visit == "BA2" & BAL.phy@sam_data$Country == "USA") #18
# sum(BAL.phy@sam_data$Visit == "BA2" & BAL.phy@sam_data$Country == "AUS") #32
# sum(BAL.phy@sam_data$Visit == "BA7" & BAL.phy@sam_data$Country == "USA") #14
# sum(BAL.phy@sam_data$Visit == "BA7" & BAL.phy@sam_data$Country == "AUS") #43
# sum(BAL.phy@sam_data$Visit == "BA7")

######################## Part 0: demo ######################## 
# USABA2 = BAL.phy@sam_data$Subject_ID[BAL.phy@sam_data$Visit == "BA2" & BAL.phy@sam_data$Country == "USA"]
# USABA7 = BAL.phy@sam_data$Subject_ID[BAL.phy@sam_data$Visit == "BA7" & BAL.phy@sam_data$Country == "USA"]
# AUSBA2 = BAL.phy@sam_data$Subject_ID[BAL.phy@sam_data$Visit == "BA2" & BAL.phy@sam_data$Country == "AUS"]
# AUSBA7 = BAL.phy@sam_data$Subject_ID[BAL.phy@sam_data$Visit == "BA7" & BAL.phy@sam_data$Country == "AUS"]
# BA2 = BAL.phy@sam_data$Sample_ID[BAL.phy@sam_data$Visit == "BA2"]
# BA7 = BAL.phy@sam_data$Sample_ID[BAL.phy@sam_data$Visit == "BA7"]
# IDs = BAL.phy@sam_data$Sample_ID 
# length(BAL.phy@sam_data$Sample_ID[BAL.phy@sam_data$Country == "USA"]) #32
# length(BAL.phy@sam_data$Sample_ID[BAL.phy@sam_data$Country == "AUS"]) #75


files<- list.files(paste0(str_remove(dirname, "BAL_species/Code"), 
                          "Clinical Data/Ben_clin/"), pattern = ".csv")
files1 <- paste0(paste0(str_remove(dirname, "BAL_species/Code"), 
                        "Clinical Data/Ben_clin/"), files)
data<-lapply(files1, function(i){
  read.csv(i, header = T)
})

names<-strsplit(files, split = ".csv")
names(data)<- names

culture<-full_join(data$demographics, data$culture, by = "Study_ID")
cols<- colnames(culture[,12:48])
culture[cols]<- lapply(culture[cols], factor)

cytokine<-full_join(data$demographics, data$cytokines, by = "Study_ID")
expiratory_ct<-full_join(data$demographics, data$expiratory_ct, by = "Study_ID")
inspiratory_ct<-full_join(data$demographics, data$inspiratroy_ct, by = "Study_ID")
pfts<-full_join(data$demographics, data$pfts, by = "Study_ID")

demo<- data$demographics
demo$Visit = ifelse(demo$Visit == "T1", "BA2", "BA7")
table1(~Site+Gender+Race+Hispanic+Genotype_F508+Age+Prophylaxis+Weight+Height|Country + Visit,demo[demo$Study_ID %in% IDs,])

table1(~FVCp+FEV0.5p+ FEF2575p+FRCp|Country+Visit.y,
       data= pfts[pfts$SampleID.y %in% IDs,])


table1(~No_Organisms_Detected+Oropharyngeal_Flora+Any_Bacterial_Growth+
         Any_Fungal_Growth+Any_Virus_Detected+Any_Bacterial_Fungal_Growth+
         Any_Classic_CF_Pathogen|Country, data= culture[culture$Study_ID %in% USABA2,])
table1(~No_Organisms_Detected+Oropharyngeal_Flora+Any_Bacterial_Growth+
         Any_Fungal_Growth+Any_Virus_Detected+Any_Bacterial_Fungal_Growth+
         Any_Classic_CF_Pathogen|Country, data= culture[culture$Study_ID %in% USABA7,])
table1(~No_Organisms_Detected+Oropharyngeal_Flora+Any_Bacterial_Growth+
         Any_Fungal_Growth+Any_Virus_Detected+Any_Bacterial_Fungal_Growth+
         Any_Classic_CF_Pathogen|Country, data= culture[culture$Study_ID %in% AUSBA2,])
table1(~No_Organisms_Detected+Oropharyngeal_Flora+Any_Bacterial_Growth+
         Any_Fungal_Growth+Any_Virus_Detected+Any_Bacterial_Fungal_Growth+
         Any_Classic_CF_Pathogen|Country, data= culture[culture$Study_ID %in% AUSBA7,])


table1(~NE_corrected_ug_ml+ IL8_corrected_pg_ml+ TNFa_corrected_pg_ml+
         IL6_corrected_pg_ml+IL1B_corrected_pg_ml+percent_neutrophils+
         total_cell_count+ abs_neutrophil|
         Country,
       data= cytokine[cytokine$Visit== "1B" & cytokine$Country == "USA",])


table1(~Healthy_perc+Bronchiectasis_perc+MucusPlugging_perc+
         Abnormal_perc+Atelectasis_perc+Disease_perc
       |
         Country,
       data= inspiratory_ct[inspiratory_ct$Country == "USA",])

table1(~Healthy_perc+Bronchiectasis_perc+MucusPlugging_perc+
         Abnormal_perc+Atelectasis_perc+Disease_perc
       |
         Country,
       data= inspiratory_ct[inspiratory_ct$Country == "Aus",])

table1(~Healthy_perc+ Airtrapping_perc
       |
         Country,
       data= expiratory_ct)
table1(~FVC_ZScore+FEV0.5_ZScore+FEF75_ZScore+ FEF2575_Z.Score+
         FRC_ZScore
       |
         Country,
       data= pfts[pfts$Visit.y=="Visit 1A" & pfts$Country == "USA",])
table1(~FVC_ZScore+FEV0.5_ZScore+FEF75_ZScore+ FEF2575_Z.Score+
         FRC_ZScore
       |
         Country,
       data= pfts[pfts$Visit.y=="Visit 1A" & pfts$Country == "Aus",])
table1(~FVC_ZScore+FEV0.5_ZScore+FEF75_ZScore+ FEF2575_Z.Score+
         FRC_ZScore
       |
         Country,
       data= pfts[pfts$Visit=="Visit 5A" & pfts$Country == "USA",])
table1(~FVC_ZScore+FEV0.5_ZScore+FEF75_ZScore+ FEF2575_Z.Score+
         FRC_ZScore
       |
         Country,
       data= pfts[pfts$Visit=="Visit 5A" & pfts$Country == "Aus",])


