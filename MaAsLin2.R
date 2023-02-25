library(Maaslin2)
library(microbiome)
library(phyloseq)
library(maditr)
library(gautils2)
library(metagMisc)
library(wesanderson)
BAL.phy = readRDS("../BAL_phy.rds")
IDs = c("IU120", "IU121", "IU123", "IU124", "IU125", "IU127", "IU128", "IU129", "IU130",
        "IU131", "IU132", "WU119", "WU121", "WU123", "WU124", "M119",  "M120",  "M121",  
        "M122",  "M124",  "M126",  "M129", "M130",  "M133",  "M134",  "M135",  "M136",  
        "M137",  "M140",  "M141",  "M142", "M143",  "M145",  "M146",  "M148",  "M149",  
        "M151",  "M152",  "M153",  "M154", "M155",  "M160",  "M161",  "P216",  "P240")  
BAL.phy = subset_samples(BAL.phy, Subject_ID %in% IDs)
load("ABX_freq_accumulate_n.RData")
meta = sample_data(BAL.phy)
ABX.freq.accumulate.T1 = ABX.freq.accumulate.n[,1:3]
ABX.freq.accumulate.T2 = ABX.freq.accumulate.n[,c(1:2,4)]
ABX.freq.accumulate.PFT.T1 = ABX.freq.accumulate.n[,c(1:2,5)]
ABX.freq.accumulate.PFT.T2 = ABX.freq.accumulate.n[,c(1:2,6)]

ABX.freq.accumulate.PFT.T1.wide = dcast(ABX.freq.accumulate.PFT.T1, Subject~Type, value.var = "n.PFT.T1")
ABX.freq.accumulate.PFT.T2.wide = dcast(ABX.freq.accumulate.PFT.T2, Subject~Type, value.var = "n.PFT.T2")

ABX.freq.accumulate.PFT.T1.wide$Subject = paste0(ABX.freq.accumulate.PFT.T1.wide$Subject, "BA2")
ABX.freq.accumulate.PFT.T2.wide$Subject = paste0(ABX.freq.accumulate.PFT.T2.wide$Subject, "BA7")
ABX.freq.accumulate.wide = rbind.data.frame(ABX.freq.accumulate.PFT.T1.wide, ABX.freq.accumulate.PFT.T2.wide)
ABX.freq.accumulate.wide[is.na(ABX.freq.accumulate.wide)] = 0
ABX.freq.accumulate.wide$Subject
meta = cbind.data.frame(meta, ABX.freq.accumulate.wide[match(meta$Sample_ID, ABX.freq.accumulate.wide$Subject),])

colnames(meta)[36:37] = c("PXStaphylococcus", "PXPseudomonas")

# USABA7 as reference level
# meta$Group = relevel(meta$Group, ref = "USABA7")

BAL.phy = phy_substitute_metadata(BAL.phy, meta)

genus.data = aggregate_taxa(BAL.phy, "Genus")
BAL.phy.genus.relabun <- transform_sample_counts(genus.data, function(OTU) OTU/sum(OTU) * 100)


aus = meta$Sample_ID[which(meta$Country == "AUS")]
usa = meta$Sample_ID[which(meta$Country == "USA")]
otu = otu_table(BAL.phy.genus.relabun)
otu = as.data.frame(otu)
otu.aus = otu[,which(colnames(otu) %in% aus)]
otu.usa = otu[,which(colnames(otu) %in% usa)]
which(rownames(otu) == "Pseudomonas")
rowMeans(otu.aus)[342] # 3.725
rowMeans(otu.usa)[342] # 1.99
                  

input.metadata = meta(BAL.phy.genus.relabun)
input.data = as(otu_table(BAL.phy.genus.relabun), "matrix")

fit_data_filter <- Maaslin2(
  input.data, input.metadata,'MaAsLin2_output_no_adjustment', transform = "AST",
  fixed_effects = c('Group'),
  random_effects = c('Subject_ID'),
  normalization = 'TSS',
  reference = 'Group,USABA2',
  standardize = FALSE)

# ### adjust for covariates 
fit_data_filter2 <- Maaslin2(
  input.data, input.metadata,'MaAsLin2_output_adjustment', transform = "AST",
  fixed_effects = c('Group', "PE", "AZITH",           
                    "PXPseudomonas",
                    "ERAD", "Other"),
  random_effects = c('Subject_ID'),
  normalization = 'TSS',
  reference = 'Group,USABA2',
  standardize = FALSE)

# ### adjust for covariates 
fit_data_filter2 <- Maaslin2(
  input.data, input.metadata,'MaAsLin2_output_adjustment1', transform = "AST",
  fixed_effects = c('Group', "Gender", "F508", "PE", "AZITH",           
                    "PXPseudomonas",
                    "ERAD", "Other"),
  random_effects = c('Subject_ID'),
  normalization = 'TSS',
  reference = 'Group,USABA2',
  standardize = FALSE)

# ### adjust for covariates 
fit_data_filter2 <- Maaslin2(
  input.data, input.metadata,'MaAsLin2_output_adjustment2', transform = "AST",
  fixed_effects = c('Group', "Gender", "F508", "PE", "AZITH",           
                    "ERAD", "Other"),
  random_effects = c('Subject_ID'),
  normalization = 'TSS',
  reference = 'Group,USABA2',
  standardize = FALSE)



##### Phylum 
phylum.data = aggregate_taxa(BAL.phy, "Phylum")
BAL.phy.phylum.relabun <- transform_sample_counts(phylum.data, function(OTU) OTU/sum(OTU) * 100)


input.metadata.phylum = meta(BAL.phy.phylum.relabun)
input.data.phylum = as(otu_table(BAL.phy.phylum.relabun), "matrix")

fit_data_filter <- Maaslin2(
  input.data.phylum, input.metadata.phylum,
  'MaAsLin2_output_no_adjustment_phylum', transform = "AST",
  fixed_effects = c('Group'),
  random_effects = c('Subject_ID'),
  normalization = 'TSS',
  reference = 'Group,USABA2',
  standardize = FALSE)

# ### adjust for covariates 
fit_data_filter2 <- Maaslin2(
  input.data.phylum, input.metadata.phylum,
  'MaAsLin2_output_adjustment_phylum', transform = "AST",
  fixed_effects = c('Group', "PE", "AZITH",           
                    "PXPseudomonas",
                    "ERAD", "Other"),
  random_effects = c('Subject_ID'),
  normalization = 'TSS',
  reference = 'Group,USABA2',
  standardize = FALSE)

# ### adjust for covariates 
fit_data_filter2 <- Maaslin2(
  input.data.phylum, input.metadata.phylum,
  'MaAsLin2_output_adjustment_phylum1', transform = "AST",
  fixed_effects = c('Group', "Gender", "F508", "PE", "AZITH",           
                    "PXPseudomonas",
                    "ERAD", "Other"),
  random_effects = c('Subject_ID'),
  normalization = 'TSS',
  reference = 'Group,USABA2',
  standardize = FALSE)

# ### adjust for covariates 
fit_data_filter2 <- Maaslin2(
  input.data.phylum, input.metadata.phylum,
  'MaAsLin2_output_adjustment_phylum2', transform = "AST",
  fixed_effects = c('Group', "Gender", "F508", "PE", "AZITH",           
                    "ERAD", "Other"),
  random_effects = c('Subject_ID'),
  normalization = 'TSS',
  reference = 'Group,USABA2',
  standardize = FALSE)

