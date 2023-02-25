library(rstudioapi)
library(dplyr)
library(tidyr)
library(glmmLasso)
library(reshape2)

setwd(dirname(getActiveDocumentContext()$path))
load(file = "variable_selection.RData")
load(file = "../ABX_freq_accumulate_n.RData")
# load(file = "../ABX_freq.RData")

ABX.freq.accumulate.n.PFT1 = ABX.freq.accumulate.n %>% 
  dplyr::select(Subject, Type, n.PFT.T1)
ABX.freq.accumulate.n.PFT2 = ABX.freq.accumulate.n %>% 
  dplyr::select(Subject, Type, n.PFT.T2)
ABX.freq.accumulate.n.wide.T1 = dcast(ABX.freq.accumulate.n.PFT1, 
                                      Subject ~ Type, value.var="n.PFT.T1")
ABX.freq.accumulate.n.wide.T2 = dcast(ABX.freq.accumulate.n.PFT2, 
                                      Subject ~ Type, value.var="n.PFT.T2")
ABX.freq.accumulate.n.wide.T1$Subject = paste0(ABX.freq.accumulate.n.wide.T1$Subject,
                                               "BA2")
ABX.freq.accumulate.n.wide.T2$Subject = paste0(ABX.freq.accumulate.n.wide.T2$Subject,
                                               "BA7")

ABX.freq.accumulate.n.wide.T1[is.na(ABX.freq.accumulate.n.wide.T1)] = 0
ABX.freq.accumulate.n.wide.T2[is.na(ABX.freq.accumulate.n.wide.T2)] = 0
ABX.freq.accumulate.n.wide = rbind.data.frame(ABX.freq.accumulate.n.wide.T1,
                                              ABX.freq.accumulate.n.wide.T2)
colnames(ABX.freq.accumulate.n.wide)[1] = "Sample_ID"
OTU = as.data.frame(OTU)
dim(OTU)
OTU.adj = OTU + 1
## normalize OTU
z_OTU.adj = log(OTU.adj)
clrx_OTU <- apply(z_OTU.adj, 2, function(x) x - rowMeans(z_OTU.adj))

OTU = clrx_OTU
OTU = as.data.frame(OTU)
OTU$Sample_ID = rownames(OTU)
clinic = OTU %>% inner_join(cytokine, by = "Sample_ID") %>%
  inner_join(PFT, by = "Sample_ID") 

clinic = clinic[,-grep(".y", colnames(clinic), fixed = T)]
colnames(clinic)[grep(".x", colnames(clinic), fixed = T)] = str_remove(colnames(clinic)[grep(".x", colnames(clinic), fixed = T)],
                                                                       ".x")

#### select taxa that are associated with fev0.5
clinic$SampleID = factor(clinic$SampleID)
bics = rep(NA, 300)
ks = c(seq(0,1,0.01), seq(2,200,1))
## use the taxa that are significantly associated with FEV0.5 with p-value < 0.05
for (k in 1:300){
  lm = glmmLasso(fix=Pfev0.5~Streptococcus + Gemella + Granulicatella +
                   Rothia + `Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium` +
                   Hephaestia + Novosphingobium +
                      Cupriavidus + Bergeyella + 
                   Enhydrobacter, 
                  rnd=list(SampleID=~1), data = clinic,
                  lambda = ks[k],
                  switch.NR=T, final.re=T, control = list()) 
  bics[k] = lm$bic
}
names(bics) = c(seq(0,1,0.01), seq(2,200,1))
which.min(bics) # lambda = 0.97
lm0= glmmLasso(fix=Pfev0.5~Streptococcus + Gemella + Granulicatella +
                 Rothia + `Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium` +
                 Hephaestia + Novosphingobium +#Actinomyces +
                 Cupriavidus + Bergeyella + #Sphingomonas +
                 Enhydrobacter, 
               rnd=list(SampleID=~1), data = clinic,
               lambda = 0.97,
               switch.NR=T, final.re=T, control = list())
temp1 = lm0$coefficients
