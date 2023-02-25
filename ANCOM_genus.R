library(ANCOMBC)
library(microbiome)
library(DT)
library(tibble)
library(grid)
library(maditr)
library(gautils2)
library(S4Vectors)
library(TreeSummarizedExperiment)
library(mia)
library(lme4)


BAL.phy = readRDS("../Data/BAL_phy.rds")
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

genus_data = aggregate_taxa(BAL.phy, "Genus")

tse = makeTreeSummarizedExperimentFromPhyloseq(BAL.phy)


### ancombc2 
## tutorial http://127.0.0.1:31970/library/ANCOMBC/doc/ANCOMBC2.html
set.seed(123)
output = ancombc2(data = tse, assay_name = "counts", tax_level = "Genus",
                  fix_formula = "Group", rand_formula = NULL,
                  p_adj_method = "holm", pseudo = 0, pseudo_sens = TRUE,
                  prv_cut = 0.10, lib_cut = 1000, s0_perc = 0.05,
                  group = "Group", struc_zero = TRUE, neg_lb = TRUE,
                  alpha = 0.05, n_cl = 2, verbose = TRUE,
                  global = TRUE, pairwise = FALSE, dunnet = TRUE, trend = FALSE,
                  iter_control = list(tol = 1e-2, max_iter = 20, 
                                      verbose = TRUE),
                  em_control = list(tol = 1e-5, max_iter = 100),
                  lme_control = lmerControl(),
                  mdfdr_control = list(fwer_ctrl_method = "holm", B = 100))


tab_zero = output$zero_ind
tab_zero %>%
  datatable(caption = "The detection of structural zeros")

tab_sens = output$pseudo_sens_tab
tab_sens %>%
  datatable(caption = "Sensitivity Scores") %>%
  formatRound(colnames(tab_sens), digits = 2)

res_prim = output$res
res_prim$taxon = rownames(res_prim)
df_GroupAUSBA2 = res_prim %>%
  dplyr::select(taxon, ends_with("GroupAUSBA2")) 
df_fig_GroupAUSBA2 = df_GroupAUSBA2 %>%
  filter(diff_GroupAUSBA2 == 1) %>% 
  arrange(desc(lfc_GroupAUSBA2)) %>%
  mutate(direct = ifelse(lfc_GroupAUSBA2 > 0, "Positive LFC", "Negative LFC"))
df_fig_GroupAUSBA2$taxon = factor(df_fig_GroupAUSBA2$taxon, levels = df_fig_GroupAUSBA2$taxon)
df_fig_GroupAUSBA2$direct = factor(df_fig_GroupAUSBA2$direct, 
                           levels = c("Positive LFC", "Negative LFC"))

fig_GroupAUSBA2 = df_fig_GroupAUSBA2 %>%
  ggplot(aes(x = taxon, y = lfc_GroupAUSBA2, fill = direct)) + 
  geom_bar(stat = "identity", width = 0.7, color = "black", 
           position = position_dodge(width = 0.4)) +
  geom_errorbar(aes(ymin = lfc_GroupAUSBA2 - se_GroupAUSBA2, ymax = lfc_GroupAUSBA2 + se_GroupAUSBA2), 
                width = 0.2, position = position_dodge(0.05), color = "black") + 
  labs(x = NULL, y = "Log fold change", 
       title = "Log fold changes AUS BA2 vs. USA BA2") + 
  scale_fill_discrete(name = NULL) +
  scale_color_discrete(name = NULL) +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(angle = 60, hjust = 1))
fig_GroupAUSBA2


### USA BA7
df_GroupUSABA7 = res_prim %>%
  dplyr::select(taxon, ends_with("GroupUSABA7")) 
df_fig_GroupUSABA7 = df_GroupUSABA7 %>%
  filter(diff_GroupUSABA7 == 1) %>% 
  arrange(desc(lfc_GroupUSABA7)) %>%
  mutate(direct = ifelse(lfc_GroupUSABA7 > 0, "Positive LFC", "Negative LFC"))
df_fig_GroupUSABA7$taxon = factor(df_fig_GroupUSABA7$taxon, levels = df_fig_GroupUSABA7$taxon)
df_fig_GroupUSABA7$direct = factor(df_fig_GroupUSABA7$direct, 
                                   levels = c("Positive LFC", "Negative LFC"))

fig_GroupUSABA7 = df_fig_GroupUSABA7 %>%
  ggplot(aes(x = taxon, y = lfc_GroupUSABA7, fill = direct)) + 
  geom_bar(stat = "identity", width = 0.7, color = "black", 
           position = position_dodge(width = 0.4)) +
  geom_errorbar(aes(ymin = lfc_GroupUSABA7 - se_GroupUSABA7, ymax = lfc_GroupUSABA7 + se_GroupUSABA7), 
                width = 0.2, position = position_dodge(0.05), color = "black") + 
  labs(x = NULL, y = "Log fold change", 
       title = "Log fold changes AUS BA2 vs. USA BA2") + 
  scale_fill_discrete(name = NULL) +
  scale_color_discrete(name = NULL) +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(angle = 60, hjust = 1))
fig_GroupUSABA7


### AUS BA7
df_GroupAUSBA7 = res_prim %>%
  dplyr::select(taxon, ends_with("GroupAUSBA7")) 
df_fig_GroupAUSBA7 = df_GroupAUSBA7 %>%
  filter(diff_GroupAUSBA7 == 1) %>% 
  arrange(desc(lfc_GroupAUSBA7)) %>%
  mutate(direct = ifelse(lfc_GroupAUSBA7 > 0, "Positive LFC", "Negative LFC"))
df_fig_GroupAUSBA7$taxon = factor(df_fig_GroupAUSBA7$taxon, levels = df_fig_GroupAUSBA7$taxon)
df_fig_GroupAUSBA7$direct = factor(df_fig_GroupAUSBA7$direct, 
                                   levels = c("Positive LFC", "Negative LFC"))

fig_GroupAUSBA7 = df_fig_GroupAUSBA7 %>%
  ggplot(aes(x = taxon, y = lfc_GroupAUSBA7, fill = direct)) + 
  geom_bar(stat = "identity", width = 0.7, color = "black", 
           position = position_dodge(width = 0.4)) +
  geom_errorbar(aes(ymin = lfc_GroupAUSBA7 - se_GroupAUSBA7, ymax = lfc_GroupAUSBA7 + se_GroupAUSBA7), 
                width = 0.2, position = position_dodge(0.05), color = "black") + 
  labs(x = NULL, y = "Log fold change", 
       title = "Log fold changes AUS BA2 vs. USA BA2") + 
  scale_fill_discrete(name = NULL) +
  scale_color_discrete(name = NULL) +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(angle = 60, hjust = 1))
fig_GroupAUSBA7

beta_val = beta_val * (q_val < 0.05) 
# Choose the maximum of beta's as the effect size
beta_pos = apply(abs(beta_val), 2, which.max) 
beta_max = vapply(seq_along(beta_pos), function(i) 
  beta_val[beta_pos[i], i], FUN.VALUE = double(1))
# Number of taxa except structural zeros
n_taxa = ifelse(is.null(out$zero_ind), 
                nrow(feature_table), 
                sum(apply(out$zero_ind, 1, sum) == 0))
# Cutoff values for declaring differentially abundant taxa
cut_off = 0.7 * (n_taxa - 1)

df_fig_w = res %>%
  dplyr::mutate(beta = beta_max,
                direct = case_when(
                  detected_0.7 == TRUE & beta > 0 ~ "Positive",
                  detected_0.7 == TRUE & beta <= 0 ~ "Negative",
                  TRUE ~ "Not Significant")) %>%
  dplyr::arrange(W)
df_fig_w$taxon_id = factor(df_fig_w$taxon_id, levels = df_fig_w$taxon_id)
df_fig_w$W = replace(df_fig_w$W, is.infinite(df_fig_w$W), n_taxa - 1)
df_fig_w$direct = factor(df_fig_w$direct, 
                         levels = c("Negative", "Positive", "Not Significant"))

p_w = df_fig_w %>%
  ggplot(aes(x = taxon_id, y = W, color = direct)) +
  geom_point(size = 2, alpha = 0.6) +
  labs(x = "Taxon", y = "W") +
  scale_color_discrete(name = NULL) + 
  geom_hline(yintercept = cut_off, linetype = "dotted", 
             color = "blue", size = 1.5) +
  geom_text(aes(x = 2, y = cut_off + 0.5, label = "W[0.7]"), 
            size = 5, vjust = -0.5, hjust = 0, color = "orange", parse = TRUE) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.major = element_blank())
p_w

## LFC
tab_lfc = res$lfc
colnames(tab_lfc)
col_name = c("USABA2 - USABA7", "AUSBA2 - USABA7", "AUSBA7 - USABA7")
colnames(tab_lfc) = col_name
tab_lfc %>% datatable(caption = "Log Fold Changes from the Primary Result") %>% formatRound(col_name, digits = 2)

## SE
tab_se = res$se
colnames(tab_se) = col_name
tab_se %>% datatable(caption = "SEs from the Primary Result") %>% formatRound(col_name, digits = 2)

## test statistic
tab_w = res$W
colnames(tab_w) = col_name
tab_w %>% datatable(caption = "Test Statistics from the Primary Result") %>% formatRound(col_name, digits = 2)

## adjusted p-value
tab_q = res$q
colnames(tab_q) = col_name
tab_q %>% datatable(caption = "Adjusted p-values from the Primary Result") %>% formatRound(col_name, digits = 2)

##Differentially abundant taxa
tab_diff = res$diff_abn
colnames(tab_diff) = col_name
tab_diff %>% datatable(caption = "Differentially Abundant Taxa from the Primary Result")

df_lfc = data.frame(tab_lfc * tab_diff, check.names = FALSE) %>% 
  rownames_to_column("taxon_id")
df_se = data.frame(tab_se * tab_diff, check.names = FALSE) %>% 
  rownames_to_column("taxon_id")
colnames(df_se)[-1] = paste0(colnames(df_se)[-1], "SE")

#################### "AUSBA7 - USABA7" ####################
df_fig_AUSBA7_USABA7 = 
  df_lfc %>% 
  dplyr::left_join(df_se, by = "taxon_id") %>%
  dplyr::transmute(taxon_id, `AUSBA7 - USABA7` , `AUSBA7 - USABA7SE`)%>%
  dplyr::filter(`AUSBA7 - USABA7` != 0) %>% 
  dplyr::arrange(desc(`AUSBA7 - USABA7`)) %>%
  dplyr::mutate(direct = ifelse(`AUSBA7 - USABA7` > 0, "Positive LFC", "Negative LFC"))
df_fig_AUSBA7_USABA7$taxon_id = factor(df_fig_AUSBA7_USABA7$taxon_id, 
                                       levels = df_fig_AUSBA7_USABA7$taxon_id)
df_fig_AUSBA7_USABA7$direct = factor(df_fig_AUSBA7_USABA7$direct, 
                                     levels = c("Positive LFC", "Negative LFC"))


#####################################################################
## adjust for gender and genotype of F508
out.adj = ancombc(phyloseq = genus_data, formula = "Group + Gender + 
                  Genotype_F508 + PE + AZITH + PXStaphylococcus + 
                  PXPseudomonas + ERAD + Other", 
                  p_adj_method = "holm", prv_cut = 0.10, lib_cut = 1000, 
                  group = "Group", struc_zero = TRUE, neg_lb = TRUE, tol = 1e-5, 
                  max_iter = 100, conserve = TRUE, alpha = 0.05, global = TRUE)

res.adj = out.adj$res
res.adj_global = out.adj$res_global

## LFC
tab.adj_lfc = res.adj$lfc
colnames(tab.adj_lfc)
col_name.adj = c("USABA2 - USABA7", "AUSBA2 - USABA7", 
                 "AUSBA7 - USABA7", "Male - Female", 
                 "Homozygous - Heterozygous",
                 "Other - Heterozygous",
                 "PE", "AZITH",
                 "PXStaphylococcus", "PXPseudomonas","ERAD", "Other")
colnames(tab.adj_lfc) = col_name.adj
tab.adj_lfc %>% datatable(caption = "Log Fold Changes from the Primary Result") %>% formatRound(col_name.adj, digits = 2)

## SE
tab.adj_se = res.adj$se
colnames(tab.adj_se) = col_name.adj
tab.adj_se %>% datatable(caption = "SEs from the Primary Result") %>% formatRound(col_name.adj, digits = 2)

## test statistic
tab.adj_w = res.adj$W
colnames(tab.adj_w) = col_name.adj
tab.adj_w %>% datatable(caption = "Test Statistics from the Primary Result") %>% formatRound(col_name.adj, digits = 2)

## adjusted p-value
tab.adj_q = res.adj$q
colnames(tab.adj_q) = col_name.adj
tab.adj_q %>% datatable(caption = "Adjusted p-values from the Primary Result") %>% formatRound(col_name.adj, digits = 2)

##Differentially abundant taxa
tab.adj_diff = res.adj$diff_abn
colnames(tab.adj_diff) = col_name.adj
tab.adj_diff %>% datatable(caption = "Differentially Abundant Taxa from the Primary Result")

df.adj_lfc = data.frame(tab.adj_lfc * tab.adj_diff, check.names = FALSE) %>% 
  rownames_to_column("taxon_id")
df.adj_se = data.frame(tab.adj_se * tab.adj_diff, check.names = FALSE) %>% 
  rownames_to_column("taxon_id")
colnames(df.adj_se)[-1] = paste0(colnames(df.adj_se)[-1], "SE")

df_fig_AUSBA7_USABA7.adj = 
  df.adj_lfc %>% 
  dplyr::left_join(df.adj_se, by = "taxon_id") %>%
  dplyr::transmute(taxon_id, `AUSBA7 - USABA7` , `AUSBA7 - USABA7SE`)%>%
  dplyr::filter(`AUSBA7 - USABA7` != 0) %>% 
  dplyr::arrange(desc(`AUSBA7 - USABA7`)) %>%
  dplyr::mutate(direct = ifelse(`AUSBA7 - USABA7` > 0, "Positive LFC", "Negative LFC"))
df_fig_AUSBA7_USABA7.adj$taxon_id = factor(df_fig_AUSBA7_USABA7.adj$taxon_id, 
                                           levels = df_fig_AUSBA7_USABA7.adj$taxon_id)
df_fig_AUSBA7_USABA7.adj$direct = factor(df_fig_AUSBA7_USABA7.adj$direct, 
                                         levels = c("Positive LFC", "Negative LFC"))


df_fig_AUSBA7_USABA7.all = rbind.data.frame(df_fig_AUSBA7_USABA7, 
                                            df_fig_AUSBA7_USABA7.adj)
df_fig_AUSBA7_USABA7.all$type = c(rep("crude", dim(df_fig_AUSBA7_USABA7)[1]), 
                                  rep("adjusted", dim(df_fig_AUSBA7_USABA7.adj)[1]))
colnames(df_fig_AUSBA7_USABA7.all) = c("taxon_id", "AUSBA7 - USABA7", 
                                       "AUSBA7 - USABA7SE", "direct", 
                                       "model_type")
text_left <- textGrob("Increase in USA", gp=gpar(fontsize=10))
text_right <- textGrob("Increase in AUS", gp=gpar(fontsize=10))
text_xlab <- textGrob("Log fold change", gp=gpar(fontsize=14))

# crude yellow; adjusted model blue
# png("Figures/AUSBA7_vs_USABA7_ANCOMBC_phyla.png", height = 5, width = 8, units = "in", res = 300)
ggplot(data = df_fig_AUSBA7_USABA7.all, 
       aes(y = taxon_id, x = `AUSBA7 - USABA7`, 
           color = `model_type`, shape = `model_type`)) + 
  geom_pointrange(aes(xmin=`AUSBA7 - USABA7`- `AUSBA7 - USABA7SE`, 
                      xmax=`AUSBA7 - USABA7`+ `AUSBA7 - USABA7SE`),
                  position = position_dodge2(.5)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  ggtitle("AUS vs. USA at T2") +
  ylab(" ") +
  xlab(" ") +
  scale_color_manual(values=c("#0072B2", "#E69F00")) +
  xlim(c(-6,5)) + 
  theme_classic() + 
  theme(#legend.position = "none",
        plot.margin = unit(c(4,1,2,1), "lines"),
        plot.title = element_text(hjust = 0.75, face = "bold", size = 16),
        axis.title.x = element_text(vjust=-4)) +
  geom_segment(aes(x =-0.1, xend = -.9, y = -0.8, yend = -0.8),
               arrow = arrow(length = unit(0.2, "cm")), color = "black") +
  geom_segment(aes(x =0.1, xend = .9, y = -0.8, yend = -0.8),
               arrow = arrow(length = unit(0.2, "cm")), color = "black") +
  annotation_custom(text_left,xmin=-3.5,xmax=-1.5,ymin=-0.8,ymax=-0.8) + 
  annotation_custom(text_right,xmin=1.5,xmax=3.5,ymin=-0.8,ymax=-0.8) + 
  annotation_custom(text_xlab,xmin=-1,xmax=1,ymin=-1,ymax=-1) + 
  coord_cartesian(ylim = c(0, nlevels(unique(df_fig_AUSBA7_USABA7.all$taxon_id))), clip = "off")
# dev.off() 
#########################################################################

