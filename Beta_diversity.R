library(dplyr)
library(ggplot2)
library(readxl)
library(wesanderson)
library(phyloseq)
library(gautils2)

dirname = dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(dirname)
phy = readRDS("../Data/phy_noncontam.rds")
meta = meta(phy)
tax = tax_table(phy)
otu = otu_table(phy)

# meta$Visit = ifelse(meta$Visit == "BA2", "T1", 
#                     ifelse(meta$Visit == "BA7", "T2", ""))
BAL.BA2 = meta %>% dplyr::filter(Sample_type == "BAL" & Visit == "BA2") #51
BAL.BA7 = meta %>% dplyr::filter(Sample_type == "BAL" & Visit == "BA7") #58
dim(meta %>% dplyr::filter(Study == "ECF" & Sample_type == "BAL" & Country == "AUS" ))
dim(meta %>% dplyr::filter(Study == "ECF" & Sample_type == "BAL" & Country == "USA" ))
meta$Group = paste0(meta$Country, meta$Visit)
meta$Group = factor(meta$Group, 
                       levels = c("USABA2", "AUSBA2",
                                  "USABA7", "AUSBA7"))
phy = phy_substitute_metadata(phy, meta)

BA2_BA7_ID = intersect(BAL.BA2$Subject_ID, BAL.BA7$Subject_ID) 
# 45; we keep these 45 matched BAL samples. 
# Now we find the matched prewash samples 

phy = subset_samples(phy, Subject_ID %in% BA2_BA7_ID)
phy = subset_samples(phy, Sample_type == "BAL")
meta = sample_data(phy)
random_tree <- rtree(ntaxa(phy), rooted = TRUE, tip.label = taxa_names(phy))
physeq = merge_phyloseq(phy, random_tree)

# tutorial: 
#https://david-barnett.github.io/evomics-material/exercises/microViz-2-dissimilarity-exercises#PERMANOVA:
physeq %>% 
  ps_mutate(reads = sample_sums(physeq)) %>% 
  samdat_tbl() %>% 
  ggplot(aes(x = reads)) + 
  geom_freqpoly(bins = 500) +
  geom_rug(alpha = 0.5) +
  scale_x_log10(labels = scales::label_number()) +
  labs(x = 'Number of classified reads', y = NULL) +
  theme_bw()

## Proportional min_prevalence given: 0.01 --> min 33/1644 samples.
physeq %>% 
  tax_filter(min_prevalence = 1 / 100, min_total_abundance = 1000) %>% 
  ntaxa()

physeqTaxaStats <- tibble(
  taxon = taxa_names(physeq),
  prevalence = microbiome::prevalence(physeq),
  total_abundance = taxa_sums(physeq)
)
TAX = as.data.frame(as.matrix(tax_table(physeq)))
physeqTaxaStats$taxon = TAX$Genus[match(physeqTaxaStats$taxon, rownames(TAX))]

p <- physeqTaxaStats %>%
  ggplot(aes(total_abundance, prevalence)) +
  geom_point(alpha = 0.5) +
  geom_rug(alpha = 0.1) +
  scale_x_continuous(
    labels = scales::label_number(), name = "Total Abundance"
  ) +
  scale_y_continuous(
    labels = scales::label_percent(), breaks = scales::breaks_pretty(n = 9),
    name = "Prevalence (%)",
    sec.axis = sec_axis(
      trans = ~ . * nsamples(physeq), breaks = scales::breaks_pretty(n = 9),
      name = "Prevalence (N samples)"
    )
  ) +
  theme_bw()
p

p + ggrepel::geom_text_repel(
  data = function(df) filter(df, total_abundance > 1e4 | prevalence > 0.2),
  mapping = aes(label = taxon), size = 2.5, min.segment.length = 0, force = 15
)

#Now let’s zoom in on the less abundant taxa by log-transforming the axes. 
#We’ll also add lines indicating the thresholds of 2% prevalence and 10000 
#reads abundance.

physeqTaxaStats %>%
  ggplot(aes(x = total_abundance, y = prevalence)) +
  geom_vline(xintercept = 1000, color = "red", linetype = "dotted") +
  geom_hline(yintercept = 1 / 100, color = "red", linetype = "dotted") +
  geom_point(alpha = 0.5) +
  geom_rug(alpha = 0.1) +
  scale_x_log10(labels = scales::label_number(), name = "Total Abundance") +
  scale_y_log10(
    labels = scales::label_percent(), breaks = scales::breaks_log(n = 9),
    name = "Prevalence (%)",
    sec.axis = sec_axis(
      trans = ~ . * nsamples(shao19), breaks = scales::breaks_log(n = 9),
      name = "Prevalence (N samples)"
    )
  ) +
  theme_bw()

# jaccard
physeq %>%
  tax_filter(min_prevalence = 1 / 100) %>%
  tax_fix() %>%
  tax_agg(rank = "Genus") %>%
  tax_transform("binary") %>% # converts counts to absence/presence: 0/1
  dist_calc(dist = "jaccard")
distances <- physeq %>%
  tax_filter(min_prevalence = 1 / 100, verbose = FALSE) %>%
  tax_fix() %>%
  tax_agg(rank = "Genus") %>%
  tax_transform("binary") %>%
  dist_calc(dist = "jaccard") %>%
  dist_get()
as.matrix(distances)[1:5, 1:5]
range(distances)

physeq %>%
  tax_filter(min_prevalence = 1 / 100, verbose = FALSE) %>%
  tax_fix() %>% 
  # ps_calc_dominant(rank = "Genus", none = "Mixed", other = "Other") %>%
  tax_agg(rank = "Genus") %>%
  dist_calc(dist = "bray") %>%
  ord_calc(method = "PCoA") %>%
  ord_plot(color = "Group", alpha = 0.5, size = 2) +
  scale_color_manual(values = c("#ABDDA4", "#D7191C", "#FDAE61", "#2B83BA")) + 
  theme_classic(12) +
  coord_fixed(0.7)

# # shiny
# physeq %>%
#   tax_filter(min_prevalence = 1 / 100, verbose = FALSE) %>%
#   tax_fix() %>% 
#   # calculate new sample variables with dominant taxon (optional)
#   # ps_calc_dominant(rank = "Genus", none = "Mixed", other = "Other") %>%
#   # launch a Shiny app in your web browser!
#   ord_explore()

# png("Figures/bray.png", width = 10, height = 10, units = "in", res=300)
physeq %>%
  tax_filter(min_prevalence = 1 / 100, verbose = FALSE) %>%
  tax_fix() %>% 
  tax_agg(rank = "Genus") %>%
  dist_calc(dist = "bray") %>%
  ord_calc(method = "PCoA") %>%
  ord_plot(alpha = 0.5, size = 2, color = "Group") +
  theme_classic() +
  ylab("PCo2 [15.8%]") +
  xlab("PCo1 [19.2%]") + 
  ggtitle("Bray-Curtis") + 
  coord_fixed(0.7) +
  stat_ellipse(aes(color = Group), size = 2) +
  scale_color_manual(values = c("#ABDDA4", "#D7191C", "#FDAE61", "#2B83BA")) + 
  theme(axis.title = element_text(size = 22, face = "bold"),
        axis.text = element_text(size = 20),
        plot.title = element_text(size = 20, face = "bold"),
        legend.title = element_blank(),
        legend.text = element_text(size = 20, face = "bold")) 
# dev.off()


physeq %>%
  tax_filter(min_prevalence = 1 / 100, verbose = FALSE) %>%
  tax_fix() %>% 
  tax_agg(rank = "Genus") %>%
  dist_calc(dist = "bray") %>%
  dist_permanova(variables = "Group", n_perms = 999, seed = 123) %>%
  perm_get()
# Group     3   2.9376 0.11258 3.6365  0.001 ***
phy.genus = physeq %>%
  tax_filter(min_prevalence = 1 / 100, verbose = FALSE) %>%
  tax_fix() %>% 
  tax_agg(rank = "Genus")
otu = otu_get(phy.genus)
set.seed(123)
adonis2(otu~Group,data=meta(meta), permutations=999, method="bray")
pairwise.adonis2(otu~Group, data = meta(meta))
# $AUSBA2_vs_AUSBA7
# Df SumOfSqs      R2      F Pr(>F)    
# Group     1    1.438 0.08811 5.6044  0.001 ***
# $AUSBA2_vs_USABA2
# Df SumOfSqs      R2      F Pr(>F)    
# Group     1   1.4926 0.12608 6.2036  0.001 ***
# $AUSBA2_vs_USABA7
# Df SumOfSqs      R2      F Pr(>F)    
# Group     1   1.5489 0.12853 6.3419  0.001 ***
dist = phy.genus %>%
  dist_calc(dist = "bray") %>%  
  dist_get()
permutest(betadisper(dist, meta$Group), pairwise = TRUE) # 0.211



### unifrac 
# png("Figures/unifrac.png", width = 10, height = 10, units = "in", res=300)
physeq %>%
  tax_filter(min_prevalence = 1 / 100, verbose = FALSE) %>%
  tax_fix() %>% 
  # tax_agg(rank = "Genus") %>%
  dist_calc(dist = "unifrac") %>% 
  ord_calc(method = "PCoA") %>%
  ord_plot(alpha = 0.5, size = 2, color = "Group") +
  theme_classic() +
  ylab("PCo2 [5.2%]") +
  xlab("PCo1 [11.3%]") +
  ggtitle("UniFrac") + 
  coord_fixed(0.7) +
  stat_ellipse(aes(color = Group), size = 2) +
  scale_color_manual(values = c("#ABDDA4", "#D7191C", "#FDAE61", "#2B83BA")) + 
  theme(axis.title = element_text(size = 22, face = "bold"),
        axis.text = element_text(size = 20),
        plot.title = element_text(size = 20, face = "bold"),
        legend.title = element_blank(),
        legend.text = element_text(size = 20, face = "bold")) 
# dev.off()

physeq %>%
  tax_filter(min_prevalence = 1 / 100, verbose = FALSE) %>%
  tax_fix() %>% 
  dist_calc(dist = "unifrac") %>%
  dist_permanova(variables = "Group", n_perms = 999, seed = 123) %>%
  perm_get()
# Group     3   1.6793 0.06378 1.9529  0.001 ***
phy.genus = physeq %>%
  tax_filter(min_prevalence = 1 / 100, verbose = FALSE) %>%
  tax_fix() 
# otu = otu_get(phy.genus)
set.seed(123)
dist = phy.genus %>%
  dist_calc(dist = "unifrac") %>%  
  dist_get()
dist.unifrac = dist
adonis2(dist~Group,data=meta(meta), permutations=999)
pairwise.adonis2(dist~Group, data = meta(meta))
# $AUSBA2_vs_AUSBA7
# Df SumOfSqs      R2      F Pr(>F)    
# Group     1   0.8185 0.04656 2.8322  0.001 ***
#   $AUSBA2_vs_USABA2
# Df SumOfSqs      R2      F Pr(>F)    
# Group     1   0.7019 0.05323 2.4176  0.001 ***
#   $AUSBA2_vs_USABA7
# Df SumOfSqs      R2      F Pr(>F)    
# Group     1   0.6451 0.04812 2.1737  0.001 ***
#   $AUSBA7_vs_USABA2
# Df SumOfSqs      R2      F Pr(>F)  
# Group     1   0.3507 0.02866 1.2685  0.075 .
# $AUSBA7_vs_USABA7
# Df SumOfSqs      R2      F Pr(>F)  
# Group     1   0.3411 0.02727 1.2054  0.094 .

permutest(betadisper(dist, meta$Group), pairwise = TRUE) # 0.378


## weighted unifrac
# png("Figures/wunifrac.png", width = 8, height = 8, units = "in", res=300)
physeq %>%
  tax_filter(min_prevalence = 1 / 100, verbose = FALSE) %>%
  tax_fix() %>% 
  # tax_agg(rank = "Genus") %>%
  dist_calc(dist = "wunifrac") %>% 
  ord_calc(method = "PCoA") %>%
  ord_plot(alpha = 0.5, size = 2, color = "Group") +
  theme_classic() +
  ylab("PCo2 [13.8%]") +
  xlab("PCo1 [20.0%]") +
  coord_fixed(0.7) +
  stat_ellipse(aes(color = Group), size = 2) +
  scale_color_manual(values = c("#ABDDA4", "#D7191C", "#FDAE61", "#2B83BA")) + 
  theme(axis.title = element_text(size = 22, face = "bold"),
        axis.text = element_text(size = 20),
        plot.title = element_text(size = 20, face = "bold"),
        legend.title = element_blank(),
        legend.text = element_text(size = 20, face = "bold")) 
# dev.off()

physeq %>%
  tax_filter(min_prevalence = 1 / 100, verbose = FALSE) %>%
  tax_fix() %>% 
  dist_calc(dist = "wunifrac") %>%
  dist_permanova(variables = "Group", n_perms = 999, seed = 123) %>%
  perm_get()
# Group     3   1.3724 0.07189 2.2204  0.001 ***
phy.genus = physeq %>%
  tax_filter(min_prevalence = 1 / 100, verbose = FALSE) %>%
  tax_fix() 
set.seed(123)
dist = phy.genus %>%
  dist_calc(dist = "wunifrac") %>%  
  dist_get()
dist.wunifrac = dist
adonis2(dist.wunifrac~Group,data=meta(meta), permutations=999)
pairwise.adonis2(dist.wunifrac~Group, data = meta(meta))
# $AUSBA2_vs_AUSBA7
# Df SumOfSqs      R2      F Pr(>F)   
# Group     1   0.6514 0.04825 2.9405  0.003 **
#   $AUSBA2_vs_USABA2
# Df SumOfSqs      R2      F Pr(>F)   
# Group     1   0.6562 0.06995 3.2342  0.005 **
#   $AUSBA2_vs_USABA7
# Df SumOfSqs      R2      F Pr(>F)  
# Group     1   0.5182 0.05245 2.3803  0.013 *
#   $AUSBA7_vs_USABA2
# Df SumOfSqs      R2      F Pr(>F)  
# Group     1   0.2951 0.03411 1.5185   0.09 .
  
permutest(betadisper(dist.wunifrac, meta$Group), pairwise = TRUE) 
# 0.049 violated 


## generalized unifrac alpha = 0.5
# png("Figures/gunifrac.png", width = 8, height = 8, units = "in", res=300)
physeq %>%
  tax_filter(min_prevalence = 1 / 100, verbose = FALSE) %>%
  tax_fix() %>% 
  # tax_agg(rank = "Genus") %>%
  dist_calc(dist = "gunifrac") %>% 
  ord_calc(method = "PCoA") %>%
  ord_plot(alpha = 0.5, size = 2, color = "Group") +
  theme_classic() +
  ylab("PCo2 [7.5%]") +
  xlab("PCo1 [11.3%]") +
  coord_fixed(0.7) +
  stat_ellipse(aes(color = Group), size = 2) +
  scale_color_manual(values = c("#ABDDA4", "#D7191C", "#FDAE61", "#2B83BA")) + 
  theme(axis.title = element_text(size = 22, face = "bold"),
        axis.text = element_text(size = 20),
        plot.title = element_text(size = 20, face = "bold"),
        legend.title = element_blank(),
        legend.text = element_text(size = 20, face = "bold")) 
# dev.off()


physeq %>%
  tax_filter(min_prevalence = 1 / 100, verbose = FALSE) %>%
  tax_fix() %>% 
  dist_calc(dist = "gunifrac") %>%
  dist_permanova(variables = "Group", n_perms = 999, seed = 123) %>%
  perm_get()
# Group     3   1.8847 0.06897 2.1235  0.001 ***
phy.genus = physeq %>%
  tax_filter(min_prevalence = 1 / 100, verbose = FALSE) %>%
  tax_fix() 
set.seed(123)
dist = phy.genus %>%
  dist_calc(dist = "gunifrac") %>%  
  dist_get()
dist.gunifrac = dist
adonis2(dist.gunifrac~Group,data=meta(meta), permutations=999)
pairwise.adonis2(dist.gunifrac~Group, data = meta(meta))
# $AUSBA2_vs_AUSBA7
# Df SumOfSqs      R2      F Pr(>F)    
# Group     1   0.8635 0.04648 2.8275  0.001 ***
#   $AUSBA2_vs_USABA2
# Df SumOfSqs      R2      F Pr(>F)    
# Group     1   0.8653 0.06384 2.9322  0.001 ***
#   $AUSBA2_vs_USABA7
# Df SumOfSqs      R2      F Pr(>F)    
# Group     1   0.7208 0.05173 2.3459  0.001 ***
#   $AUSBA7_vs_USABA2
# Df SumOfSqs      R2     F Pr(>F)  
# Group     1   0.4178 0.03303 1.469  0.036 *
#   $AUSBA7_vs_USABA7
# Df SumOfSqs      R2      F Pr(>F)  
# Group     1   0.3864 0.02941 1.3028  0.083 .
#   


permutest(betadisper(dist.wunifrac, meta$Group), pairwise = TRUE) 
# 0.049 violated 
