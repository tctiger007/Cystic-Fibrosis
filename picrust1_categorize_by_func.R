### Reproducing the categorize by function (level 3) functionality in plain-text tables.
### Doing this because adding a column of KEGG Pathways to a table and then converting
### that table to BIOM is difficult.
### referencing picrust2 author

categorize_by_function_l3 <- function(in_ko, kegg_brite_mapping) {
  # Function to create identical output as categorize_by_function.py script,
  # but with R objects instead of BIOM objects in Python.
  # Input KO table is assumed to have rownames as KOs and sample names as columns.

  out_pathway <- data.frame(matrix(NA, nrow=0, ncol=(ncol(in_ko) + 1)))

  colnames(out_pathway) <- c("pathway", colnames(in_ko))

  for(ko in rownames(in_ko)) {
    
    # Skip KO if not in KEGG BRITE mapping df
    # (this occurs with newer KOs that weren't present in PICRUSt1).
    if(! ko %in% rownames(kegg_brite_mapping)) {
      next
    }
    
    pathway_list <- strsplit(kegg_brite_mapping[ko, "metadata_KEGG_Pathways"], "\\|")[[1]]
    
    for(pathway in pathway_list) {
      
      pathway <- strsplit(pathway, ";")[[1]][3]
      
      new_row <- data.frame(matrix(c(NA, as.numeric(in_ko[ko,])), nrow=1, ncol=ncol(out_pathway)))
      colnames(new_row) <- colnames(out_pathway)
      new_row$pathway <- pathway
      out_pathway <- rbind(out_pathway, new_row)
    }
    
  }
  
  out_pathway = data.frame(aggregate(. ~ pathway, data = out_pathway, FUN=sum))
  
  rownames(out_pathway) <- out_pathway$pathway
  
  out_pathway <- out_pathway[, -which(colnames(out_pathway) == "pathway")]
  
  if(length(which(rowSums(out_pathway) == 0)) > 0) {
    out_pathway <- out_pathway[-which(rowSums(out_pathway) == 0), ]
  }
  
  return(out_pathway)
  
}


### Example commands:
### Read in BRITE hierarchy per KO.
kegg_brite_map <- read.table("picrust1_KO_BRITE_map.tsv",
                             header=TRUE, sep="\t", quote = "", stringsAsFactors = FALSE, comment.char="", row.names=1)


## When reading in tab-delimited file of KO predictions (PICRUSt2 output):
test_ko <- read.table("picrust2_out/KO_metagenome_out/pred_metagenome_unstrat.tsv", 
                      header=TRUE, sep="\t", row.names=1)
# test_ec = read.table("picrust2_out/EC_metagenome_out/pred_metagenome_unstrat.tsv", 
#                      header=TRUE, sep="\t", row.names=1)
# test_pw = read.table("picrust2_out/pathways_out/path_abun_unstrat.tsv",
#                      header=TRUE, sep="\t", row.names=1)


## Run function to categorize all KOs by level 3 in BRITE hierarchy.
test_ko_L3 <- categorize_by_function_l3(test_ko, kegg_brite_map)
write.csv(test_ko_L3, "KO_L3.csv")

