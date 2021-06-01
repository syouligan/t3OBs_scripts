#!/usr/bin/Rscript

# --------------------------------------------------------------------------
# Venn diagram of overlap between differentially expressed genes in each group
# --------------------------------------------------------------------------

# Set adaptive working directory
if(dir.exists("/Users/mac/cloudstor/")) {
  setwd("/Users/mac/cloudstor/imperial_projects/THcircuit/") # Uses practice data (5% of cells from each sample) if running locally
  place <- "local"
} else {
  setwd("/share/ScratchGeneral/scoyou/imperial_projects/THcircuit/")
  place <- "wolfpack"
}

library('gplots')

# Load gene annotation info
geneActivity <- read.csv('project_results/t3OBs/Gene_annotation_w_activity.csv', row.names = 1, header = TRUE)

# Make venn diagram of each LFC threshold
# --------------------------------------------------------------------------

for(LFC in c(0)) {
# for(LFC in c(0, 0.15, 0.3, 0.45)) {
  
  # Make list of vectors of differentially expressed genes for each LFC threshold
  WTV1_WTT3 <- read.csv(paste0("project_results/t3OBs/DEGs_WTV1_WTT3_", LFC, ".csv"), header = TRUE)
  WTV1_WTT3_DGE <- WTV1_WTT3[WTV1_WTT3$adj.P.Val < 0.05 & (!is.na(WTV1_WTT3$GeneSymbol)), ]
  
  WTV1_ALPHAV1 <- read.csv(paste0("project_results/t3OBs/DEGs_WTV1_ALPHAV1_", LFC, ".csv"), header = TRUE)
  WTV1_ALPHAV1_DGE <- WTV1_ALPHAV1[WTV1_ALPHAV1$adj.P.Val < 0.05 & (!is.na(WTV1_ALPHAV1$GeneSymbol)), ]
  
  WTT3_ALPHAT3 <- read.csv(paste0("project_results/t3OBs/DEGs_WTT3_ALPHAT3_", LFC, ".csv"), header = TRUE)
  WTT3_ALPHAT3_DGE <- WTT3_ALPHAT3[WTT3_ALPHAT3$adj.P.Val < 0.05 & (!is.na(WTT3_ALPHAT3$GeneSymbol)), ]
  
  ALPHAV1_ALPHAT3 <- read.csv(paste0("project_results/t3OBs/DEGs_ALPHAV1_ALPHAT3_", LFC, ".csv"), header = TRUE)
  ALPHAV1_ALPHAT3_DGE <- ALPHAV1_ALPHAT3[ALPHAV1_ALPHAT3$adj.P.Val < 0.05 & (!is.na(ALPHAV1_ALPHAT3$GeneSymbol)), ]
  
  DGE_list <- list("WTV1_WTT3" = WTV1_WTT3_DGE$GeneSymbol, "WTV1_ALPHAV1" = WTV1_ALPHAV1_DGE$GeneSymbol, "WTT3_ALPHAT3" = WTT3_ALPHAT3_DGE$GeneSymbol, "ALPHAV1_ALPHAT3" = ALPHAV1_ALPHAT3_DGE$GeneSymbol)
  
  # Venn diagram overlap DEGs between conditions
  pdf(paste0("project_results/t3OBs/DGE_venn", LFC, ".pdf"))
  ItemsList <- venn(DGE_list, show.plot = TRUE)
  dev.off()
  
  lists <- attributes(ItemsList)$intersections
  print(lists)
  
  # Make gene lists of intersections
  for(int in names(lists)){
    genes_tmp <- lists[[int]]
    df_tmp <- geneActivity[geneActivity$GeneSymbol %in% genes_tmp, ]
    write.csv(df_tmp, paste0("project_results/t3OBs/Venn_partition_", int, "_LFC_", LFC, ".csv"))
    
    geneActivity[,paste0(int, "_", LFC)] <- geneActivity$GeneSymbol %in% genes_tmp
  }
  
  # Make a gene list of GOIs
  tmp_df <- geneActivity[, c(paste0("WTV1_WTT3:WTT3_ALPHAT3_", LFC), paste0("WTV1_WTT3_", LFC), paste0("WTV1_WTT3:WTV1_ALPHAV1_", LFC), paste0("WTV1_WTT3:WTV1_ALPHAV1:WTT3_ALPHAT3_", LFC))]
  GOIs <- geneActivity[rowSums(tmp_df) > 0, ]
  write.csv(df_tmp, paste0("project_results/t3OBs/GOIs_", LFC, ".csv"))

  geneActivity[,paste0("GOI_", LFC)] <- geneActivity$GeneSymbol %in% GOIs$GeneSymbol
  }

write.csv(geneActivity, "project_results/t3OBs/Gene_annotation_w_activity_DGE_Venn.csv")



