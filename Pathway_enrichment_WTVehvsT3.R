#!/usr/bin/Rscript

# --------------------------------------------------------------------------
#! Perform patwhay enrichment in WT Veh vs T3 treated osteoblasts
# --------------------------------------------------------------------------

# Set adaptive working directory
if(dir.exists("/Users/mac/cloudstor/")) {
  setwd("/Users/mac/cloudstor/imperial_projects/THcircuit/") # Uses practice data (5% of cells from each sample) if running locally
  place <- "local"
} else {
  setwd("/share/ScratchGeneral/scoyou/imperial_projects/THcircuit/")
  place <- "wolfpack"
}

# Libraries
library('dplyr')
library('tidyr')
library('ggplot2')
library('readr')
library('Matrix')
library('clusterProfiler')
library('ReactomePA')
library('org.Mm.eg.db')
library('gplots')
library('msigdbr')
library('gtools')
library('ggcorrplot')

# Load gene annotation info
geneActivity <- read.csv('project_results/t3OBs/Gene_annotation_w_activity_DGE_Venn.csv', row.names = 1, header = TRUE)

# set gene universe as any genes expressed in WT osteoblasts (Veh or T3 treated)
universe <- as.character(geneActivity[geneActivity$WT_any, "GeneSymbol"])
universe <- universe[!is.na(universe)]
universe_entrez <- as.character(geneActivity[geneActivity$WT_any, "EntrezID"])
universe_entrez <- universe_entrez[!is.na(universe_entrez)]

# Makes HPO geneset for enrichment testing
h_df <- msigdbr(species = "Mus musculus", category = "C5", subcategory = "HPO")
h_t2g <- h_df %>% dplyr::select(gs_name, gene_symbol) %>% as.data.frame()

# Makes TF geneset for enrichment testing
htf_df <- msigdbr(species = "Mus musculus", category = "C3", subcategory = 'TFT:GTRD')
htf_t2g <- htf_df %>% dplyr::select(gs_name, gene_symbol) %>% as.data.frame()

# Makes PID geneset for enrichment testing
hPID_df <- msigdbr(species = "Mus musculus", category = "C2", subcategory = "CP:PID")
hPID_t2g <- hPID_df %>% dplyr::select(gs_name, gene_symbol) %>% as.data.frame()

# Make list of genes of interest
DGE <- read.csv("~/cloudstor/projects/THRA_KO/project_results/DGE/WTV1WTT3_DEG_0.15LFC.csv", header = TRUE, row.names = 1)
DGE <- DGE[DGE$adj.P.Val < 0.05, ]
DGE_up <- DGE[DGE$logFC > 0, ]
DGE_down <- DGE[DGE$logFC < 0, ]

DGE_list <- list("WTVehvsT3ob_DEG" = DGE, "WTVehvsT3ob_DEGUP" = DGE_up, "WTVehvsT3ob_DEGDOWN" = DGE_down)

# Perform pathway enrichment analysis
for(i in names(DGE_list)){
  interesting <- DGE_list[[i]]
  
  BPenrich <- enrichGO(interesting[, "GeneSymbol"],
                       OrgDb = org.Mm.eg.db,
                       keyType = "SYMBOL",
                       ont = "BP",
                       pvalueCutoff = 0.05,
                       pAdjustMethod = "fdr",
                       minGSSize = 10,
                       maxGSSize = 500,
                       # readable = TRUE,
                       universe = universe)
  GOBP_GOI <- as.data.frame(BPenrich@result)
  write.csv(GOBP_GOI, paste0("project_results/t3OBs/Enriched_GOBP_", i, ".csv"))
  
  MFenrich <- enrichGO(interesting[, "GeneSymbol"],
                       OrgDb = org.Mm.eg.db,
                       keyType = "SYMBOL",
                       ont = "MF",
                       pvalueCutoff = 0.05,
                       pAdjustMethod = "fdr",
                       minGSSize = 10,
                       maxGSSize = 500,
                       # readable = TRUE,
                       universe = universe)
  GOMF_GOI <- as.data.frame(MFenrich@result)
  write.csv(GOMF_GOI, paste0("project_results/t3OBs/Enriched_GOMF_", i, ".csv"))
  
  CCenrich <- enrichGO(interesting[, "GeneSymbol"],
                       OrgDb = org.Mm.eg.db,
                       keyType = "SYMBOL",
                       ont = "CC",
                       pvalueCutoff = 0.05,
                       pAdjustMethod = "fdr",
                       minGSSize = 10,
                       maxGSSize = 500,
                       # readable = TRUE,
                       universe = universe)
  GOCC_GOI <- as.data.frame(CCenrich@result)
  write.csv(GOCC_GOI, paste0("project_results/t3OBs/Enriched_GOCC_", i, ".csv"))
  
  # Perform HPO enrichment analysis
  HPOEnrich <- enricher(gene = interesting[, "GeneSymbol"], 
                        TERM2GENE = h_t2g,
                        pvalueCutoff = 1,
                        qvalueCutoff = 1,
                        pAdjustMethod = "fdr",
                        minGSSize = 15,
                        maxGSSize = 500,
                        universe = universe)
  HPO_GOI <- as.data.frame(HPOEnrich@result)
  write.csv(HPO_GOI, paste0("project_results/t3OBs/Enriched_HPO_", i, ".csv"))
  
  TFEnrich <- enricher(gene = interesting[, "GeneSymbol"], 
                       TERM2GENE = htf_t2g,
                       pvalueCutoff = 1,
                       qvalueCutoff = 1,
                       pAdjustMethod = "fdr",
                       minGSSize = 15,
                       maxGSSize = 500,
                       universe = universe)
  TF_GOI <- as.data.frame(TFEnrich@result)
  write.csv(TF_GOI, paste0("project_results/t3OBs/Enriched_TF_", i, ".csv"))
  
  PIDEnrich <- enricher(gene = interesting[, "GeneSymbol"], 
                        TERM2GENE = hPID_t2g,
                        pvalueCutoff = 1,
                        qvalueCutoff = 1,
                        pAdjustMethod = "fdr",
                        minGSSize = 15,
                        maxGSSize = 500,
                        universe = universe)
  PID_GOI <- as.data.frame(PIDEnrich@result)
  write.csv(PID_GOI, paste0("project_results/t3OBs/Enriched_PID_", i, ".csv"))
  
  
  # Perform KEGG enrichment analysis
  KEGGenricMmig <- enrichKEGG(interesting[, "EntrezID"],
                              organism = "mmu",
                              keyType = "kegg",
                              pvalueCutoff = 0.05,
                              pAdjustMethod = "fdr",
                              universe = universe_entrez)
  write.csv(KEGGenricMmig, paste0("project_results/t3OBs/Enriched_KEGG_", i, ".csv"))
  
  # Perform REACTOME enrichment analysis
  REACTOMEenricMmig <- enrichPathway(interesting[, "EntrezID"],
                                     organism = "mouse",
                                     pvalueCutoff = 0.05,
                                     pAdjustMethod = "fdr",
                                     # readable = TRUE,
                                     universe = universe_entrez)
  write.csv(REACTOMEenricMmig, paste0("project_results/t3OBs/Enriched_REACTOME_", i, ".csv"))
}