#!/usr/bin/Rscript

# --------------------------------------------------------------------------
#! Perform patwhay enrichment in the chondrocyte signature
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
geneActivity <- read.csv('data/annotations/Gene_annotation_mouse_ensembl_2015.csv', row.names = 1, header = TRUE)

# Perform enrichment analysis on each element of list
universe <- as.character(geneActivity[geneActivity$Biotype == "protein_coding", "GeneSymbol"])
universe <- universe[!is.na(universe)]
universe_entrez <- as.character(geneActivity[geneActivity$Biotype == "protein_coding", "EntrezID"])
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

# Make list of genes discordant in Lung sample 3
topsets <- read.csv("summary/Genes in top enriched sets.csv")
topsets <- topsets$All.256.DGE
topsets <- topsets[grepl("^ENS", topsets)]
sig <- geneActivity[geneActivity$Biotype == "protein_coding" & (rownames(geneActivity) %in% topsets), ]

DGE_list <- list("T3ob_DEG" = sig)

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
  write.csv(GOBP_GOI, paste0("project_results/t3OBs/Enriched_GOBP_", i, "_pairwise_list.csv"))
  
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
  write.csv(GOMF_GOI, paste0("project_results/t3OBs/Enriched_GOMF_", i, "_pairwise_list.csv"))
  
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
  write.csv(GOCC_GOI, paste0("project_results/t3OBs/Enriched_GOCC_", i, "_pairwise_list.csv"))
  
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
  write.csv(HPO_GOI, paste0("project_results/t3OBs/Enriched_HPO_", i, "_pairwise_list.csv"))
  
  TFEnrich <- enricher(gene = interesting[, "GeneSymbol"], 
                       TERM2GENE = htf_t2g,
                       pvalueCutoff = 1,
                       qvalueCutoff = 1,
                       pAdjustMethod = "fdr",
                       minGSSize = 15,
                       maxGSSize = 500,
                       universe = universe)
  TF_GOI <- as.data.frame(TFEnrich@result)
  write.csv(TF_GOI, paste0("project_results/t3OBs/Enriched_TF_", i, "_pairwise_list.csv"))
  
  PIDEnrich <- enricher(gene = interesting[, "GeneSymbol"], 
                        TERM2GENE = hPID_t2g,
                        pvalueCutoff = 1,
                        qvalueCutoff = 1,
                        pAdjustMethod = "fdr",
                        minGSSize = 15,
                        maxGSSize = 500,
                        universe = universe)
  PID_GOI <- as.data.frame(PIDEnrich@result)
  write.csv(PID_GOI, paste0("project_results/t3OBs/Enriched_PID_", i, "_pairwise_list.csv"))
  
  
  # Perform KEGG enrichment analysis
  KEGGenricMmig <- enrichKEGG(interesting[, "EntrezID"],
                              organism = "mmu",
                              keyType = "kegg",
                              pvalueCutoff = 1,
                              pAdjustMethod = "fdr",
                              universe = universe_entrez)
  write.csv(KEGGenricMmig, paste0("project_results/t3OBs/Enriched_KEGG_", i, "_pairwise_list.csv"))
  
  # Perform REACTOME enrichment analysis
  REACTOMEenricMmig <- enrichPathway(interesting[, "EntrezID"],
                                     organism = "mouse",
                                     pvalueCutoff = 1,
                                     pAdjustMethod = "fdr",
                                     # readable = TRUE,
                                     universe = universe_entrez)
  write.csv(REACTOMEenricMmig, paste0("project_results/t3OBs/Enriched_REACTOME_", i, "_pairwise_list.csv"))
}