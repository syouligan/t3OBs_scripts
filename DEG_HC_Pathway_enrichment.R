#!/usr/bin/Rscript

# --------------------------------------------------------------------------
# Generate table of normalised counts for plotting known pathways
# --------------------------------------------------------------------------

# Set adaptive working directory
if(dir.exists("/Users/mac/cloudstor/")) {
  setwd("/Users/mac/cloudstor/imperial_projects/THcircuit/") # Uses practice data (5% of cells from each sample) if running locally
  place <- "local"
} else {
  setwd("/share/ScratchGeneral/scoyou/imperial_projects/THcircuit/")
  place <- "wolfpack"
}

library("limma")
library("edgeR")
library("ggplot2")
library("reshape")
library("dplyr")
library("tidyr")
library('viridis')
library('clusterProfiler')
library('ReactomePA')
library('org.Mm.eg.db')
library('gplots')
library('msigdbr')

# Load and prepare data
# --------------------------------------------------------------------------

# Load data files and remove trailing decimal from gene ids
gene_count <- read.csv('data/t3OBs/RSEM_count.csv', row.names = 1, header = TRUE)
rownames(gene_count) <- sub("(.*?)\\..*", "\\1", rownames(gene_count))

# Remove control samples
gene_count <- gene_count[, !grepl("Standard", colnames(gene_count))]

# Make a design table
sampleName <- colnames(gene_count)
sampleType <- factor(c(gsub( "\\_[0-9]*$", "", sampleName)))
sampleRep <- factor(c(gsub( "^.*?_","", sampleName)))
sampleTreatment <- ifelse(grepl("T3", sampleName), "T3", "V1")
sampleBatch <- ifelse(grepl("1|2|3", sampleRep), "Batch1", "Batch2")
sampleBackground <- ifelse(grepl("ALPHA", sampleName), "ALPHA", "WT")
sampleTable <- data.frame(sampleName, sampleType, sampleRep, sampleTreatment, sampleBatch, sampleBackground)

# Subset to batch 2 as there are dominant batch effects
sampleBatch <- ifelse(grepl("1|2|3", sampleRep), "Batch1", "Batch2")
gene_count <- gene_count[,sampleBatch == 'Batch2']

# Load gene annotation info
geneActivity <- read.csv('project_results/t3OBs/Gene_annotation_w_activity_DGE_Venn.csv', row.names = 1, header = TRUE)

# Perform enrichment analysis on each element of list
universe <- as.character(geneActivity[, "GeneSymbol"])
universe <- universe[!is.na(universe)]
universe_entrez <- as.character(geneActivity[, "EntrezID"])
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

# Load differential gene expression at LFC 0.15
DEGs_0.15LFC <-read.csv("summary/DEGs_LFC_0.15_pairwise.csv", header = TRUE)
DEGs_0LFC <-read.csv("summary/DEGs_LFC_0_pairwise.csv", header = TRUE)

# Make dataframe of normalised counts
# --------------------------------------------------------------------------

# Make a design table
sampleName <- colnames(gene_count)
sampleType <- factor(c(gsub( "\\_[0-9]*$", "", sampleName)))
sampleRep <- factor(c(gsub( "^.*?_","", sampleName)))
sampleTreatment <- ifelse(grepl("T3", sampleName), "T3", "V1")
sampleBatch <- ifelse(grepl("1|2|3", sampleRep), "Batch1", "Batch2")
sampleBackground <- ifelse(grepl("ALPHA", sampleName), "ALPHA", "WT")
sampleTable <- data.frame(sampleName, sampleType, sampleRep, sampleTreatment, sampleBatch, sampleBackground)

# Define experiemental design and contrasts
design <- model.matrix(~0 + sampleType + sampleRep, data = sampleTable)

contrasts <- makeContrasts( WTV1_WTT3 = sampleTypeWTT3-sampleTypeWTV1,
                            WTV1_ALPHAV1 = sampleTypeALPHAV1-sampleTypeWTV1,
                            WTV1_ALPHAT3 = sampleTypeALPHAT3-sampleTypeWTV1,
                            WTT3_ALPHAV1 = sampleTypeALPHAV1-sampleTypeWTT3,
                            WTT3_ALPHAT3 = sampleTypeALPHAT3-sampleTypeWTT3,
                            ALPHAV1_ALPHAT3 = sampleTypeALPHAT3-sampleTypeALPHAV1,
                            levels = design )

# Subset to genes expressed in any of the treatment groups
genemap <- geneActivity[geneActivity$Obs_any, ]
expressed <- gene_count[geneActivity$Obs_any, ]

# Make DGEList object and normalise using TMM
y <- DGEList(counts = expressed, group=factor(c(gsub( "\\_[0-9]*$", "", colnames(gene_count)))), genes = genemap)
y <- calcNormFactors(y, method="TMM")

# Voom transform and retrieve counts
y <- voom(y, design)
CPMvoom <- data.frame(y$E)

# Add geneSymbols
idx <- match(rownames(CPMvoom), rownames(genemap))
CPMvoom$GeneSymbol <- genemap$GeneSymbol [idx]
CPMvoom <- CPMvoom[,c("GeneSymbol", "WTV1_4", "WTV1_5", "WTV1_6", "WTT3_4", "WTT3_5", "WTT3_6", "ALPHAV1_4", "ALPHAV1_5", "ALPHAV1_6", "ALPHAT3_4", "ALPHAT3_5", "ALPHAT3_6")]

# Convert to zscores
CPMzscores <- data.frame(y$E)
CPMzscores <- data.frame(t(scale(t(CPMzscores))))
CPMzscores <- CPMzscores[,c("WTV1_4", "WTV1_5", "WTV1_6", "WTT3_4", "WTT3_5", "WTT3_6", "ALPHAV1_4", "ALPHAV1_5", "ALPHAV1_6", "ALPHAT3_4", "ALPHAT3_5", "ALPHAT3_6")]

# Make list of genes of interest
idx <- match(rownames(CPMzscores), rownames(genemap))
CPMzscores$GeneSymbol <- genemap$GeneSymbol [idx]
CPMzscores <- CPMzscores[!(CPMzscores$GeneSymbol %in% c("Nav1", "Mical3")),] # remove duplicated genes

DEGs_0LFC <- unique(DEGs_0LFC[, "GeneSymbol"])
DEGs_0LFC <- DEGs_0LFC[!(DEGs_0LFC == "" | is.na(DEGs_0LFC))]

DEGs_0.15LFC <- unique(DEGs_0.15LFC[, "GeneSymbol"])
DEGs_0.15LFC <- DEGs_0.15LFC[!(DEGs_0.15LFC == "" | is.na(DEGs_0.15LFC))]

pathways <- list("DEGs_0.15LFC" = DEGs_0.15LFC)

# Perform heirarchical clustering and pathway enrichment on each cluster
# --------------------------------------------------------------------------

for(i in names(pathways)) {
  
  zscores_GOI <- as.matrix(CPMzscores[CPMzscores$GeneSymbol %in% pathways[[i]], 1:12])

  # Perform heirarchical clustering and return gene cluster membership
  hr <- hclust(dist(zscores_GOI))
  mycl <- cutree(hr, k = 4)
  clusterCols <- viridis(4)
  myClusterSideBar <- clusterCols[mycl]

  pdf(paste0("project_results/t3OBs/DGE_HC/", i,"_total_genes_HC_expression_heatmap.pdf"))
  heatmap.2(zscores_GOI, Rowv=TRUE, Colv="none", offsetRow = 0.01, labCol = colnames(zscores_GOI), labRow = FALSE, dendrogram="row", symm=FALSE, trace = "none", density.info = "none", breaks = seq(-1.75, 1.75, 0.5), col = hcl.colors(7, 'RdBu', rev = TRUE), RowSideColors= myClusterSideBar)
  dev.off()
  
  # Extract gene ids for each cluster
  tmp <- data.frame(mycl)
  tmp$Ensembl <- rownames(tmp)
  tmp <- tmp %>%
    group_split(mycl)
  
  # Make box plot of each cluster
  for(clust in 1:length(tmp)) {
    cluster <- clust
    genes <- tmp[[clust]]$Ensembl
    
    expression <- CPMzscores[rownames(CPMzscores) %in% genes, 1:12]
    gather(expression) %>%
      mutate(key = factor(key, levels = c("WTV1_4", "WTV1_5", "WTV1_6", "WTT3_4", "WTT3_5", "WTT3_6", "ALPHAV1_4", "ALPHAV1_5", "ALPHAV1_6", "ALPHAT3_4", "ALPHAT3_5", "ALPHAT3_6"))) %>%
      ggplot(aes(x=key, y=value)) +
      geom_hline(yintercept = 0) +
      geom_boxplot() +
      scale_fill_viridis_d(option = "inferno") +
      theme_classic() +
      ggsave(paste0("project_results/t3OBs/DGE_HC/HC_", i, "_cluster_", clust, "_boxplot.pdf"), useDingbats = FALSE)
  }
  
  # Perform pathway enrichment on each cluster
  for(clust in 1:length(tmp)) {
    cluster <- clust
    genes <- tmp[[clust]]$Ensembl
    
    interesting <- geneActivity[rownames(geneActivity) %in% genes, ]

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
    write.csv(GOBP_GOI, paste0("project_results/t3OBs/DGE_HC/HC_", i, "_cluster_", clust, "_Enriched_GOBP.csv"))
    
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
    write.csv(GOMF_GOI, paste0("project_results/t3OBs/DGE_HC/HC_", i, "_cluster_", clust, "_Enriched_GOMF.csv"))
    
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
    write.csv(GOCC_GOI, paste0("project_results/t3OBs/DGE_HC/HC_", i, "_cluster_", clust, "_Enriched_GOCC.csv"))
    
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
    write.csv(HPO_GOI, paste0("project_results/t3OBs/DGE_HC/HC_", i, "_cluster_", clust, "_Enriched_HPO.csv"))
    
    TFEnrich <- enricher(gene = interesting[, "GeneSymbol"], 
                         TERM2GENE = htf_t2g,
                         pvalueCutoff = 1,
                         qvalueCutoff = 1,
                         pAdjustMethod = "fdr",
                         minGSSize = 15,
                         maxGSSize = 500,
                         universe = universe)
    TF_GOI <- as.data.frame(TFEnrich@result)
    write.csv(TF_GOI, paste0("project_results/t3OBs/DGE_HC/HC_", i, "_cluster_", clust, "_Enriched_TF.csv"))
    
    PIDEnrich <- enricher(gene = interesting[, "GeneSymbol"], 
                          TERM2GENE = hPID_t2g,
                          pvalueCutoff = 1,
                          qvalueCutoff = 1,
                          pAdjustMethod = "fdr",
                          minGSSize = 15,
                          maxGSSize = 500,
                          universe = universe)
    PID_GOI <- as.data.frame(PIDEnrich@result)
    write.csv(PID_GOI, paste0("project_results/t3OBs/DGE_HC/HC_", i, "_cluster_", clust, "_Enriched_PID.csv"))
    
    
    # Perform KEGG enrichment analysis
    KEGGenricMmig <- enrichKEGG(interesting[, "EntrezID"],
                                organism = "mmu",
                                keyType = "kegg",
                                pvalueCutoff = 0.05,
                                pAdjustMethod = "fdr",
                                universe = universe_entrez)
    write.csv(KEGGenricMmig, paste0("project_results/t3OBs/DGE_HC/HC_", i, "_cluster_", clust, "_Enriched_KEGG.csv"))
    
    # Perform REACTOME enrichment analysis
    REACTOMEenricMmig <- enrichPathway(interesting[, "EntrezID"],
                                       organism = "mouse",
                                       pvalueCutoff = 0.05,
                                       pAdjustMethod = "fdr",
                                       # readable = TRUE,
                                       universe = universe_entrez)
    write.csv(REACTOMEenricMmig, paste0("project_results/t3OBs/DGE_HC/HC_", i, "_cluster_", clust, "_Enriched_REACTOME.csv"))
  }
  
}



