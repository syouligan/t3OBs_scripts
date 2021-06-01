#!/usr/bin/Rscript

# --------------------------------------------------------------------------
# Find genes differentially expressed WT and ALPHA KO osteoblasts treated with T3
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

# Load and prepare data
# --------------------------------------------------------------------------

# Load data files and remove trailing decimal from gene ids
gene_count <- read.csv('data/t3OBs/RSEM_count.csv', row.names = 1, header = TRUE)
rownames(gene_count) <- sub("(.*?)\\..*", "\\1", rownames(gene_count))

# Remove control samples
gene_count <- gene_count[, !grepl("Standard", colnames(gene_count))]

# Subset to batch 2 as there are dominant batch effects
sampleBatch <- ifelse(grepl("_1|_2|_3", colnames(gene_count)), "Batch1", "Batch2")
gene_count <- gene_count[,sampleBatch == 'Batch2']

# Load gene annotation info
geneActivity <- read.csv('project_results/t3OBs/Gene_annotation_w_activity.csv', row.names = 1, header = TRUE)

# Make a design table
sampleName <- colnames(gene_count)
sampleType <- factor(c(gsub( "\\_[0-9]*$", "", sampleName)))
sampleRep <- factor(c(gsub( "^.*?_","", sampleName)))
sampleTreatment <- ifelse(grepl("T3", sampleName), "T3", "V1")
sampleBatch <- ifelse(grepl("1|2|3", sampleRep), "Batch1", "Batch2")
sampleBackground <- ifelse(grepl("ALPHA", sampleName), "ALPHA", "WT")
sampleTable <- data.frame(sampleName, sampleType, sampleRep, sampleTreatment, sampleBatch, sampleBackground)

# Establish experimental design
# --------------------------------------------------------------------------

# Define experiemental design and contrasts
design <- model.matrix(~0 + sampleType + sampleRep, data = sampleTable)

contrasts <- makeContrasts( WTV1_WTT3 = sampleTypeWTT3-sampleTypeWTV1,
                            WTV1_ALPHAV1 = sampleTypeALPHAV1-sampleTypeWTV1,
                            WTT3_ALPHAT3 = sampleTypeALPHAT3-sampleTypeWTT3,
                            ALPHAV1_ALPHAT3 = sampleTypeALPHAT3-sampleTypeALPHAV1,
                            levels = design )

# Perform differential gene expression with conditional (comparison dependent) background
# --------------------------------------------------------------------------

# Set LFC for testing with treat
for(LFC in c(0, 0.15, 0.3, 0.45)) {
  for(cont in colnames(contrasts)) {
    # Subset data to genes active in either group
    group1 <- strsplit2(cont, "_")[1]
    group2 <- strsplit2(cont, "_")[2]
    
    # Subset to active genes
    genemap_tmp <- geneActivity[geneActivity[,paste0(group1, "_activity")] == 3 | geneActivity[,paste0(group2, "_activity")] == 3, ]
    expression_tmp <- gene_count[rownames(gene_count) %in% rownames(genemap_tmp),]
    
    # Make DGEList object and normalise using TMM
    y <- DGEList(counts = expression_tmp, group=sampleType, genes = genemap_tmp)
    y <- calcNormFactors(y, method="TMM")
    
    # Voom transform, fit linear models and caluclate differential expression
    y <- voom(y, design)
    fit <- lmFit(y, design)
    fit <- contrasts.fit(fit, contrasts)
    fit <- treat(fit, lfc = LFC)
    DEG1 <- decideTests(fit)
    print(summary(DEG1))
    
    # Save results and add to geneActivity dataframe
    res <- data.frame(topTreat(fit, n = Inf, p.value = Inf, coef = cont))
    write.csv(res, paste0("project_results/t3OBs/DEGs_", cont,"_", LFC,".csv"))
    
    # Write gene annotation csv with gene activity information
    idx <- match(rownames(geneActivity), rownames(res))
    geneActivity[,paste0(cont,"_", LFC, "_LFC")] <- res$logFC [idx]
    geneActivity[,paste0(cont,"_", LFC, "_t")] <- res$t [idx]
    geneActivity[,paste0(cont,"_", LFC, "_adj.P.Val")] <- res$adj.P.Val [idx]
    geneActivity[,paste0(cont,"_", LFC, "_significant")] <- geneActivity[,paste0(cont,"_", LFC, "_adj.P.Val")] < 0.05
  }
}

write.csv(geneActivity, "project_results/t3OBs/Gene_annotation_w_activity_DGE.csv")
