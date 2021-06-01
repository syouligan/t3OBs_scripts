#!/usr/bin/Rscript

# --------------------------------------------------------------------------
# Find genes differentially expressed in Astrocytes_Morte
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

# Set
set = "Astrocytes_Morte"

# Load data files and remove trailing decimal from gene ids
gene_count <- read.csv(paste0('data/public_data/', set, '/Total_counts.csv'), row.names = 1, header = TRUE)

# Load gene annotation info
geneActivity <- read.csv('project_results/t3OBs/Gene_annotation_w_activity.csv', row.names = 1, header = TRUE)

# Make a design table
sampleName <- colnames(gene_count)
sampleType <- factor(c(gsub( "\\_[0-9]*$", "", sampleName)))
sampleRep <- factor(c(gsub( "^.*?_","", sampleName)))
sampleTreatment <- ifelse(grepl("T3", sampleName), "T3", "V1")
sampleExtra <- ifelse(grepl("CHX", sampleName), "CHX", "Control")
sampleTime <- ifelse(grepl("6Hr", sampleName), "6Hr", "24Hr")
# sampleTable <- data.frame(sampleName, sampleType, sampleRep, sampleTreatment, sampleBatch, sampleBackground)
sampleTable <- data.frame(sampleName, sampleType, sampleRep, sampleTreatment, sampleTime, sampleExtra)

# Establish experimental design
# --------------------------------------------------------------------------

# Define experiemental design and contrasts
design <- model.matrix(~0 + sampleType, data = sampleTable)

contrasts <- makeContrasts(X6hrV1_X6hrT3 = sampleTypeX6hrT3-sampleTypeX6hrV1,
                           X24hrV1_X24hrT3 = sampleTypeX24hrT3-sampleTypeX24hrV1,
                           X6hrCHXV1_X6hrCHXT3 = sampleTypeX6hrCHXT3-sampleTypeX6hrCHXV1,
                            levels = design )

# Filter to genes with at least 3 counts in any groups
binarised <- data.frame(gene_count > 1)
groups <- unique(sampleTable$sampleType)
for(grp in groups){
  cols <- grepl(grp, colnames(binarised))
  binarised[, paste0(grp, "_activity")] <- rowSums(binarised[, cols])
  
  idx <- match(rownames(geneActivity), rownames(binarised))
  geneActivity[,paste0(grp, "_activity")] <- binarised[,paste0(grp, "_activity")] [idx]
}

# Perform differential gene expression with conditional (comparison dependent) background
# --------------------------------------------------------------------------

# Set LFC for testing with treat
LFC = 0

for(cont in colnames(contrasts)) {
  # Subset data to genes active in either group
  group1 <- strsplit2(cont, "_")[1]
  group2 <- strsplit2(cont, "_")[2]
  
  # Subset to active genes
  genemap_tmp <- geneActivity[which(geneActivity[,paste0(group1, "_activity")] == sum(grepl(group1, sampleTable$sampleType)) | geneActivity[,paste0(group2, "_activity")] == sum(grepl(group2, sampleTable$sampleType))), ]
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
  write.csv(res, paste0("data/public_data/", set, "/DEGs_", cont,"_", LFC, "_", set, ".csv"))
  
  # Write gene annotation csv with gene activity information
  idx <- match(rownames(geneActivity), rownames(res))
  geneActivity[,paste0(set, "_", cont, "_LFC")] <- res$logFC [idx]
  geneActivity[,paste0(set, "_", cont, "_t")] <- res$t [idx]
  geneActivity[,paste0(set, "_", cont, "_adj.P.Val")] <- res$adj.P.Val [idx]
}


