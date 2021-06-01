#!/usr/bin/Rscript

# --------------------------------------------------------------------------
#! Analyse known ob/oc crosstalk pathways
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
library("limma")
library("edgeR")
library("ggplot2")
library("reshape")
library("pathview")
library('KEGGREST')

# Makes HPO geneset for enrichment testing
h_df <- msigdbr(species = "Mus musculus", category = "C5", subcategory = "GO:BP")
h_t2g <- h_df %>% dplyr::select(gs_name, gene_symbol) %>% as.data.frame()

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
geneActivity <- read.csv('project_results/t3OBs/Gene_annotation_w_activity.csv', row.names = 1, header = TRUE)

# Establish experimental design
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

# Perform differential gene expression with conditional (comparison dependent) background
# --------------------------------------------------------------------------

# Subset to active genes
genemap_tmp <- geneActivity[geneActivity$Obs_any, ]
expression_tmp <- gene_count[rownames(gene_count) %in% rownames(genemap_tmp),]

# Make a design table
sampleName <- colnames(expression_tmp)
sampleType <- factor(c(gsub( "\\_[0-9]*$", "", sampleName)))
sampleRep <- factor(c(gsub( "^.*?_","", sampleName)))
sampleTable_tmp <- data.frame(sampleName, sampleType, sampleRep)
design_tmp <- model.matrix(~0 + sampleType + sampleRep, data = sampleTable_tmp)
    
# Make DGEList object and normalise using TMM
y <- DGEList(counts = expression_tmp, group=factor(c(gsub( "\\_[0-9]*$", "", colnames(expression_tmp)))), genes = genemap_tmp)
y <- calcNormFactors(y, method="TMM")
    
# Voom transform, fit linear models and caluclate differential expression
y <- voom(y, design_tmp)
fit <- lmFit(y, design_tmp)
fit <- contrasts.fit(fit, contrasts)
fit <- treat(fit, lfc = LFC)
DEG1 <- decideTests(fit)
print(summary(DEG1))
    
# Save results 
res <- data.frame(topTreat(fit, n = Inf, p.value = Inf))
write.csv(res, "project_results/t3OBs/DEGs_all_contrasts_0LFC.csv")

# Plot key kegg pathways on KEGG diagram
# --------------------------------------------------------------------------

pathways <- data.frame('Description' = keggList("pathway", 'mmu'))
rownames(pathways) <- as.character(gsub(".*:", "", rownames(pathways)))
pathways$Description <- as.character(gsub("/", "_", pathways$Description))

# Subset to genes that are active in every bone type
pathLFC <- as.matrix(res[, "adj.P.Val", drop = FALSE])

# Plot on Kegg graph and export expanded list of genes from each node (this is then curated in illustrator)
KEGGids <- "mmu04928"
KEGGnames <- "Parathyroid hormone synthesis, secretion and action"
AG_classic <- pathview(gene.data = pathLFC, gene.idtype = "ensembl", gene.annotpkg = "org.Mm.eg.db", pathway.id = KEGGids, species = "mmu", out.suffix = KEGGnames, kegg.native = TRUE, same.layer = FALSE, node.sum = "max", limit = c(0, 1), bins = list(gene = 14), both.dirs = list(gene = TRUE), new.signature = FALSE)
AG_classic_data <- AG_classic$plot.data.gene
AG_classic_data$all.mapped <- gsub(",", "/", AG_classic_data$all.mapped)
write.csv(AG_classic_data, "PTH_classic_data.csv")

# Export gene mapping in the KEGG pathway diagram
AG_all <- pathview(gene.data = pathLFC, gene.idtype = "ensembl", gene.annotpkg = "org.Mm.eg.db", pathway.id = KEGGids, species = "mmu", out.suffix = KEGGnames, kegg.native = FALSE, map.symbol = TRUE, same.layer = FALSE, split.group = TRUE, expand.node = TRUE, node.sum = "max", limit = c(0, 1), bins = list(gene = ), both.dirs = list(gene = TRUE), low = "#e6f598", mid = "#e6f598", high = "#3288bd", new.signature = FALSE)
write.csv(AG_all$plot.data.gene, "PTH_all_expanded_data.csv")
AG_all_data <- AG_all$plot.data.gene
AG_all_data <- AG_all_data[! is.na(AG_all_data$ge1), ]
AG_all_data$labels <- as.factor(as.character(AG_all_data$labels, levels = AG_all_data$labels))

ggplot(data = AG_all_data) +
  geom_col(aes(x = labels, y = ge1, fill = ge1)) +
  scale_fill_gradientn(colours = c("#fff7b2", "#fff7b2", "#aad4b2", "#7dc7c1", "#4db2d2", "#2c8cbe", "#0b589e"), values = c(0, 0.34, 0.42, 0.5, 0.58, 0.66 ,1), breaks = seq(-3, 7, 1), limits = c(-3, 7)) +
  scale_y_continuous(breaks = seq(-5, 7, 1), limits = c(-5, 7)) +
  theme_classic() +
  labs(y = "-log10 p-value") +
  coord_flip() +
  ggsave("PTH_pathway_barplot.pdf")



# Arrange data for plotting
# --------------------------------------------------------------------------

# Save normalised gene expression formatted for graphpad
norm_expression <- data.frame(y$E)
norm_expression$Genes <- rownames(norm_expression)
idx <- match(norm_expression$Genes, rownames(geneActivity))
norm_expression$GeneSymbol <- geneActivity$GeneSymbol [idx]

# Create gene sets of interest
SMAD_BMP <- geneActivity[grepl('^Bmp[0-9]|^Smad', geneActivity$GeneSymbol) & geneActivity$Obs_any, ]
EPHRIN <- geneActivity[grepl('^Efn', geneActivity$GeneSymbol) & geneActivity$Obs_any, ]
SEMA <- geneActivity[grepl('^Sema', geneActivity$GeneSymbol) & geneActivity$Obs_any, ]
RANKL <- geneActivity[geneActivity$GeneSymbol %in% c('Tnfrsf11a', 'Tnfsf11', 'Tnfrsf11b') & geneActivity$Obs_any, ]
OTHERS <- geneActivity[geneActivity$GeneSymbol %in% c('Cthrc1', 'Csf1', 'Ablim1') & geneActivity$Obs_any, ]

GOI_lists <- list('SMAD_BMP' = SMAD_BMP, 'EPHRIN' = EPHRIN, 'SEMA' = SEMA, 'RANKL' = RANKL, 'OTHERS' = OTHERS)

