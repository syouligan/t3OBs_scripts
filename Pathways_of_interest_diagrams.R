#!/usr/bin/Rscript

# --------------------------------------------------------------------------
# Plot kegg pathways overlaid with fold change between T3 treated WT and Alpha-null OBs
# --------------------------------------------------------------------------

# Set adaptive working directory
if(dir.exists("/Users/mac/cloudstor/")) {
  setwd("/Users/mac/cloudstor/imperial_projects/THcircuit/") # Uses practice data (5% of cells from each sample) if running locally
  place <- "local"
} else {
  setwd("/share/ScratchGeneral/scoyou/imperial_projects/THcircuit/")
  place <- "wolfpack"
}

library('pathview')
library('dplyr')
library('KEGGREST')
library('org.Mm.eg.db')
library('reshape2')
library('ggplot2')

# Load data files and remove trailing decimal from gene ids
gene_count <- read.csv('project_results/t3OBs/Voom_transformed_CPM_active_genes.csv', row.names = 1, header = TRUE)
gene_count_long <- gene_count[,3:14]
gene_count_long <- data.frame(t(scale(t(gene_count_long))))
gene_count_long$GeneSymbol <- gene_count$GeneSymbol
gene_count_long <- melt(gene_count_long)
colnames(gene_count_long) <- c('GeneSymbol', 'Sample', 'Expression')

# Cap zscores at 1.75, -1.75
gene_count_long$Expression[gene_count_long$Expression < -1.75] <- -1.75
gene_count_long$Expression[gene_count_long$Expression > 1.75] <- 1.75

# Load gene activity annotation
geneActivity <- read.csv('project_results/t3OBs/Gene_annotation_w_activity_DGE_Venn.csv', row.names = 1, header = TRUE)
geneActivity$Ensembl <- rownames(geneActivity)

# Load differential expression info and annotate to geneActivity
WTV1_ALPHAV1 <- read.csv('project_results/t3OBs/DEGs_WTV1_ALPHAV1_0_pairwise.csv', row.names = 1, header = TRUE)
WTV1_WTT3 <- read.csv('project_results/t3OBs/DEGs_WTV1_WTT3_0_pairwise.csv', row.names = 1, header = TRUE)
WTT3_ALPHAT3 <- read.csv('project_results/t3OBs/DEGs_WTT3_ALPHAT3_0_pairwise.csv', row.names = 1, header = TRUE)
ALPHAV1_ALPHAT3 <- read.csv('project_results/t3OBs/DEGs_ALPHAV1_ALPHAT3_0_pairwise.csv', row.names = 1, header = TRUE)
DEGs <- list('WTV1_WTT3' = WTV1_WTT3, 'WTV1_ALPHAV1' = WTV1_ALPHAV1, 'WTT3_ALPHAT3' = WTT3_ALPHAT3, 'ALPHAV1_ALPHAT3' = ALPHAV1_ALPHAT3)

for(i in names(DEGs)){
  tmp <- DEGs[[i]]
  idx <- match(rownames(geneActivity), rownames(tmp))
  geneActivity[,paste0("DEG_LFC_", i)] <- tmp$logFC [idx]
}

# Load all pathways in KEGG
pathways_all <- data.frame('Description' = keggList("pathway", 'mmu'))
rownames(pathways_all) <- as.character(gsub(".*:", "", rownames(pathways_all)))
pathways_all$Description <- as.character(gsub("/", "_", pathways_all$Description))

# Keep pathways of interest
POIs <- c("mmu04310", "mmu04360", "mmu04350", "mmu04928", "mmu04918", "mmu04919", "mmu04915", 'mmu00010', 'mmu00190', 'mmu00020', 'mmu04064', 'mmu04512', 'mmu04514', 'mmu04060', "mmu00590", "mmu00564")
pathways <- pathways_all[rownames(pathways_all) %in% POIs, , drop = FALSE]

# Plot pathway diagrams and heatmaps for each POI
setwd("/Users/mac/cloudstor/imperial_projects/THcircuit/project_results/t3OBs/pathview/")

for (path in rownames(pathways)) {
  KEGGids <- path
  KEGGnames <- as.character(pathways[KEGGids, , drop = TRUE])
  
  for (comp in names(DEGs)) {
    pathLFC <- as.matrix(geneActivity[, paste0('DEG_LFC_', comp), drop = FALSE])
    rownames(pathLFC) <- geneActivity$Ensembl
    pathLFC <- pathLFC[!is.na(pathLFC), , drop = FALSE]
    # Plot on Kegg graph and export expanded list of genes from each node (this is then curated in illustrator)
    AG_classic <- pathview(gene.data = pathLFC, gene.idtype = "ensembl", gene.annotpkg = "org.Mm.eg.db", pathway.id = KEGGids, species = "mmu", out.suffix = paste0(KEGGnames, '_DEGs_', comp), limit = c(-1.75, 1.75), kegg.native = TRUE, same.layer = FALSE, node.sum = "max", both.dirs = list(gene = TRUE), low = "#023FA5", mid = "#E2E2E2", high = "#8E063B", new.signature = FALSE)
    AG_classic_data <- AG_classic$plot.data.gene
  }
  
  # Plot heatmap of all genes in pathway
  all.genes <- unique(unlist(strsplit(AG_classic_data$all.mapped, ",")))
  all.genes_symbol <- geneActivity[geneActivity$EntrezID %in% all.genes, 'GeneSymbol']

  hc <- hclust(dist(gene_count[gene_count$GeneSymbol %in% all.genes_symbol, 3:14]))$order
  hc_order <- gene_count$GeneSymbol[gene_count$GeneSymbol %in% all.genes_symbol][hc]
  
  gene_count_long %>%
    filter(GeneSymbol %in% all.genes_symbol) %>%
    mutate(Sample = factor(Sample, levels = c("WTV1_4", "WTV1_5", "WTV1_6", "ALPHAV1_4", "ALPHAV1_5", "ALPHAV1_6", "WTT3_4", "WTT3_5", "WTT3_6", "ALPHAT3_4", "ALPHAT3_5", "ALPHAT3_6"))) %>%
    mutate(GeneSymbol = factor(GeneSymbol, levels = rev(hc_order))) %>%
    ggplot(aes(x=Sample, y=GeneSymbol, fill = Expression)) + 
    geom_tile() +
    theme_classic() +
    coord_fixed() +
    theme(axis.text.y = element_text(color = 'black', face = 'italic'),
          axis.text.x = element_text(color = 'black', hjust = 1, angle = 45),
          axis.line=element_blank(),
          axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank()) +
    scale_fill_gradientn(colours = hcl.colors(7, 'RdBu', rev = TRUE), breaks = c(-1.5, 0, 1.5), limits = c(-1.75, 1.75)) +
    ggsave(paste0(KEGGids, ".", KEGGnames, ".pdf"))
}


