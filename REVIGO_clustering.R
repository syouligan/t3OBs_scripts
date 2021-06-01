#!/usr/bin/Rscript

# --------------------------------------------------------------------------
#! Perform patwhay enrichment in the chondrocyte signature
# --------------------------------------------------------------------------

# Working directory
if(dir.exists("/Users/mac/cloudstor/")) {
  setwd("/Users/mac/cloudstor/imperial_projects/THcircuit/") # Uses practice data (5% of cells from each sample) if running locally
  place <- "local"
} else {
  setwd("/share/ScratchGeneral/scoyou/imperial_projects/THcircuit/")
  place <- "wolfpack"
}

library("clusterProfiler")
library("org.Mm.eg.db")
library("mclust")
library("tm")
library("ggplot2")

# Load chondrocyte gene activity data
geneActivity <- read.csv('project_results/t3OBs/Gene_annotation_w_activity_DGE.csv', row.names = 1, header = TRUE)

# Identify semantically similar clusters of significantly enriched GO terms. Removal of redundant terms and MDS scaling was performed using http://revigo.irb.hr/ (Supek, Fran, et al. "REVIGO summarizes and visualizes long lists of gene ontology terms." PloS one 6.7 (2011): e21800.)
# --------------------------------------------------------------------------
for(GO_set in c("BP", "MF", "CC")) {
  # Load REVIGO output (calculatd from GO ID and p-values from enrichGO function, C=0.9, database = Mus musculus, similarity measure = SimRel)
  revigoGO <- read.csv(paste0("project_results/t3OBs/Revigo_GO", GO_set,".csv"), header = TRUE)
  
  # Load significantly enriched GO terms
  GOsig <- read.csv(paste0("project_results/t3OBs/Enriched_GO", GO_set,"_T3ob_DEG_pairwise_list.csv"), header = TRUE)
  
  # Cluster terms based on REVIGO semantic space values and plot BIC, classification and uncertainty plots
  revigoGO <- revigoGO[! revigoGO$PlotX == " null",]
  semValsGO <- revigoGO[,c("PlotX", "PlotY")]
  semValsGO$PlotX <- as.numeric( as.character(semValsGO$PlotX) )
  semValsGO$PlotY <- as.numeric( as.character(semValsGO$PlotY) )
  
  clusterGO <- Mclust(semValsGO, modelNames = c("VII", "VVI", "VVV"), G = 1:20)
  capture.output(summary(clusterGO, parameters = TRUE), file = paste0("REVIGO_GO", GO_set,"_clustering_summary.txt"))
  
  # Make clustering QC-plots
  pdf(paste0("project_results/t3OBs/REVIGO_GO_", GO_set,"_clusters.pdf"))
  par(mfrow=c(2,2))
  plot(clusterGO, what = "BIC")
  plot(clusterGO, what = "classification")
  plot(clusterGO, what = "density")
  plot(clusterGO, what = "density", type = "persp", col = "dodgerblue3")
  dev.off()
  
  # Make dataframe with cluster classification data
  classGO <- data.frame(revigoGO)
  classGO$classification <- clusterGO$classification
  classgenesGO <- classGO[ ,c("TermID", "Name", "classification", "Value", "PlotX", "PlotY")]
  idx <- match(classgenesGO$TermID, GOsig$ID)
  classgenesGO$Genes <- GOsig$geneID [ idx ]
  
  # Make dataframe of data for plotting
  GO_class <- data.frame(classgenesGO)
  GO_class$Name <- as.character(GO_class$Name)
  GO_class$PlotX <- as.numeric(as.character(GO_class$PlotX))
  GO_class$PlotY <- as.numeric(as.character(GO_class$PlotY))
  GO_class$classification <- as.factor(as.character(GO_class$classification))
  GO_class$classification <- factor(GO_class$classification, sort(as.numeric(levels(GO_class$classification))))
  
  # Plot clusters with 50% CI and size adjusted to uncertainty
  ggplot(data = GO_class) +
    stat_ellipse(data = GO_class, aes( PlotX, PlotY, fill = classification), alpha = 0.3, geom = "polygon", type = "norm", level = 0.5, color = "magenta") +
    geom_point(aes(PlotX, PlotY, fill = classification), shape = 21, colour = "magenta", size = 6, alpha = 0.4) +
    scale_fill_brewer("Cluster", type = "qual", palette = "Set3") +
    scale_size("Certainty", range = c(1, 5), limits = c(0.5, 1)) +
    labs(x = "Semantic space X", y = "Semantic space Y") +
    theme_classic() +
    ggsave(paste0("project_results/t3OBs/GO_", GO_set,"_cluster_classification.pdf"))
  
  # Make bar plots for the top 5 GO terms in each cluster ranked by p-value (ascending)
  GO_class$Value <- abs(as.numeric(as.character(GO_class$Value)))
  all_top <- data.frame()
  
  rm(Top_genes)
  for (i in rev(levels(GO_class$classification))) {
    
    clust_top <- GO_class[GO_class$classification == i,]
    
    if(nrow(clust_top) > 2) {
      clust_top <- clust_top[with(clust_top, order(-Value)),]
      all_top <- rbind(clust_top[1:ifelse(nrow(clust_top) > 4, 5, nrow(clust_top)),], all_top)
      
      # Get top genes
      genes_tmp <- data.frame(table(unlist(strsplit(paste(clust_top$Genes, sep="", collapse="/"), "/"))))
      genes_tmp <- as.character(genes_tmp[order(genes_tmp$Freq, decreasing = TRUE), "Var1"][1:5])
      
      if(exists("Top_genes")){
        Top_genes <- rbind(Top_genes, data.frame('Label' = i, 'Genes' = paste(genes_tmp, sep="", collapse=" ")))
      } else {
        Top_genes <- data.frame('Label' = i, 'Genes' = paste(genes_tmp, sep="", collapse=" "))
      }
      
    } else {
      print("Too few")
    }
  }
  
  idx <- match(all_top$classification, Top_genes$Label)
  all_top$Genes <- Top_genes$Genes [idx]
  
  
  all_top$Label <- paste0(all_top$TermID, ": ", all_top$Name)
  all_top$Label <- factor(all_top$Label, levels = rev(all_top$Label))
  
  ggplot(data = all_top) +
    geom_col(aes(x = Label, y = Value, fill = Genes)) +
    scale_fill_brewer("Cluster", type = "qual", palette = "Set3") +
    theme_classic() +
    labs(x = "-log10 p-value") +
    scale_x_discrete(labels = function(x) str_wrap(x, width = 30)) +
    theme(axis.title.y = element_blank(), axis.text.x=element_text(angle=0, hjust=1)) +
    coord_flip() +
    ggsave(paste0("project_results/t3OBs/ClusterGO", GO_set,"_barplot.pdf"))
 }
