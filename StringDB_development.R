#!/usr/bin/Rscript

# --------------------------------------------------------------------------
# StringDB analysis of top hits
# --------------------------------------------------------------------------

# Set adaptive working directory
if(dir.exists("/Users/mac/cloudstor/")) {
  setwd("/Users/mac/cloudstor/imperial_projects/THcircuit/") # Uses practice data (5% of cells from each sample) if running locally
  place <- "local"
} else {
  setwd("/share/ScratchGeneral/scoyou/imperial_projects/THcircuit/")
  place <- "wolfpack"
}

library('dplyr')
library('org.Mm.eg.db')
library('reshape2')
library('ggplot2')
library('STRINGdb')
library('igraph')

# Load gene activity annotation
geneActivity <- read.csv('project_results/t3OBs/Gene_annotation_w_activity_DGE_Venn.csv', row.names = 1, header = TRUE)
geneActivity$Ensembl <- rownames(geneActivity)

# Look here for igraph example https://www.bioconductor.org/packages/release/bioc/vignettes/netSmooth/inst/doc/buildingPPIsFromStringDB.html
# Look here for vignette https://www.bioconductor.org/packages/release/bioc/vignettes/STRINGdb/inst/doc/STRINGdb.pdf

string_db <- STRINGdb$new(version="11", species = 10090, score_threshold=400, input_directory="")
mouse_graph <- string_db$get_graph()
edge.scores <- E(mouse_graph)$combined_score
ninetyth.percentile <- quantile(edge.scores, 0.9)
thresh <- data.frame(name='90th percentile', val=ninetyth.percentile)
mouse_graph <- subgraph.edges(mouse_graph, E(mouse_graph)[combined_score > ninetyth.percentile])

hr <- string_db$mp( "hr" )
neighbours <- string_db$get_neighbors(hr)
string_db$plot_network( neighbours )
clustersList = string_db$get_clusters(neighbours)
par(mfrow=c(2,2))
for(i in seq(1:4)){
  string_db$plot_network(clustersList[[i]])
}
