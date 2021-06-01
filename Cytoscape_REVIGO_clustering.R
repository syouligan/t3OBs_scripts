#!/usr/bin/Rscript

# --------------------------------------------------------------------------
#! Make I-graph from cytoscape
# --------------------------------------------------------------------------

# Set adaptive working directory
if(dir.exists("/Users/mac/cloudstor/")) {
  setwd("/Users/mac/cloudstor/imperial_projects/THcircuit/") # Uses practice data (5% of cells from each sample) if running locally
  place <- "local"
} else {
  setwd("/share/ScratchGeneral/scoyou/imperial_projects/THcircuit/")
  place <- "wolfpack"
}

library("xml2")
library("stringr")
library("tidyverse")
library("igraph")
library('network')
library("ggnet")
library("ggplot2")
library("intergraph")
library("RColorBrewer")
library('pals')

# Load significantly enriched GO terms
BPsig <- read.csv("project_results/t3OBs/Enriched_GOBP_T3ob_DEG_pairwise_list.csv", header = TRUE)
BPsig <- BPsig[BPsig$p.adjust < 0.05,]

# Make iGraph in R from xgmml cytoscape file (from revigo). From https://gist.github.com/bkutlu/5a6f3b144d88169916f586cd2080d106
# --------------------------------------------------------------------------

# Import Revigo.xgmml
x1 <- read_xml("project_results/t3OBs/Revigo_BP.xgmml")

# Get the nodes with xpath expression
xpath_nodes <- str_replace_all("/graph/node", '/', '/d1:')
nodes <- data.frame(xml_find_all(x1, xpath_nodes, xml_ns(x1)) %>%
  xml_attrs() %>%
  do.call(rbind, .))
nodes$GOID <- nodes$id
nodes$id <- paste(nodes$GOID, nodes$label, sep = "-")
nodes <- nodes[!grepl("obsolete", nodes$id),]


# Get the edges with xpath expression
xpath_edges <- str_replace_all("/graph/edge",'/','/d1:')
edges <- xml_find_all(x1, xpath_edges, xml_ns(x1)) %>% 
  xml_attrs() %>%
  do.call(rbind, .)

# Curate into edges dataframe
edges <- data.frame(edges[,c(2:3)])
idx <- match(edges$source, nodes$GOID)
edges$from <- nodes$id [idx]
idx <- match(edges$target, nodes$GOID)
edges$to <- nodes$id [idx]
edges$connection <- paste(edges$from, edges$to, sep = ":")
edges <- edges[,c(3:5, 1:2)]

# Create iGraph object from dataframe
ig_gobp <- graph_from_data_frame(d = edges, vertices = nodes, directed = F)

# Plot top terms in each cluster
# --------------------------------------------------------------------------
BPsig$Name <- paste(BPsig$ID, BPsig$Description, sep = "-")
idx <- match(BPsig$Name, names(igraph::components(ig_gobp)$membership))
BPsig$membership <- igraph::components(ig_gobp)$membership [idx]
BPsig <- BPsig[!is.na(BPsig$membership),]

# Set singletons to the same group
idx <- match(BPsig$membership, 1:length(igraph::components(ig_gobp)$csize))
BPsig$csize <- igraph::components(ig_gobp)$csize [idx]
BPsig$membership[BPsig$csize == 1] <- 'Singletons'

# Annotate with membership information
net <- asNetwork(ig_gobp)
idx <- match(network.vertex.names(net), names(igraph::components(ig_gobp)$membership))
net %v% "membership" <- igraph::components(ig_gobp)$membership [idx]

# Set singletons to the same group
idx <- match(get.vertex.attribute(net, "membership"), 1:length(igraph::components(ig_gobp)$csize))
net %v% "csize" <- igraph::components(ig_gobp)$csize [idx]
sing <- net %v% "csize" == 1
membership <- net %v% "membership"
membership[sing] <- 'Singletons'
net %v% "membership" <- membership

# Annotate with GO information
idx <- match(network.vertex.names(net), BPsig$Name)
net %v% "Count" <- BPsig$Count [idx]
net %v% "p.adjust" <- -log10(BPsig$p.adjust) [idx]
net %v% "Description" <- BPsig$Description [idx]
net %v% "Name" <- BPsig$Name [idx]

# Label with top terms per group and top 6 genes (by count among terms)
rm(Top_terms)
rm(Top_genes)
for(i in unique(BPsig$membership)){
  # Get top terms
  tmp <- BPsig[BPsig$membership == i, ]
  num <- ifelse(nrow(tmp) < 10 | i == 'Singletons', 2, 4)
  tmp_vec <- tmp[order(tmp$p.adjust), "Name"][1:num]
  
  if(exists("Top_terms")){
    Top_terms <- c(Top_terms, tmp_vec)
  } else {
    Top_terms <- tmp_vec
  }
  
  # Get top genes
  genes_tmp <- data.frame(table(unlist(strsplit(paste(tmp$geneID, sep="", collapse="/"), "/"))))
  genes_tmp <- as.character(genes_tmp[order(genes_tmp$Freq, decreasing = TRUE), "Var1"][1:6])
  
  if(exists("Top_genes")){
    Top_genes <- rbind(Top_genes, data.frame('Cluster' = i, 'Genes' = paste(genes_tmp, sep="", collapse=" ")))
  } else {
    Top_genes <- data.frame('Cluster' = i, 'Genes' = paste(genes_tmp, sep="", collapse=" "))
  }
}

# Label clusters with unique genes
idx <- match(net %v% "membership", Top_genes$Cluster)
net %v% "membership" <- Top_genes$Genes [idx]

idx <- match(BPsig$membership, Top_genes$Cluster)
BPsig$Top_genes <- Top_genes$Genes [idx]

# Plot the GOBP network
getPalette <- cols25(n = length(unique(net %v% "membership")))
names(getPalette) <- unique(net %v% "membership")
ggnet2(net, color = "membership", size = "p.adjust", alpha = 0.5, size.cut = TRUE, color.legend = "Cluster", palette = getPalette, label = Top_terms, label.size = 2, label.color = "grey20") +
  ggsave("project_results/t3OBs/Revigo_GOBP_clustering_cytoscape.pdf")




