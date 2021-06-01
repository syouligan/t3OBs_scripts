#!/usr/bin/Rscript

# --------------------------------------------------------------------------
#! Plot significant REACTOME terms
# --------------------------------------------------------------------------

# Working directory
if(dir.exists("/Users/mac/cloudstor/")) {
  setwd("/Users/mac/cloudstor/imperial_projects/THcircuit/") # Uses practice data (5% of cells from each sample) if running locally
  place <- "local"
} else {
  setwd("/share/ScratchGeneral/scoyou/imperial_projects/THcircuit/")
  place <- "wolfpack"
}

library(ggplot2)
library(dplyr)
library(tidyr)
library(forcats)

# Load and prepare data
# --------------------------------------------------------------------------

# Load chondrocyte gene activity data
geneActivity <- read.csv('project_results/t3OBs/Gene_annotation_w_activity_DGE.csv', row.names = 1, header = TRUE)

# Load REACTOME enrichment data
REACTOME <- read.csv("project_results/t3OBs/Enriched_REACTOME_T3ob_DEG_pairwise_list.csv", header = TRUE)

# Extract number of genes in each group
REACTOME <- separate(data = REACTOME, col = BgRatio, into = c("Group", "Total"), sep = "\\/")
REACTOME$Group <- as.numeric(REACTOME$Group)
REACTOME$Percent <- as.integer(REACTOME$Count / REACTOME$Group *100)

# Subset to non-redundant terms, removing terms with 100% gene overlap keeping top term (ranked by p-value)
REACTOME <- data.frame(REACTOME %>%
  group_by(geneID) %>%
  top_n(-1, p.adjust))

# Plot non-redundant terms
REACTOME %>%
  mutate(Description = fct_reorder(Description, Percent)) %>%
  ggplot() +
  geom_vline(xintercept = seq(0, 20, 10), linetype = "dashed", color="black") +
  geom_point(aes(x = Percent, y = Description, fill = -log10(p.adjust), size = Group), shape = 21, color = "white") +
  scale_fill_gradientn(aes(fill = -log10(p.adjust)), colors = c("orange", "purple"), limits = c(1.30103, max(-log10(REACTOME$p.adjust))), values = NULL, space = "Lab", na.value = "grey80", guide = "colourbar", aesthetics = "fill") +
  scale_size_continuous("Group size", range = c(1,9)) +
  scale_y_discrete(labels = function(y) str_wrap(y, width = 30)) +
  scale_x_continuous(breaks = seq(0, 25, 5), labels = seq(0, 25, 5), position = "top", expand = c(0,0.5)) +
  labs(x = "% of REACTOME genes", x = NA) +
  theme_classic() +
  ggsave("project_results/t3OBs/REACTOME_signature_percent_dotplot.pdf", useDingbats = FALSE)

REACTOME %>%
  mutate(Description = fct_reorder(Description, Percent)) %>%
  top_n(50, Percent) %>%
  ggplot() +
  geom_vline(xintercept = seq(0, 20, 10), linetype = "dashed", color="black") +
  geom_point(aes(x = Percent, y = Description, fill = -log10(p.adjust), size = Group), shape = 21, color = "white") +
  scale_fill_gradientn(aes(fill = -log10(p.adjust)), colors = c("orange", "purple"), limits = c(1.30103, max(-log10(REACTOME$p.adjust))), values = NULL, space = "Lab", na.value = "grey80", guide = "colourbar", aesthetics = "fill") +
  scale_size_continuous("Group size", range = c(1,9)) +
  scale_x_continuous(breaks = seq(0, 25, 5), labels = seq(0, 25, 5), position = "top", expand = c(0,0.5)) +
  labs(x = "% of REACTOME genes", x = NA) +
  theme_classic() +
  ggsave("project_results/t3OBs/REACTOME_signature_percent_dotplot_no_wrapping.pdf", useDingbats = FALSE)

