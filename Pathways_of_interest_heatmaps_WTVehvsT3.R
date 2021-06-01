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

# Load differential gene expression at LFC 0.15
DGE <- read.csv("~/cloudstor/projects/THRA_KO/project_results/DGE/WTV1WTT3_DEG_0.15LFC.csv", header = TRUE, row.names = 1)
DGE <- DGE[DGE$adj.P.Val < 0.05, ]

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
CPMvoom$GOI <- genemap$GOI_0 [idx]

CPMvoom <- CPMvoom[,c( "GOI", "GeneSymbol", "WTV1_4", "WTV1_5", "WTV1_6", "WTT3_4", "WTT3_5", "WTT3_6", "ALPHAV1_4", "ALPHAV1_5", "ALPHAV1_6", "ALPHAT3_4", "ALPHAT3_5", "ALPHAT3_6")]
write.csv(CPMvoom, "project_results/t3OBs/Voom_transformed_CPM_active_genes.csv", row.names = TRUE)

# Make heatmaps of genes/pathways of interest
# --------------------------------------------------------------------------

# Convert to zscores
CPMzscores <- data.frame(y$E)
CPMzscores <- CPMzscores[, grepl("WT", colnames(CPMzscores))]
CPMzscores <- data.frame(t(scale(t(CPMzscores))))

# Transform to long format
idx <- match(rownames(CPMzscores), rownames(genemap))
CPMzscores$GeneSymbol <- genemap$GeneSymbol [idx]
CPMzscores <- CPMzscores[!(CPMzscores$GeneSymbol %in% c("Nav1", "Mical3")),] # remove duplicated genes

CPMzscores_long <- melt(CPMzscores)
colnames(CPMzscores_long) <- c('GeneSymbol', 'Sample', 'Expression')

# List pathways of interest
RANK_RANKL_OPG <- c("Tnfsf11", "Tnfrsf11a", "Tnfrsf11b")
WNT <- c("Wnt2b", "Wnt4", "Wnt5a", "Wnt5b", "Wnt6", "Wnt9a", "Wnt10a", "Wnt10b", "Wnt11", "Wnt16", "Sfrp1", "Sfrp2", "Sfrp4", "Notum", "Wif1", "Dkk2", "Dkk3", "Serpinf1", "Ctnnb1")
TGFB_BMP_SMAD <- c("Tgfb1", "Tgfb2", "Tgfb3", "Bmp1", "Bmp2", "Bmpr1a", "Bmpr1b", "Bmpr2", "Bmp2k", "Bmp3", "Bmp4", "Bmp5", "Bmp6", "Gdf6", "Smad1", "Smad2", "Smad3", "Smad4", "Smad5", "Smad6", "Smad7", "Smad9")
EPHRIN_EPH <- c("Efna1", "Efna2", "Efna3", "Efna4", "Efna5", "Efnb1", "Efnb2", "Efnb3", "Epha1", "Epha2", "Epha3", "Epha4", "Epha5", "Epha7", "Ephb2", "Ephb3", "Ephb4", "Ephb6", "Ephx1", "Ephx2")
ESTROGEN <- c("Esr1", "Src", "Ncoa1", "Ncoa2", "Ncoa3", "Fos", "Jun", "Sp1", "Spns2", "Tgfa")
SEMAPHORIN_PLEXIN <- c("Sema3a", "Sema3b", "Sema3c", "Sema3d", "Sema3e", "Sema3f", "Sema4b", "Sema4c", "Sema4d", "Sema4g", "Sema5a", "Sema6a", "Sema6c", "Sema6d", "Sema7a", "Plxna1", "Plxna2", "Plxna3", "Plxna4", "Plxnb1", "Plxnb2", "Plxnc1", "Plxnd1")
GLYCOLYSIS <- c("Eno2", "Eno3", "Eno1", "Pgk2", "Pgk1", "Pgam1", "Pgam2", "Gapdhs", "Gapdh", "Tpi1", "Gpi", "Aldoa", "Aldob", "Aldoc", "Pfkl", "Pfkm", "Pfkp", "Adpgk", "Gck", "Gckr", "Nup98", "Seh1l", "Nup107", "Nup43", "Nup160", "Nup37", "Nup85", "Nup133", "Sec13", "Tpr", "Nup153", "Nup88", "Nup188", "Nup214", "Ndc1", "Nup210", "Pom121", "Pom121c", "Nup35", "Nup93", "Nup155", "Nup205", "Seh1l", "Rae1", "Nup98", "Nup98", "Nup62", "Nup58", "Nup58", "Nup54", "Nup50", "Aaas", "Ranbp2", "Nup42", "Gnpda1", "Gnpda2", "Hk1", "Hk2", "Hk3", "Pkm", "Pklr", "Pklr", "Pkm", "Bpgm", "Pgm2l1", "Pgp", "Ppp2r1a", "Ppp2r1b", "Ppp2ca", "Ppp2cb", "Ppp2r5d", "Pfkfb1", "Pfkfb3", "Pfkfb4", "Pfkfb2", "Prkaca", "Prkacb", "Prkacg")
PYRUVATE <- c("Dld", "Glo1", "Hagh", "Mpc1", "Mpc2", "Pdhb", "Pdha2", "Pdha1", "Dlat", "Pdhx", "Ldhal6b", "Bsg", "Slc16a3", "Slc16a8", "Slc16a1", "Ldhal6a", "Ldhb", "Ldha", "Ldhc", "Me1", "Gstz1", "Pdk2", "Pdk4", "Pdk3", "Pdk1", "Rxra", "Ppard", "Pdp1", "Pdpr", "Pdp2", "Vdac1")
TCA <- c("Sdha", "Sdhb", "Sdhd", "Sdhc", "Mdh2", "Me2", "Fh", "Ogdh", "Dlst", "Dld", "Idh2", "Sucla2", "Suclg1", "Me3", "Idh3b", "Idh3a", "Idh3g", "Aco2", "Fahd1", "Nnt", "Cs", "Suclg2")
RESP_TRANSPORT <- c("Ndufs5", "Ndufs6", "Ndufv1", "Ndufv3", "Ndufa12", "Ndufs4", "Ndufs1", "Ndufv2", "mt-Nd4", "mt-Nd5", "mt-Nd1", "Ndufa9", "Ndufs3", "Ndufs8", "Ndufs7", "Ndufs2", "Ndufb8", "Ndufa13", "Ndufb10", "Ndufb7", "Ndufa10", "Ndufc1", "Ndufb2", "Ndufb5", "Ndufb9", "Ndufa11", "Ndufa2", "Ndufa7", "Ndufa8", "Ndufa1", "Ndufa3", "Ndufc2", "Ndufb1", "Ndufab1", "Ndufa5", "Ndufb3", "Ndufa6", "Ndufb4", "Ndufb11", "Coa1", "mt-Nd2", "Ndufb6", "mt-Nd3", "Tmem186", "mt-Nd6", "Coq10a", "Coq10b", "Etfdh", "Etfb", "Etfa", "Ndufaf2", "Ndufaf7", "Ndufaf3", "Ndufaf4", "Timmdc1", "Tmem126b", "Ecsit", "Ndufaf1", "Acad9", "Ndufaf5", "Ndufaf6", "Nubpl", "Uqcrfs1", "mt-Cyb", "Uqcr10", "Uqcrq", "Uqcrc2", "Cyc1", "Uqcrh", "Uqcrb", "Uqcrc1", "Uqcr11", "Cycs", "Cox6c", "Cox6b1", "Ndufa4", "mt-Co2", "Cox6a1", "Cox7c", "Cox5b", "Cox4i1", "Cox7a2l", "Cox5a", "Cox7b", "mt-Co1", "mt-Co3", "Cox8a", "Sco2", "Taco1", "Sco1", "Surf1", "Lrpprc", "Cox19", "Cox18", "Cox14", "Cox16", "Cox20", "Cox11", "Sdha", "Sdhb", "Sdhd", "Sdhc", "Trap1")
MITO_ATP <- c("Atp5po", "Atp5mc3", "mt-Atp6", "Atp5me", "Atp5mf", "mt-Atp8", "Atp5mc1", "Atp5mg", "Atp5pf", "Atp5pd", "Atp5pb", "Atp5mc2", "Atp5f1b", "Atp5f1c", "Atp5f1d", "Atp5f1e", "Atp5f1a", "Dmac2l")
MITO_UNCOUPLING <- c("Hr", "Vdr", "Ucp1", "Ucp2", "Ucp3", "Pm20d1", "Slc25a27", "Slc25a14")
GOIs <- geneActivity[geneActivity$GOI_0, "GeneSymbol"]
DGE <- DGE[, "GeneSymbol"]
DGE <- DGE[!(is.na(DGE) | DGE == "")]

GLYCOLYSIS_LEE <- c("Slc2a1", "Slc2a3", "Slc2a6", "Slc2a8", "Slc2a10", "Slc2a13", "Hk1", "Hk2", "Gpi1", "Pfkfb2", "Pfkfb3", "Pfkfb4", "Pfkl", 'Pfkm', 'Pfkp', 'Aldoa', 'Aldoc', 'Tpi1', 'Gapdh', 'Pgk1', 'Pgam1', 'Pgam2', 'Eno1', 'Eno2', 'Eno3', 'Pkm')
TCA_LEE <- c("Cs", "Aco1", "Aco2", "Idh2", "Idh3a", "Idh3b", "Idh3g", "Ogdh", "Ogdhl", "Dlst", "Dld", "Sucla2", "Suclg1", 'Suclg2', 'Sdha', 'Sdhb', 'Sdhc', 'Sdhd', 'Fh1', 'Mdh2')
FATTYACID_LEE <- c("Acat1", "Acat2", "Hadh", "Hadha", "Acaa1a", "Acaa2", "Echs1", "Acox3", "Acadm", "Acads", "Acadsb", "Acadvl", "Gcdh", 'Acsl1', 'Acsl3', 'Acsl4', 'Cpt1a', 'Cpt1b', 'Cpt2', 'Eci1', 'Eci2', 'Adh1', 'Adh5', 'Adh7', 'Aldh2', 'Aldh3a2', 'Aldh7a1', 'Aldh9a1')
MITOGENESIS_LEE <- c("Ppargc1b", "Nrf1", "Nfo2l2", "Tfam", "Dnm1l", "Mfn1", "Mfn2", "Opa1")
ETCC1_LEE <- c("mt-Nd1", "mt-Nd2", "mt-Nd4", "mt-Nd5", "mt-Nd6", "Ndufs1", "Ndufs2", "Ndufs3", "Ndufs4", "Ndufs5", "Ndufs6", "Ndufs7", "Ndufs8", 'Ndufv1', 'Ndufv2', 'Ndufv3', 'Ndufa1', 'Ndufa2', 'Ndufa3', 'Ndufa4', 'Ndufa4l2', 'Ndufa5', 'Ndufa6', 'Ndufa7', 'Ndufa8', 'Ndufa9',  'Ndufa10', 'Ndufa11', 'Ndufa12', 'Ndufa13', 'Ndufb2', 'Ndufb3', 'Ndufb4', 'Ndufb5', 'Ndufb6', 'Ndufb7', 'Ndufb8', 'Ndufb9', 'Ndufb10', 'Ndufb11',  'Ndufc1',  'Ndufc2')
ETCC2_LEE <- c("Sdha", "Sdhb", "Sdhc", "Sdhd")
ETCC3_LEE <- c("Cyc1", "mt-Cycb", "Uqcrb", "Uqcrc1", "Uqcrc2", "Uqcrfs1", "Uqcrh", "Uqcrq", "Uqcr10", "Uqcr11")
ETCC4_LEE <- c("mt-Co1", "Cox4l1", "Cox4l2", "Cox5a", "Coxb", "Cox6a1", "Cox6a2", "Cox6b1", "Cox6b2", "Cox6c", "Cox7a1", "Cox7a2", "Cox7a2l", 'Cox8a', 'Cox10', 'Cox11', 'Cox15', 'Cox17')
ETCC5_LEE <- c("Atp5a1", "Atp5b", "Atp5c1", "Atp5d", "Atp5e", "Atp5f1", "Atp5g1", "Atp5g3", "Atp5h", "Atp5j", "Atp5j2", "Atp5k", "Atp5l", 'Atp5o', 'Atp5s')

pathways <- list('RANK_RANKL_OPG' = RANK_RANKL_OPG, 'WNT' = WNT, 'TGFB_BMP_SMAD' = TGFB_BMP_SMAD, 'EPHRIN_EPH' = EPHRIN_EPH, 'ESTROGEN' = ESTROGEN, 'SEMAPHORIN_PLEXIN' = SEMAPHORIN_PLEXIN, "GLYCOLYSIS" = GLYCOLYSIS, "PYRUVATE" = PYRUVATE, "TCA" = TCA, "RESP_TRANSPORT" = RESP_TRANSPORT, "MITO_ATP" = MITO_ATP, "MITO_UNCOUPLING" = MITO_UNCOUPLING, 'GOIs' = GOIs, "DGE" = DGE)
pathways <- list('GLYCOLYSIS_LEE' = GLYCOLYSIS_LEE, 'TCA_LEE' = TCA_LEE, 'FATTYACID_LEE' = FATTYACID_LEE, 'MITOGENESIS_LEE' = MITOGENESIS_LEE, 'ETCC1_LEE' = ETCC1_LEE, 'ETCC2_LEE' = ETCC2_LEE, "ETCC3_LEE" = ETCC3_LEE, "ETCC4_LEE" = ETCC4_LEE, "ETCC5_LEE" = ETCC5_LEE, "DGE" = DGE)

# Cap zscores at 1.75, -1.75
CPMzscores_long$Expression[CPMzscores_long$Expression < -1.5] <- -1.5
CPMzscores_long$Expression[CPMzscores_long$Expression > 1.5] <- 1.5

# Plot pathways of interest
for(i in names(pathways)) {
  
  hc <- hclust(dist(CPMzscores[CPMzscores$GeneSymbol %in% pathways[[i]], 1:6]))$order
  hc_order <- CPMzscores$GeneSymbol[CPMzscores$GeneSymbol %in% pathways[[i]]][hc]
  
  CPMzscores_long %>%
    filter(GeneSymbol %in% pathways[[i]]) %>%
    mutate(Sample = factor(Sample, levels = c("WTV1_4", "WTV1_5", "WTV1_6", "WTT3_4", "WTT3_5", "WTT3_6"))) %>%
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
    scale_fill_gradientn(colours = hcl.colors(7, 'RdBu', rev = TRUE), breaks = c(-1.5, 0, 1.5), limits = c(-1.5, 1.5)) +
    ggsave(paste0("project_results/t3OBs/Pathway_", i, "_expression_WTVehvsT3.pdf"))
}


# Plot heatmap with dendrogram of DEGs
zscores_GOI <- as.matrix(CPMzscores[CPMzscores$GeneSymbol %in% DGE, 1:6])
hr <- hclust(dist(zscores_GOI))
mycl <- cutree(hr, k = 2)
clusterCols <- viridis(2)
myClusterSideBar <- clusterCols[mycl]

col <- hcl.colors(14, 'RdBu', rev = TRUE)
col[7:8] <- 'white'

pdf(paste0("project_results/t3OBs/WTVehvsT3_DGE_HC_expression_heatmap.pdf"))
heatmap.2(zscores_GOI, Rowv=TRUE, Colv="none", offsetRow = 0.01, labCol = colnames(zscores_GOI), labRow = FALSE, dendrogram="row", symm=FALSE, trace = "none", density.info = "none", breaks = seq(-1.75, 1.75, 0.25), col = col, RowSideColors= myClusterSideBar)
dev.off()

