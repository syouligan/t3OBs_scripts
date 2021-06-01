#!/usr/bin/Rscript

# --------------------------------------------------------------------------
# Find correlation and overlap of DEGs with V1vsT3 response in osteoblasts and other tissues
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
library("ggpubr")
library("ggrepel")
library("dplyr")

# Osteoblasts DEGs
OB_WTV1_WTT3 <- read.csv("project_results/t3OBs/DEGs_WTV1_WTT3_0_pairwise.csv")
OB_WTV1_WTT3_DEGs <- OB_WTV1_WTT3[OB_WTV1_WTT3$adj.P.Val < 0.05, ]
OB_WTV1_WTT3$OB_DEGs <- OB_WTV1_WTT3$adj.P.Val < 0.05

# Calculate correlation with Muscle_Nicolaisen LFC
# --------------------------------------------------------------------------

# Muscle Nicolaisen
MUSCLE_WTV1_WTT3 <- read.csv("data/public_data/Muscle_Nicolaisen/DEGs_WTV1_WTT3_0_Muscle_Nicolaisen.csv")
MUSCLE_WTV1_WTT3_DEGs <- MUSCLE_WTV1_WTT3[MUSCLE_WTV1_WTT3$adj.P.Val < 0.05, ]

idx <- match(OB_WTV1_WTT3$X, MUSCLE_WTV1_WTT3$X)
OB_WTV1_WTT3$MUSCLE_logFC <- MUSCLE_WTV1_WTT3$logFC [idx]
OB_WTV1_WTT3$MUSCLE_adj.P.Val <- MUSCLE_WTV1_WTT3$adj.P.Val [idx]
OB_WTV1_WTT3$MUSCLE_DEGs <- OB_WTV1_WTT3$MUSCLE_adj.P.Val < 0.05

OB_WTV1_WTT3 <- OB_WTV1_WTT3 %>% 
  mutate(OB_MUSCLE_DEGs =
           case_when(
             MUSCLE_DEGs & OB_DEGs ~ "Both",
             MUSCLE_DEGs & OB_DEGs == FALSE ~ "MUSCLE_only",
             MUSCLE_DEGs == FALSE & OB_DEGs ~ "Osteoblasts_only",
             MUSCLE_DEGs == FALSE & OB_DEGs == FALSE ~ "Neither")
  )

# Calculate correlation and plot relationship
OB_WTV1_WTT3 %>%
  filter(OB_DEGs | MUSCLE_DEGs) %>%
  ggplot(aes(x=logFC, y=MUSCLE_logFC)) +
  geom_point(alpha=0.3) +
  geom_label_repel(data = OB_WTV1_WTT3[1:20,], aes(x=logFC, y=MUSCLE_logFC, label = GeneSymbol, fill = OB_MUSCLE_DEGs)) +
  stat_cor(method = "spearman", cor.coef.name = "rho", alternative = "two.sided") +
  geom_smooth(method = 'lm', se = TRUE) +
  theme_minimal() +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  theme(aspect.ratio = 1) +
  ggsave("data/public_data/Muscle_Nicolaisen/Correlation_with_Obs_V1vsT3.png")

# Calculate correlation with Striatum_Richard LFC
# --------------------------------------------------------------------------

# Striatum Richard
STRIATUM_WTV1_WTT3 <- read.csv("data/public_data/Striatum_Richard/DEGs_PTUOnly_PTUT3_0_Striatum_Richard.csv")
STRIATUM_WTV1_WTT3_DEGs <- STRIATUM_WTV1_WTT3[STRIATUM_WTV1_WTT3$adj.P.Val < 0.05, ]

idx <- match(OB_WTV1_WTT3$GeneSymbol, STRIATUM_WTV1_WTT3$X)
OB_WTV1_WTT3$STRIATUM_logFC <- STRIATUM_WTV1_WTT3$logFC [idx]
OB_WTV1_WTT3$STRIATUM_adj.P.Val <- STRIATUM_WTV1_WTT3$adj.P.Val [idx]
OB_WTV1_WTT3$STRIATUM_DEGs <- OB_WTV1_WTT3$STRIATUM_adj.P.Val < 0.05

OB_WTV1_WTT3 <- OB_WTV1_WTT3 %>% 
  mutate(OB_STRIATUM_DEGs =
           case_when(
             STRIATUM_DEGs & OB_DEGs ~ "Both",
             STRIATUM_DEGs & OB_DEGs == FALSE ~ "STRIATUM_only",
             STRIATUM_DEGs == FALSE & OB_DEGs ~ "Osteoblasts_only",
             STRIATUM_DEGs == FALSE & OB_DEGs == FALSE ~ "Neither")
  )

# Calculate correlation and plot relationship
OB_WTV1_WTT3 %>%
  filter(OB_DEGs | STRIATUM_DEGs) %>%
  ggplot(aes(x=logFC, y=STRIATUM_logFC)) +
  geom_point(alpha=0.3) +
  geom_label_repel(data = OB_WTV1_WTT3[1:20,], aes(x=logFC, y=STRIATUM_logFC, label = GeneSymbol, fill = OB_STRIATUM_DEGs)) +
  stat_cor(method = "spearman", cor.coef.name = "rho", alternative = "two.sided") +
  geom_smooth(method = 'lm', se = TRUE) +
  theme_minimal() +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  theme(aspect.ratio = 1) +
  ggsave("data/public_data/Striatum_Richard/Correlation_with_Obs_V1vsT3.png")

# Calculate correlation with Liver_Finan LFC
# --------------------------------------------------------------------------

# Liver_Finan
LIVER_WTV1_WTT3 <- read.csv("data/public_data/Liver_Finan/DEGs_V1_T3_0_Liver_Finan.csv")
LIVER_WTV1_WTT3_DEGs <- LIVER_WTV1_WTT3[LIVER_WTV1_WTT3$adj.P.Val < 0.05, ]

idx <- match(OB_WTV1_WTT3$GeneSymbol, LIVER_WTV1_WTT3$X)
OB_WTV1_WTT3$LIVER_logFC <- LIVER_WTV1_WTT3$logFC [idx]
OB_WTV1_WTT3$LIVER_adj.P.Val <- LIVER_WTV1_WTT3$adj.P.Val [idx]
OB_WTV1_WTT3$LIVER_DEGs <- OB_WTV1_WTT3$LIVER_adj.P.Val < 0.05

OB_WTV1_WTT3 <- OB_WTV1_WTT3 %>% 
  mutate(OB_LIVER_DEGs =
           case_when(
             LIVER_DEGs & OB_DEGs ~ "Both",
             LIVER_DEGs & OB_DEGs == FALSE ~ "LIVER_only",
             LIVER_DEGs == FALSE & OB_DEGs ~ "Osteoblasts_only",
             LIVER_DEGs == FALSE & OB_DEGs == FALSE ~ "Neither")
  )

# Calculate correlation and plot relationship
OB_WTV1_WTT3 %>%
  filter(OB_DEGs | LIVER_DEGs) %>%
  ggplot(aes(x=logFC, y=LIVER_logFC)) +
  geom_point(alpha=0.3) +
  geom_label_repel(data = OB_WTV1_WTT3[1:20,], aes(x=logFC, y=LIVER_logFC, label = GeneSymbol, fill = OB_LIVER_DEGs)) +
  stat_cor(method = "spearman", cor.coef.name = "rho", alternative = "two.sided") +
  geom_smooth(method = 'lm', se = TRUE) +
  theme_minimal() +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  theme(aspect.ratio = 1) +
  ggsave("data/public_data/Liver_Finan/Correlation_with_Obs_V1vsT3.png")

# Calculate correlation with Cerebrocortical_Gil_Ibanez LFC
# --------------------------------------------------------------------------

# Cerebrocortical_Gil_Ibanez
CEREBROCORTICAL_WTV1_WTT3 <- read.csv("data/public_data/Cerebrocortical_Gil_Ibanez/DEGs_V1_T3_0_Cerebrocortical_Gil_Ibanez.csv")
CEREBROCORTICAL_WTV1_WTT3_DEGs <- CEREBROCORTICAL_WTV1_WTT3[CEREBROCORTICAL_WTV1_WTT3$adj.P.Val < 0.05, ]

idx <- match(OB_WTV1_WTT3$X, CEREBROCORTICAL_WTV1_WTT3$X)
OB_WTV1_WTT3$CEREBROCORTICAL_logFC <- CEREBROCORTICAL_WTV1_WTT3$logFC [idx]
OB_WTV1_WTT3$CEREBROCORTICAL_adj.P.Val <- CEREBROCORTICAL_WTV1_WTT3$adj.P.Val [idx]
OB_WTV1_WTT3$CEREBROCORTICAL_DEGs <- OB_WTV1_WTT3$CEREBROCORTICAL_adj.P.Val < 0.05

OB_WTV1_WTT3 <- OB_WTV1_WTT3 %>% 
  mutate(OB_CEREBROCORTICAL_DEGs =
           case_when(
             CEREBROCORTICAL_DEGs & OB_DEGs ~ "Both",
             CEREBROCORTICAL_DEGs & OB_DEGs == FALSE ~ "CEREBROCORTICAL_only",
             CEREBROCORTICAL_DEGs == FALSE & OB_DEGs ~ "Osteoblasts_only",
             CEREBROCORTICAL_DEGs == FALSE & OB_DEGs == FALSE ~ "Neither")
  )

# Calculate correlation and plot relationship
OB_WTV1_WTT3 %>%
  filter(OB_DEGs | CEREBROCORTICAL_DEGs) %>%
  ggplot(aes(x=logFC, y=CEREBROCORTICAL_logFC)) +
  geom_point(alpha=0.3) +
  geom_label_repel(data = OB_WTV1_WTT3[1:20,], aes(x=logFC, y=CEREBROCORTICAL_logFC, label = GeneSymbol, fill = OB_CEREBROCORTICAL_DEGs)) +
  stat_cor(method = "spearman", cor.coef.name = "rho", alternative = "two.sided") +
  geom_smooth(method = 'lm', se = TRUE) +
  theme_minimal() +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  theme(aspect.ratio = 1) +
  ggsave("data/public_data/Cerebrocortical_Gil_Ibanez/Correlation_with_Obs_V1vsT3.png")

# Calculate correlation with Astrocytes_Morte LFC
# --------------------------------------------------------------------------

# Astrocytes_Morte
ASTROCYTES_WTV1_WTT3 <- read.csv("data/public_data/Astrocytes_Morte/DEGs_X24hrV1_X24hrT3_0_Astrocytes_Morte.csv")
ASTROCYTES_WTV1_WTT3_DEGs <- ASTROCYTES_WTV1_WTT3[ASTROCYTES_WTV1_WTT3$adj.P.Val < 0.05, ]

idx <- match(OB_WTV1_WTT3$X, ASTROCYTES_WTV1_WTT3$X)
OB_WTV1_WTT3$ASTROCYTES_logFC <- ASTROCYTES_WTV1_WTT3$logFC [idx]
OB_WTV1_WTT3$ASTROCYTES_adj.P.Val <- ASTROCYTES_WTV1_WTT3$adj.P.Val [idx]
OB_WTV1_WTT3$ASTROCYTES_DEGs <- OB_WTV1_WTT3$ASTROCYTES_adj.P.Val < 0.05

OB_WTV1_WTT3 <- OB_WTV1_WTT3 %>% 
  mutate(OB_ASTROCYTES_DEGs =
           case_when(
             ASTROCYTES_DEGs & OB_DEGs ~ "Both",
             ASTROCYTES_DEGs & OB_DEGs == FALSE ~ "ASTROCYTES_only",
             ASTROCYTES_DEGs == FALSE & OB_DEGs ~ "Osteoblasts_only",
             ASTROCYTES_DEGs == FALSE & OB_DEGs == FALSE ~ "Neither")
  )

# Calculate correlation and plot relationship
OB_WTV1_WTT3 %>%
  filter(OB_DEGs | ASTROCYTES_DEGs) %>%
  ggplot(aes(x=logFC, y=ASTROCYTES_logFC)) +
  geom_point(alpha=0.3) +
  geom_label_repel(data = OB_WTV1_WTT3[1:20,], aes(x=logFC, y=ASTROCYTES_logFC, label = GeneSymbol, fill = OB_ASTROCYTES_DEGs)) +
  stat_cor(method = "spearman", cor.coef.name = "rho", alternative = "two.sided") +
  geom_smooth(method = 'lm', se = TRUE) +
  theme_minimal() +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  theme(aspect.ratio = 1) +
  ggsave("data/public_data/Astrocytes_Morte/Correlation_with_Obs_V1vsT3.png")

# Calculate correlation with Intestine_Godart LFC
# --------------------------------------------------------------------------

# Intestine_Godart
INTESTINE_WTV1_WTT3 <- read.csv("data/public_data/Intestine_Godart/DEGs_V1_T3_0_Intestine_Godart.csv")
INTESTINE_WTV1_WTT3_DEGs <- INTESTINE_WTV1_WTT3[INTESTINE_WTV1_WTT3$adj.P.Val < 0.05, ]

idx <- match(OB_WTV1_WTT3$GeneSymbol, INTESTINE_WTV1_WTT3$X)
OB_WTV1_WTT3$INTESTINE_logFC <- INTESTINE_WTV1_WTT3$logFC [idx]
OB_WTV1_WTT3$INTESTINE_adj.P.Val <- INTESTINE_WTV1_WTT3$adj.P.Val [idx]
OB_WTV1_WTT3$INTESTINE_DEGs <- OB_WTV1_WTT3$INTESTINE_adj.P.Val < 0.05

OB_WTV1_WTT3 <- OB_WTV1_WTT3 %>% 
  mutate(OB_INTESTINE_DEGs =
           case_when(
             INTESTINE_DEGs & OB_DEGs ~ "Both",
             INTESTINE_DEGs & OB_DEGs == FALSE ~ "INTESTINE_only",
             INTESTINE_DEGs == FALSE & OB_DEGs ~ "Osteoblasts_only",
             INTESTINE_DEGs == FALSE & OB_DEGs == FALSE ~ "Neither")
  )

# Calculate correlation and plot relationship
OB_WTV1_WTT3 %>%
  filter(OB_DEGs | INTESTINE_DEGs) %>%
  ggplot(aes(x=logFC, y=INTESTINE_logFC)) +
  geom_point(alpha=0.3) +
  geom_label_repel(data = OB_WTV1_WTT3[1:20,], aes(x=logFC, y=INTESTINE_logFC, label = GeneSymbol, fill = OB_INTESTINE_DEGs)) +
  stat_cor(method = "spearman", cor.coef.name = "rho", alternative = "two.sided") +
  geom_smooth(method = 'lm', se = TRUE) +
  theme_minimal() +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  theme(aspect.ratio = 1) +
  ggsave("data/public_data/Intestine_Godart/Correlation_with_Obs_V1vsT3.png")


