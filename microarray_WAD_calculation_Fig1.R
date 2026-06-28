# =============================================================================
# WAD Calculation
# =============================================================================

library(readr)
library(tidyverse)

# =============================================================================
# 1. Calculate Fold Change (FC) and Weighted Average Difference (WAD)
# =============================================================================

# -----------------------------------------------------------------------------
# Load patient metadata and mRNA expression data
# -----------------------------------------------------------------------------
df <- read_csv("data/patients_and_t_log2_mRNA_data.csv")

df <- df %>%
  mutate(patient = substr(df$Title, 1, nchar(df$Title) - 1)) %>%
  filter(patient != "38" & patient != "48")

cluster <- read_csv("data/INFLAMMATORY_cluster_labels.csv")
colnames(cluster)[1] <- "patient"

df_cluster <- merge(df, cluster, by = "patient")
df_cluster$patient <- as.character(df_cluster$patient)

# -----------------------------------------------------------------------------
# Select cluster and time-point comparison
# -----------------------------------------------------------------------------
#
# Modify this section to calculate WAD for each comparison.
#
# Comparisons used in this study:
#
# Cluster 1
# - Baseline vs 10–14 days
# - Baseline vs 90 days
# - 10–14 days vs 90 days
#
# Cluster 2
# - Baseline vs 10–14 days
# - Baseline vs 90 days
# - 10–14 days vs 90 days
#
# Current example:
# Cluster 2, 10–14 days vs 90 days
#
cluster_and_time_point <- df_cluster %>%
  filter(biopsy != "pretreatment") %>%
  filter(Cluster == 2)

# -----------------------------------------------------------------------------
# Calculate fold change (FC)
# -----------------------------------------------------------------------------
cluster_and_time_point %>%
  group_by(biopsy) %>%
  summarise(across(6:22288, ~ mean(.x, na.rm = TRUE), .names = "{.col}")) -> FC

FC %>%
  as.data.frame() %>%
  column_to_rownames(var = "biopsy") %>%
  t() -> t_FC

str(t_FC)

t_FC %>%
  as.data.frame() %>%
  mutate(
    log2FC = `90 days` - `10-14 days`,
    log2ave = (`10-14 days` + `90 days`) / 2
  ) -> t_FC

# -----------------------------------------------------------------------------
# Calculate weighted average difference (WAD)
# -----------------------------------------------------------------------------
max <- max(t_FC$log2ave[is.finite(t_FC$log2ave)], na.rm = TRUE)
min <- min(t_FC$log2ave[is.finite(t_FC$log2ave)], na.rm = TRUE)

t_FC %>%
  mutate(w = (log2ave - min) / (max - min)) %>%
  mutate(WAD = log2FC * w) -> t_FC_WAD

write_csv(
  t_FC_WAD,
  "output/Cluster_2_2weeks_3months_WAD.csv"
)


# =============================================================================
# Downstream analysis
# =============================================================================
#
# After WAD calculation, representative probes were
# selected for each gene based on probe annotation.
#
# Genes were then ranked according to their WAD scores
# and subjected to Gene Set Enrichment Analysis (GSEA)
# using the clusterProfiler package.