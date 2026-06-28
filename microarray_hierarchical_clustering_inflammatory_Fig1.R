# =============================================================================
# Hierarchical Clustering Using HALLMARK_INFLAMMATORY_RESPONSE Gene Set
# Methods:
#   - Fold changes of HALLMARK_INFLAMMATORY_RESPONSE genes between
#     baseline → 2 weeks and 2 weeks → 3 months used as clustering features
#   - Hierarchical clustering with Euclidean distance and Ward's method
# =============================================================================


# -----------------------------------------------------------------------------
# 1. Load libraries
# -----------------------------------------------------------------------------
library(tidyverse)


# -----------------------------------------------------------------------------
# 2. Load data
# -----------------------------------------------------------------------------
patients_log2_mRNA_data <- read_csv("data/patients_and_t_log2_mRNA_data.csv")
HALLMARK_INFLAMMATORY   <- read_csv("data/HALLMARK_INFLAMMATORY.csv")
Genename                <- read_csv("data/gene_name.csv")
probes                  <- read_csv("data/after_probe_selection_gene_expression_data.csv")


# -----------------------------------------------------------------------------
# 3. Filter samples and compute log2 fold changes between consecutive timepoints
#    Timepoints per patient:
#      T1: baseline, T2: 2 weeks, T3: 3 months
#    Features used for clustering:
#      T2 - T1 (baseline → 2 weeks)
#      T3 - T2 (2 weeks  → 3 months)
# -----------------------------------------------------------------------------
expr_filtered <- patients_log2_mRNA_data %>%
  filter(!Title %in% c("38A", "38B", "48A", "48B")) %>%
  mutate(patient = substr(Title, 1, nchar(Title) - 1))

df_diff <- expr_filtered %>%
  group_by(patient) %>%
  mutate(across(5:22287, ~ .x - lag(.x), .names = "{.col}"))

# Remove lag-NA rows (T1 rows) and retain only expression columns
df_diff2 <- df_diff %>%
  dplyr::select(-22288) %>%
  na.omit() %>%
  dplyr::select(5:22288)


# -----------------------------------------------------------------------------
# 4. Select best probes and rename columns to Gene Symbols
# -----------------------------------------------------------------------------
Genename2 <- Genename %>%
  dplyr::select("ID", "Gene Symbol", "ENTREZ_GENE_ID") %>%
  mutate(
    `Gene Symbol`  = sapply(strsplit(as.character(`Gene Symbol`),  " /// "), `[`, 1),
    ENTREZ_GENE_ID = sapply(strsplit(as.character(ENTREZ_GENE_ID), " /// "), `[`, 1)
  )

# Retain only best probes identified in preprocessing
df_diff2_cleaned <- df_diff2 %>%
  dplyr::select(any_of(probes$selected_probe))

Genename2_cleaned <- Genename2 %>%
  filter(ID %in% probes$selected_probe) %>%
  arrange(`Gene Symbol`)

colnames(df_diff2_cleaned)[2:ncol(df_diff2_cleaned)] <- Genename2_cleaned$`Gene Symbol`


# -----------------------------------------------------------------------------
# 5. Extract HALLMARK_INFLAMMATORY_RESPONSE genes
# -----------------------------------------------------------------------------
hallmark_genes <- HALLMARK_INFLAMMATORY$Gene_SYMBOL[
  HALLMARK_INFLAMMATORY$Gene_SYMBOL %in% colnames(df_diff2_cleaned)
]

df_diff3 <- df_diff2_cleaned %>%
  dplyr::select(patient, all_of(hallmark_genes))


# -----------------------------------------------------------------------------
# 6. Reshape to wide format
#    Each row: one patient
#    Columns: gene_T1 (baseline→2w fold change), gene_T2 (2w→3m fold change)
# -----------------------------------------------------------------------------
feature_cols <- colnames(df_diff3)[-1]  # exclude patient column

df_wide <- df_diff3 %>%
  group_by(patient) %>%
  mutate(timepoint = row_number()) %>%
  ungroup() %>%
  pivot_wider(
    names_from  = timepoint,
    values_from = all_of(feature_cols),
    names_glue  = "{.value}_T{timepoint}"
  )

write_csv(df_wide, "output/df_for_clustering.csv")


# -----------------------------------------------------------------------------
# 7. Hierarchical clustering
#    Distance: Euclidean
#    Linkage:  Ward's method (ward.D2)
#    k = 2 clusters
# -----------------------------------------------------------------------------
df_wide_mat <- df_wide %>%
  column_to_rownames(var = "patient")

set.seed(123)
dist_matrix <- dist(df_wide_mat, method = "euclidean")
hc           <- hclust(dist_matrix, method = "ward.D2")

# Dendrogram
plot(hc,
     main = "Hierarchical Clustering: HALLMARK_INFLAMMATORY_RESPONSE",
     sub  = "", xlab = "", cex = 0.8)
rect.hclust(hc, k = 2, border = c("red", "blue"))

# Assign cluster labels (k = 2)
cluster_labels <- cutree(hc, k = 2)

result <- data.frame(
  Sample  = rownames(df_wide_mat),
  Cluster = cluster_labels
)
print(result)

write_csv(result, "output/INFLAMMATORY_cluster_labels.csv")
