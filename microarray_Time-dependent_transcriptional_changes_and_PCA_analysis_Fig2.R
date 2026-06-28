# =============================================================================
# Time-dependent transcriptional changes and PCA analysis
# =============================================================================
#
# Workflow:
# 1. Identify time-dependent genes within each cluster
# 2. Perform PCA using significant genes
# 3. Evaluate association between PC1 and treatment duration
# 4. Perform GSEA using PCA loadings


library(readr)
library(tidyverse)
library(clusterProfiler)
library(DOSE)
library(org.Hs.eg.db)
library(enrichplot)
library(ggridges)
library(ggpubr)

##################### Load data ##############################

processed <-
  read_csv("data/after_probe_selection_gene_expression_data_processed.csv")

patients_and_t_log2_mRNA_data <-
  read_csv(
    "data/patients_and_t_log2_mRNA_data.csv"
  ) %>%
  dplyr::select(
    ID,
    Title
  )

patients_cluster_information <-
  read_csv(
    "data/patients_Cluster_information.csv"
  )

##################### Sample information ##############################

exprs_mat <-
  processed[,2:156] %>%
  column_to_rownames("gene") %>%
  as.matrix()

info <-
  inner_join(
    patients_and_t_log2_mRNA_data,
    patients_cluster_information,
    by = "Title"
  )

##################### Cluster selection ##############################

# Example:
# Cluster 1 or Cluster 2

cluster_id <- 1

info_cluster <-
  info %>%
  filter(
    Cluster == cluster_id
  )

exprs_mat_cluster <-
  exprs_mat[
    ,
    info_cluster$ID
  ]

exprs_mat_cluster[
  !is.finite(exprs_mat_cluster)
] <- 0

##################### Time dependent genes ##############################

info_cluster$time_factor <-
  factor(
    info_cluster$biopsy
  )

pvals <-
  apply(
    exprs_mat_cluster,
    1,
    function(x){
      
      df_tmp <-
        data.frame(
          expr = x,
          time = info_cluster$time_factor
        )
      
      fit <-
        aov(
          expr ~ time,
          data = df_tmp
        )
      
      summary(fit)[[1]][["Pr(>F)"]][1]
    }
  )

pvals_adj <-
  p.adjust(
    pvals,
    method = "BH"
  )

sig_genes <-
  names(
    pvals_adj[
      pvals_adj < 0.05
    ]
  )

exprs_mat_cluster <-
  exprs_mat_cluster %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  filter(
    gene %in% sig_genes
  ) %>%
  column_to_rownames("gene") %>%
  as.matrix()

##################### PCA ##############################

pca_res <-
  prcomp(
    t(exprs_mat_cluster),
    scale. = TRUE
  )

pca_df <-
  as.data.frame(
    pca_res$x
  ) %>%
  rownames_to_column(
    "ID"
  ) %>%
  inner_join(
    info,
    by = "ID"
  )

pca_df$biopsy <-
  factor(
    pca_df$biopsy,
    levels = c(
      "pretreatment",
      "10-14 days",
      "90 days"
    )
  )

##################### PCA visualization ##############################

ggplot(
  pca_df,
  aes(
    PC1,
    PC2,
    color = biopsy
  )
) +
  geom_point(size = 3) +
  geom_path(
    aes(group = patient),
    linewidth = 0.8,
    alpha = 0.5,
    color = "gray50"
  ) +
  theme_bw()

##################### PC1 and treatment time points ##############################

cor_pca_df <-
  pca_df %>%
  mutate(
    time = case_when(
      biopsy == "pretreatment" ~ 0,
      biopsy == "10-14 days" ~ 1,
      biopsy == "90 days" ~ 2
    )
  ) %>%
  dplyr::select(
    PC1,
    time
  )

ggplot(
  cor_pca_df,
  aes(
    PC1,
    time
  )
) +
  geom_point(size = 3) +
  geom_smooth(method = "lm") +
  stat_cor(method = "spearman") +
  theme_bw()

cor(
  cor_pca_df$PC1,
  cor_pca_df$time
)

##################### GSEA using PC1 loadings ##############################

loadings <-
  pca_res$rotation

PC1_vec <-
  loadings[, "PC1"]

symbol_to_entrez <-
  mapIds(
    org.Hs.eg.db,
    keys = rownames(loadings),
    column = "ENTREZID",
    keytype = "SYMBOL",
    multiVals = "first"
  )

PC1_df <-
  data.frame(
    ENTREZID = symbol_to_entrez,
    loading = PC1_vec
  )

PC1_df <-
  PC1_df %>%
  filter(
    !is.na(ENTREZID)
  ) %>%
  group_by(
    ENTREZID
  ) %>%
  slice_max(
    abs(loading),
    n = 1,
    with_ties = FALSE
  ) %>%
  ungroup()

gene_list <-
  PC1_df$loading

names(gene_list) <-
  PC1_df$ENTREZID

gene_list <-
  sort(
    gene_list,
    decreasing = TRUE
  )

hallmark_gene_sets <-
  read.gmt(
    "data/h.all.v2024.1.Hs.entrez.gmt"
  )

gsea_results <-
  GSEA(
    geneList = gene_list,
    TERM2GENE = hallmark_gene_sets,
    minGSSize = 20,
    verbose = FALSE,
    pvalueCutoff = 0.1
  )

##################### Visualization ##############################

positive_gsea <-
  gsea_results

positive_gsea@result <-
  gsea_results@result[
    gsea_results@result$NES > 0,
  ]

dotplot(
  positive_gsea,
  x = "NES",
  showCategory = 10,
  color = "p.adjust",
  size = "setSize"
)

positive_gsea <-
  pairwise_termsim(
    positive_gsea
  )

emapplot(
  positive_gsea,
  showCategory = 20
)

##################### Notes ##############################

# Example shown for Cluster 2.
# The same workflow can be applied to Cluster 1.
#
# Hallmark gene sets were obtained from MSigDB.
# Gene identifiers were converted from HGNC symbols
# to ENTREZ IDs using org.Hs.eg.db.
# GSEA was performed using clusterProfiler.