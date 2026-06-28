# =============================================================================
# Fibroblast Subclustering and CAF Annotation
# Methods:
#   - Tumor-derived fibroblast clusters defined as CAFs
#   - CAFs classified into iCAF, myCAF, apCAF based on marker genes
#     (inflammation-related, ECM-related, ap-related)
#   - ECM-related CAFs further subclustered into ECM related 1 and ECM related 2
# =============================================================================


# -----------------------------------------------------------------------------
# 1. Load libraries
# -----------------------------------------------------------------------------
library(tidyverse)
library(Seurat)
library(harmony)
library(viridis)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)


# -----------------------------------------------------------------------------
# 2. Load data and subset fibroblasts
# -----------------------------------------------------------------------------
all_cells <- readRDS("data/all_cells.rds")

Fibro_obj <- subset(all_cells, subset = new_label == "Fibroblasts")


# -----------------------------------------------------------------------------
# 3. Preprocessing and dimensionality reduction (Fibroblasts)
# -----------------------------------------------------------------------------
Fibro_obj <- NormalizeData(Fibro_obj)
Fibro_obj <- FindVariableFeatures(Fibro_obj, selection.method = "vst", nfeatures = 2000)
Fibro_obj <- ScaleData(Fibro_obj)
Fibro_obj <- RunPCA(Fibro_obj, features = VariableFeatures(Fibro_obj), npcs = 40)

ElbowPlot(Fibro_obj, ndims = 40)  # confirm dims = 30

Fibro_obj <- RunHarmony(Fibro_obj, group.by.vars = "orig.ident", dims.use = 1:30)

Fibro_obj_res <- Fibro_obj %>%
  FindNeighbors(reduction = "harmony", dims = 1:30) %>%
  FindClusters(resolution = 0.1) %>%
  RunUMAP(reduction = "harmony", dims = 1:30)

DimPlot(Fibro_obj_res, reduction = "umap", raster = FALSE)
DimPlot(Fibro_obj_res, reduction = "umap", group.by = "TissueType", raster = FALSE)

# Marker gene expression to confirm fibroblast identity
FeaturePlot(Fibro_obj_res, features = c("PDGFRA", "PDGFRB"), raster = FALSE)  # fibroblast common
FeaturePlot(Fibro_obj_res, features = c("FAP", "PDPN"), raster = FALSE)       # CAF markers
FeaturePlot(Fibro_obj_res, features = c("THY1", "ACTA2", "COL1A1", "COL1A2",
                                        "IL6", "CXCL12", "IL1A", "LIF",
                                        "CD74", "SLPI", "HLA-DRA"))           # iCAF, myCAF, apCAF

DotPlot(Fibro_obj_res,
        features = c("THY1", "ACTA2", "COL1A1", "COL1A2",
                     "IL6",  "CXCL12", "IL1A",  "LIF",
                     "CD74", "SLPI",   "HLA-DRA"),
        group.by = "seurat_clusters") +
  scale_color_viridis(option = "turbo") +
  RotatedAxis()


# -----------------------------------------------------------------------------
# 4. Define CAFs: subset clusters predominantly from tumor tissue
# -----------------------------------------------------------------------------
CAFs_obj <- subset(Fibro_obj_res, subset = seurat_clusters %in% c("1", "3"))

DimPlot(CAFs_obj, reduction = "umap", raster = FALSE)


# -----------------------------------------------------------------------------
# 5. Re-cluster CAFs and classify into iCAF / myCAF / apCAF
# -----------------------------------------------------------------------------
CAFs_obj <- FindVariableFeatures(CAFs_obj, selection.method = "vst", nfeatures = 5000)

CAFs_obj <- CAFs_obj %>%
  ScaleData(features = VariableFeatures(CAFs_obj)) %>%
  RunPCA(features = VariableFeatures(CAFs_obj), npcs = 40) %>%
  RunHarmony(group.by.vars = "orig.ident", dims.use = 1:30) %>%
  FindNeighbors(reduction = "harmony", dims = 1:30) %>%
  FindClusters(resolution = 0.3) %>%
  RunUMAP(reduction = "harmony", dims = 1:30)

DimPlot(CAFs_obj, reduction = "umap", raster = FALSE, label = TRUE)

# Marker gene expression by CAF subtype
FeaturePlot(CAFs_obj,
            features = c("IL6", "CXCL12", "IL1A", "LIF"),
            cols = viridis(10, option = "turbo"))   # iCAF

FeaturePlot(CAFs_obj,
            features = c("TAGLN", "ACTA2", "COL1A1", "COL1A2"),
            cols = viridis(10, option = "turbo"))   # myCAF

FeaturePlot(CAFs_obj,
            features = c("CD74", "SLPI", "HLA-DRA"),
            cols = viridis(10, option = "turbo"))   # apCAF

FeaturePlot(CAFs_obj,
            features = c("LOX", "LOXL1", "LOXL2"),
            cols = viridis(10, option = "turbo"))   # ECM remodeling

DotPlot(CAFs_obj,
        features = c("TAGLN", "ACTA2", "COL1A1", "COL1A2",
                     "IL6",   "CXCL12", "CCL2",  "LIF",
                     "CD74",  "SLPI",   "HLA-DRA",
                     "LOX",   "LOXL1",  "LOXL2")) +
  scale_color_viridis(option = "turbo") +
  RotatedAxis()

# Module scores for CAF subtypes
apCAF_genes <- c("CD74", "SLPI", "HLA-DRA")
myCAF_genes <- c("TAGLN", "ACTA2", "COL1A1", "COL1A2")
iCAF_genes  <- c("IL6", "CXCL12", "CCL2", "LIF")

CAFs_obj <- AddModuleScore(CAFs_obj,
                           features = list(apCAF_genes, myCAF_genes, iCAF_genes),
                           name     = c("apCAFscore", "myCAFscore", "iCAFscore"))

FeaturePlot(CAFs_obj,
            features = c("apCAFscore1", "myCAFscore2", "iCAFscore3"),
            cols = viridis(10, option = "turbo"))

# -----------------------------------------------------------------------------
# 6. Remove low-quality / ambiguous cluster and re-cluster for final annotation
# -----------------------------------------------------------------------------
# Cluster 5 is excluded based on marker gene inspection above
CAFs_obj_clean <- subset(CAFs_obj, subset = seurat_clusters != "5")

CAFs_obj_clean <- FindVariableFeatures(CAFs_obj_clean,
                                       selection.method = "vst",
                                       nfeatures = 5000)

CAFs_obj_clean <- CAFs_obj_clean %>%
  ScaleData(features = VariableFeatures(CAFs_obj_clean)) %>%
  RunPCA(features = VariableFeatures(CAFs_obj_clean), npcs = 40) %>%
  RunHarmony(group.by.vars = "orig.ident", dims.use = 1:30) %>%
  FindNeighbors(reduction = "harmony", dims = 1:30) %>%
  FindClusters(resolution = 0.5) %>%
  RunUMAP(reduction = "harmony", dims = 1:30)

DimPlot(CAFs_obj_clean, reduction = "umap", raster = FALSE, label = TRUE)

# Confirm subtype identity with marker genes and violin plots
DotPlot(CAFs_obj_clean,
        features = c("TAGLN", "ACTA2", "COL1A1", "COL1A2",
                     "IL6",   "CXCL12", "CCL2",  "LIF",
                     "CD74",  "SLPI",   "HLA-DRA",
                     "LOX",   "LOXL1",  "LOXL2")) +
  scale_color_viridis(option = "turbo") +
  RotatedAxis()

FeaturePlot(CAFs_obj_clean,
            features = c("IL6", "CXCL12", "CCL2", "LIF"),
            cols = viridis(10, option = "turbo"))   # iCAF

FeaturePlot(CAFs_obj_clean,
            features = c("TAGLN", "ACTA2", "COL1A1", "COL1A2"),
            cols = viridis(10, option = "turbo"))   # myCAF

FeaturePlot(CAFs_obj_clean,
            features = c("CD74", "SLPI", "HLA-DRA"),
            cols = viridis(10, option = "turbo"))   # apCAF

FeaturePlot(CAFs_obj_clean,
            features = c("LOX", "LOXL1", "LOXL2"),
            cols = viridis(10, option = "turbo"))   # ECM remodeling


# Annotate clusters
Idents(CAFs_obj_clean) <- "seurat_clusters"
CAFs_obj_clean <- RenameIdents(CAFs_obj_clean,
                               "0" = "ECM related",
                               "1" = "ECM related",
                               "2" = "ECM related",
                               "3" = "ECM related",
                               "4" = "Inflammation related",
                               "5" = "Inflammation related",
                               "6" = "Inflammation related",
                               "7" = "Ap related",
                               "8" = "ECM related",
                               "9" = "ECM related",
                               "10" = "ECM related"
)
CAFs_obj_clean$subtype <- Idents(CAFs_obj_clean)

DimPlot(CAFs_obj_clean, reduction = "umap", group.by = "subtype", raster = FALSE)

DotPlot(CAFs_obj_clean,
        features = c("TAGLN", "ACTA2", "COL1A1", "COL1A2",
                     "IL6",   "CXCL12", "CCL2",  "LIF",
                     "CD74",  "SLPI",   "HLA-DRA")) +
  scale_color_viridis(option = "turbo") +
  RotatedAxis()


# -----------------------------------------------------------------------------
# 7. Further subcluster ECM-related CAFs → ECM related 1 / ECM related 2
# -----------------------------------------------------------------------------
CAFs_obj_ecm <- subset(CAFs_obj_clean, subset = subtype == "ECM related")

CAFs_obj_ecm <- FindVariableFeatures(CAFs_obj_ecm, selection.method = "vst", nfeatures = 5000)

CAFs_obj_ecm <- CAFs_obj_ecm %>%
  ScaleData(features = VariableFeatures(CAFs_obj_ecm)) %>%
  RunPCA(features = VariableFeatures(CAFs_obj_ecm), npcs = 40) %>%
  RunHarmony(group.by.vars = "orig.ident", dims.use = 1:30) %>%
  FindNeighbors(reduction = "harmony", dims = 1:30) %>%
  FindClusters(resolution = 0.1) %>%
  RunUMAP(reduction = "harmony", dims = 1:30)

DimPlot(CAFs_obj_ecm, reduction = "umap", raster = FALSE, label = TRUE)

# Marker gene expression
FeaturePlot(CAFs_obj_ecm,
            features = c("IL6", "CXCL12", "LIF"),
            cols = viridis(10, option = "turbo"))   # iCAF-like contamination check

FeaturePlot(CAFs_obj_ecm,
            features = c("TAGLN", "ACTA2", "COL1A1", "COL1A2"),
            cols = viridis(10, option = "turbo"))   # myCAF

FeaturePlot(CAFs_obj_ecm,
            features = c("LOX", "LOXL1", "LOXL2"),
            cols = viridis(10, option = "turbo"))   # remodeling

DotPlot(CAFs_obj_ecm,
        features = c("COL1A1", "COL1A2", "ACTA2", "LOX", "LOXL1", "LOXL2"),
        group.by = "seurat_clusters",
        scale    = FALSE) +
  scale_color_viridis(option = "turbo") +
  RotatedAxis()

Idents(CAFs_obj_ecm) <- "seurat_clusters"
CAFs_obj_ecm <- RenameIdents(CAFs_obj_ecm,
                             "0" = "ECM related 1",
                             "1" = "ECM related 2",
                             "2" = "ECM related 2"
)
CAFs_obj_ecm$subtype_ECM <- Idents(CAFs_obj_ecm)

DimPlot(CAFs_obj_ecm, reduction = "umap", group.by = "subtype_ECM", raster = FALSE)


# -----------------------------------------------------------------------------
# 8. DEG analysis and GO enrichment for ECM related 1 vs 2
# -----------------------------------------------------------------------------
Idents(CAFs_obj_ecm) <- "subtype_ECM"
ecm_markers <- FindAllMarkers(CAFs_obj_ecm,
                              only.pos        = TRUE,
                              logfc.threshold = 0.25)
ecm_markers <- ecm_markers %>% filter(p_val_adj < 0.05)
write.csv(ecm_markers, "output/DEG_ECM_related_CAFs.csv", row.names = FALSE)

# GO enrichment — ECM related 1
gene_ecm1 <- ecm_markers %>%
  filter(cluster == "ECM related 1") %>%
  pull(gene)

GO_ecm1 <- enrichGO(gene          = gene_ecm1,
                    keyType       = "SYMBOL",
                    OrgDb         = org.Hs.eg.db,
                    ont           = "BP",
                    pAdjustMethod = "BH",
                    qvalueCutoff  = 0.20)
GO_ecm1 <- pairwise_termsim(GO_ecm1)
dotplot(GO_ecm1, showCategory = 30)
emapplot(GO_ecm1, showCategory = 10)

# GO enrichment — ECM related 2
gene_ecm2 <- ecm_markers %>%
  filter(cluster == "ECM related 2") %>%
  pull(gene)

GO_ecm2 <- enrichGO(gene          = gene_ecm2,
                    keyType       = "SYMBOL",
                    OrgDb         = org.Hs.eg.db,
                    ont           = "BP",
                    pAdjustMethod = "BH",
                    qvalueCutoff  = 0.20)
GO_ecm2 <- pairwise_termsim(GO_ecm2)
dotplot(GO_ecm2, showCategory = 30)
emapplot(GO_ecm2, showCategory = 10)


# -----------------------------------------------------------------------------
# 9. ECM production / remodeling module scores for ECM-related CAFs
# -----------------------------------------------------------------------------
Production_genes <- c("COL1A1", "COL1A2", "COL3A1", "LUM", "DCN", "SPARC")
Remodeling_genes <- c("ACTA2", "LOX", "LOXL1", "LOXL2", "TAGLN", "POSTN")
ECM_genes        <- c(Production_genes, Remodeling_genes)

CAFs_obj_ecm <- AddModuleScore(CAFs_obj_ecm,
                               features = list(Production_genes,
                                               Remodeling_genes,
                                               ECM_genes),
                               name     = c("Productionscore",
                                            "Remodelingscore",
                                            "ECMscore"))

FeaturePlot(CAFs_obj_ecm,
            features = c("Productionscore1", "Remodelingscore2", "ECMscore3"),
            cols = viridis(10, option = "turbo"))

DotPlot(CAFs_obj_ecm,
        features = c("Productionscore1", "Remodelingscore2", "ECMscore3")) +
  scale_color_viridis(option = "turbo") +
  RotatedAxis()

# Scatter plot: Production vs Remodeling score (Z-scaled), faceted by subtype
umap_df <- Embeddings(CAFs_obj_ecm, "umap") %>%
  as.data.frame() %>%
  mutate(cluster           = CAFs_obj_ecm$subtype_ECM,
         Productionscore1  = scale(CAFs_obj_ecm$Productionscore1)[, 1],
         Remodelingscore2  = scale(CAFs_obj_ecm$Remodelingscore2)[, 1],
         ECMscore3         = scale(CAFs_obj_ecm$ECMscore3)[, 1])

cluster_means <- umap_df %>%
  group_by(cluster) %>%
  summarise(mean_prod = mean(Productionscore1, na.rm = TRUE),
            mean_remo = mean(Remodelingscore2, na.rm = TRUE))

ggplot(umap_df, aes(x = Productionscore1, y = Remodelingscore2, color = ECMscore3)) +
  geom_point(alpha = 0.8, size = 3) +
  scale_color_viridis(option = "turbo", name = "ECM score") +
  geom_vline(data = cluster_means,
             aes(xintercept = mean_prod),
             color = "black", linetype = "dashed", linewidth = 0.7) +
  geom_hline(data = cluster_means,
             aes(yintercept = mean_remo),
             color = "black", linetype = "dashed", linewidth = 0.7) +
  labs(x = "Production score (Z)", y = "Remodeling score (Z)") +
  facet_wrap(~cluster) +
  theme_bw() +
  theme(axis.text  = element_text(color = "black", face = "bold"),
        axis.title = element_text(face = "bold"),
        legend.text = element_text(color = "black", face = "bold"),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1))


# -----------------------------------------------------------------------------
# 10. Transfer ECM related 1/2 labels back to full clean CAF object and compute scores
# -----------------------------------------------------------------------------
CAFs_obj_clean$type <- as.character(CAFs_obj_clean$subtype)

cells_ecm1 <- colnames(CAFs_obj_ecm)[CAFs_obj_ecm$subtype_ECM == "ECM related 1"]
cells_ecm2 <- colnames(CAFs_obj_ecm)[CAFs_obj_ecm$subtype_ECM == "ECM related 2"]

CAFs_obj_clean$type[cells_ecm1] <- "ECM related 1"
CAFs_obj_clean$type[cells_ecm2] <- "ECM related 2"

Idents(CAFs_obj_clean) <- "type"
DimPlot(CAFs_obj_clean, reduction = "umap", raster = FALSE)

# Module scores on full clean CAF object
CAFs_obj_clean <- AddModuleScore(CAFs_obj_clean,
                                 features = list(Production_genes,
                                                 Remodeling_genes,
                                                 ECM_genes),
                                 name     = c("Productionscore",
                                              "Remodelingscore",
                                              "ECMscore"))

FeaturePlot(CAFs_obj_clean,
            features = c("Productionscore1", "Remodelingscore2", "ECMscore3"),
            cols = viridis(10, option = "turbo"))

# Scatter plot for full clean CAF object
umap_df_all <- Embeddings(CAFs_obj_clean, "umap") %>%
  as.data.frame() %>%
  mutate(cluster          = CAFs_obj_clean$type,
         Productionscore1 = scale(CAFs_obj_clean$Productionscore1)[, 1],
         Remodelingscore2 = scale(CAFs_obj_clean$Remodelingscore2)[, 1],
         ECMscore3        = scale(CAFs_obj_clean$ECMscore3)[, 1])

cluster_means_all <- umap_df_all %>%
  group_by(cluster) %>%
  summarise(mean_prod = mean(Productionscore1, na.rm = TRUE),
            mean_remo = mean(Remodelingscore2, na.rm = TRUE))

ggplot(umap_df_all, aes(x = Productionscore1, y = Remodelingscore2, color = ECMscore3)) +
  geom_point(alpha = 0.8, size = 3) +
  scale_color_viridis(option = "turbo", name = "ECM score") +
  geom_vline(data = cluster_means_all,
             aes(xintercept = mean_prod),
             color = "red", linetype = "dashed", linewidth = 0.7) +
  geom_hline(data = cluster_means_all,
             aes(yintercept = mean_remo),
             color = "red", linetype = "dashed", linewidth = 0.7) +
  geom_vline(aes(xintercept = 0), color = "black", linetype = "dashed") +
  geom_hline(aes(yintercept = 0), color = "black", linetype = "dashed") +
  labs(x = "Production score (Z)", y = "Remodeling score (Z)") +
  facet_wrap(~cluster) +
  xlim(c(0, 2)) + ylim(c(0, 2.5)) +
  theme_bw() +
  theme(axis.text   = element_text(color = "black", face = "bold"),
        axis.title  = element_text(face = "bold"),
        legend.text = element_text(color = "black", face = "bold"),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1))

saveRDS(CAFs_obj_clean, file = "output/CAFs_obj_final.rds")
