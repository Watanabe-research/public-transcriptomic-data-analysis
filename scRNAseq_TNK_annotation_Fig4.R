# =============================================================================
# T/NK Cell Subclustering and Annotation
# Methods: SingleR with BlueprintEncodeData and MonacoImmuneData as references.
#          Only cells assigned the same cell type in both references are retained.
# =============================================================================


# -----------------------------------------------------------------------------
# 1. Load libraries
# -----------------------------------------------------------------------------
library(tidyverse)
library(Seurat)
library(harmony)
library(SingleR)
library(celldex)
library(SingleCellExperiment)
library(viridis)


# -----------------------------------------------------------------------------
# 2. Load data and add tissue type metadata
# -----------------------------------------------------------------------------
all_cells <- readRDS("data/all_cells.rds")

all_cells$TissueType <- ifelse(grepl("^N", all_cells$orig.ident), "Normal", "Tumor")


# -----------------------------------------------------------------------------
# 3. Subset T/NK cells
# -----------------------------------------------------------------------------
TNK_obj <- subset(all_cells, subset = new_label == "T/NK cells")


# -----------------------------------------------------------------------------
# 4. Preprocessing and dimensionality reduction
# -----------------------------------------------------------------------------
TNK_obj <- NormalizeData(TNK_obj)
TNK_obj <- FindVariableFeatures(TNK_obj, selection.method = "vst", nfeatures = 2000)
TNK_obj <- ScaleData(TNK_obj)
TNK_obj <- RunPCA(TNK_obj, features = VariableFeatures(TNK_obj), npcs = 40)

ElbowPlot(TNK_obj, ndims = 40)  # check elbow to confirm dims = 30


# -----------------------------------------------------------------------------
# 5. Harmony batch correction, clustering, and UMAP
# -----------------------------------------------------------------------------
TNK_obj <- RunHarmony(TNK_obj, group.by.vars = "orig.ident", dims.use = 1:30)

TNK_obj_res <- TNK_obj %>%
  FindNeighbors(reduction = "harmony", dims = 1:30) %>%
  FindClusters(resolution = 1.5) %>%
  RunUMAP(reduction = "harmony", dims = 1:30)

DimPlot(TNK_obj_res, reduction = "umap", raster = FALSE, label = TRUE) +
  ggtitle("T/NK cells: nfeatures = 2000, dims = 30")


# -----------------------------------------------------------------------------
# 6. SingleR annotation — Monaco ImmuneData and BlueprintEncodeData
# -----------------------------------------------------------------------------

# --- Monaco ImmuneData ---
ref_mona  <- celldex::MonacoImmuneData()
results_mona <- SingleR(
  test   = as.SingleCellExperiment(TNK_obj_res),
  ref    = ref_mona,
  labels = ref_mona$label.fine
)
TNK_obj_res$singlr_labels_mona <- results_mona$labels

# --- BlueprintEncodeData ---
ref_blue  <- celldex::BlueprintEncodeData()
results_blue <- SingleR(
  test   = as.SingleCellExperiment(TNK_obj_res),
  ref    = ref_blue,
  labels = ref_blue$label.fine
)
TNK_obj_res$singlr_labels_blue <- results_blue$labels

# Visualise raw SingleR labels
DimPlot(TNK_obj_res, reduction = "umap", group.by = "singlr_labels_mona", label = FALSE) +
  ggtitle("SingleR: Monaco")
DimPlot(TNK_obj_res, reduction = "umap", group.by = "singlr_labels_blue", label = FALSE) +
  ggtitle("SingleR: Blueprint")
DimPlot(TNK_obj_res, reduction = "umap", group.by = "TissueType",         label = FALSE)


# -----------------------------------------------------------------------------
# 7. Map fine labels → broad subtypes
# -----------------------------------------------------------------------------

# --- Monaco ---
Idents(TNK_obj_res) <- "singlr_labels_mona"
TNK_obj_res <- RenameIdents(TNK_obj_res,
  "Terminal effector CD4 T cells" = "CD4 T cells",
  "Naive CD4 T cells"             = "CD4 T cells",
  "Follicular helper T cells"     = "CD4 T cells",
  "Th1/Th17 cells"                = "CD4 T cells",
  "Th1 cells"                     = "CD4 T cells",
  "Th2 cells"                     = "CD4 T cells",
  "Th17 cells"                    = "CD4 T cells",
  "Terminal effector CD8 T cells" = "CD8 T cells",
  "Central memory CD8 T cells"    = "CD8 T cells",
  "Effector memory CD8 T cells"   = "CD8 T cells",
  "Naive CD8 T cells"             = "CD8 T cells",
  "Natural killer cells"          = "NK cells",
  "T regulatory cells"            = "Treg",
  "Vd2 gd T cells"                = "Other T cells",
  "Non-Vd2 gd T cells"            = "Other T cells",
  "MAIT cells"                    = "Other T cells",
  # non-T/NK populations
  "Progenitor cells"              = "Not T_NK",
  "Myeloid dendritic cells"       = "Not T_NK",
  "Intermediate monocytes"        = "Not T_NK",
  "Plasmacytoid dendritic cells"  = "Not T_NK",
  "Exhausted B cells"             = "Not T_NK",
  "Non classical monocytes"       = "Not T_NK",
  "Plasmablasts"                  = "Not T_NK",
  "Switched memory B cells"       = "Not T_NK",
  "Non-switched memory B cells"   = "Not T_NK",
  "Naive B cells"                 = "Not T_NK",
  "Classical monocytes"           = "Not T_NK",
  "Low-density basophils"         = "Not T_NK"
)
TNK_obj_res$subtype_mona <- Idents(TNK_obj_res)

# --- Blueprint ---
Idents(TNK_obj_res) <- "singlr_labels_blue"
TNK_obj_res <- RenameIdents(TNK_obj_res,
  "CD4+ T-cells"              = "CD4 T cells",
  "CD4+ Tcm"                  = "CD4 T cells",
  "CD4+ Tem"                  = "CD4 T cells",
  "CD8+ T-cells"              = "CD8 T cells",
  "CD8+ Tcm"                  = "CD8 T cells",
  "CD8+ Tem"                  = "CD8 T cells",
  "NK cells"                  = "NK cells",
  "Tregs"                     = "Treg",
  # non-T/NK populations
  "CLP"                       = "Not T_NK",
  "CMP"                       = "Not T_NK",
  "GMP"                       = "Not T_NK",
  "MEP"                       = "Not T_NK",
  "HSC"                       = "Not T_NK",
  "Monocytes"                 = "Not T_NK",
  "Macrophages"               = "Not T_NK",
  "Macrophages M1"            = "Not T_NK",
  "Macrophages M2"            = "Not T_NK",
  "DC"                        = "Not T_NK",
  "Plasma cells"              = "Not T_NK",
  "Class-switched memory B-cells" = "Not T_NK",
  "Memory B-cells"            = "Not T_NK",
  "naive B-cells"             = "Not T_NK",
  "Megakaryocytes"            = "Not T_NK",
  "Endothelial cells"         = "Not T_NK",
  "mv Endothelial cells"      = "Not T_NK",
  "Epithelial cells"          = "Not T_NK",
  "Fibroblasts"               = "Not T_NK",
  "Adipocytes"                = "Not T_NK",
  "Mesangial cells"           = "Not T_NK",
  "Melanocytes"               = "Not T_NK",
  "Keratinocytes"             = "Not T_NK",
  "Neurons"                   = "Not T_NK",
  "Astrocytes"                = "Not T_NK",
  "Skeletal muscle"           = "Not T_NK",
  "Chondrocytes"              = "Not T_NK",
  "Erythrocytes"              = "Not T_NK"
)
TNK_obj_res$subtype_blue <- Idents(TNK_obj_res)


# -----------------------------------------------------------------------------
# 8. Consensus annotation: retain cells with matching labels in both references
# -----------------------------------------------------------------------------
TNK_obj_res$common_TNK_labels <- ifelse(
  as.character(TNK_obj_res$subtype_mona) == as.character(TNK_obj_res$subtype_blue),
  as.character(TNK_obj_res$subtype_mona),
  "not assigned"
)

DimPlot(TNK_obj_res, reduction = "umap", group.by = "common_TNK_labels", label = FALSE) +
  ggtitle("Consensus annotation (Monaco ∩ Blueprint)")


# -----------------------------------------------------------------------------
# 9. QC plots for consensus-annotated cells
# -----------------------------------------------------------------------------
target_labels <- c("CD4 T cells", "CD8 T cells", "NK cells", "Treg")

Idents(TNK_obj_res) <- "common_TNK_labels"
sub_obj <- subset(TNK_obj_res, idents = target_labels)

# Dot plot — canonical markers
DotPlot(sub_obj,
        features      = c("PTPRC", "ITGAM", "CD3D", "CD4",
                          "CD8A", "KLRD1", "KLRC1", "FOXP3", "IL2RA"),
        group.by      = "common_TNK_labels") +
  scale_color_viridis(option = "turbo") +
  RotatedAxis() +
  ggtitle("Canonical marker expression by consensus label")

# Cell-number bar chart
ggplot(as.data.frame(table(TNK_obj_res$common_TNK_labels)),
       aes(x = Var1, y = Freq, fill = Var1)) +
  geom_bar(stat = "identity") +
  ylab("Cell number") +
  xlab("") +
  theme_bw() +
  theme(legend.position = "none") +
  ggtitle("Cell counts per consensus label")


# -----------------------------------------------------------------------------
# 10. Save
# -----------------------------------------------------------------------------
saveRDS(TNK_obj_res, "output/TNK_obj_res.rds")