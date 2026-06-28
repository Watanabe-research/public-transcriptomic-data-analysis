# =============================================================================
# Microarray Data Preprocessing
# Methods:
#   - MAS5-normalized gene expression data, already log2-transformed
#   - Probes selected based on Affymetrix annotation reliability
#     (priority: "_at" > "_s_at" > "_x_at")
#   - When multiple probes map to the same gene, the probe with the highest
#     average expression across all samples is selected
# =============================================================================


# -----------------------------------------------------------------------------
# 1. Load libraries
# -----------------------------------------------------------------------------
library(tidyverse)


# -----------------------------------------------------------------------------
# 2. Load data
#    patients_and_t_log2_mRNA_data.csv:
#      - Columns 1-4 : patient metadata
#      - Columns 5-  : MAS5-normalized, log2-transformed probe expression values
#    gene_name.csv:
#      - Affymetrix probe annotation file containing probe ID and Gene Symbol
# -----------------------------------------------------------------------------
patients_and_t_log2_mRNA_data <- read_csv("data/patients_and_t_log2_mRNA_data.csv")
gene_name                      <- read_csv("data/gene_name.csv")


# -----------------------------------------------------------------------------
# 3. Separate patient metadata and expression matrix
# -----------------------------------------------------------------------------
pati_info <- patients_and_t_log2_mRNA_data %>%
  dplyr::select(1:4) %>%
  mutate(patient = substr(patients_and_t_log2_mRNA_data$Title, 1,
                          nchar(patients_and_t_log2_mRNA_data$Title) - 1))

# Expression matrix: probes as rows, samples as columns
n_meta_cols <- 4
mRNA <- patients_and_t_log2_mRNA_data[ , (n_meta_cols + 1):ncol(patients_and_t_log2_mRNA_data)] %>%
  t()
colnames(mRNA) <- pati_info$ID


# -----------------------------------------------------------------------------
# 4. Compute average expression across all samples (used for probe selection)
# -----------------------------------------------------------------------------
mRNA_df <- mRNA %>%
  as.data.frame() %>%
  mutate(ave = rowMeans(across(everything()))) %>%
  rownames_to_column(var = "probe_ID")


# -----------------------------------------------------------------------------
# 5. Merge with Affymetrix probe annotation
# -----------------------------------------------------------------------------
probe_annotation <- gene_name %>%
  dplyr::select("ID", "Gene Symbol") %>%
  rename(probe_ID    = "ID",
         Gene_Symbol = "Gene Symbol") %>%
  # When a probe maps to multiple genes (separated by " /// "), retain the first
  mutate(Gene_Symbol = sapply(strsplit(Gene_Symbol, " /// "), `[`, 1))

df <- merge(mRNA_df, probe_annotation, by = "probe_ID") %>%
  na.omit()


# -----------------------------------------------------------------------------
# 6. Select best probe per gene
#    Priority: "_at" (most reliable) > "_s_at" > "_x_at"
#    Among probes of the same priority tier, the one with the highest average
#    expression is selected.
# -----------------------------------------------------------------------------
select_best_probe <- function(data) {
  data %>%
    group_by(Gene_Symbol) %>%
    reframe(
      selected_probe = {
        if (any(grepl("_at$", probe_ID))) {
          probe_ID[grepl("_at$", probe_ID)][which.max(ave[grepl("_at$", probe_ID)])]
        } else if (any(grepl("_s_at$", probe_ID))) {
          probe_ID[grepl("_s_at$", probe_ID)][which.max(ave[grepl("_s_at$", probe_ID)])]
        } else if (any(grepl("_x_at$", probe_ID))) {
          probe_ID[grepl("_x_at$", probe_ID)][which.max(ave[grepl("_x_at$", probe_ID)])]
        } else {
          NA_character_
        }
      }
    ) %>%
    left_join(data, by = c("selected_probe" = "probe_ID"))
}

df_cleaned <- select_best_probe(df)
write_csv(df_cleaned, "output/after_probe_selection_gene_expression_data.csv")
