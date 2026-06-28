library(readr)
library(tidyverse)
library(immunedeconv)

##################### Load data ##############################

# Gene expression matrix after probe selection
#
# Input file:
# data/after_probe_selection_gene_expression_data_processed.csv
#
# Rows: genes
# Columns: samples

Mat <-
  read_csv(
    "data/after_probe_selection_gene_expression_data_processed.csv"
  ) %>%
  column_to_rownames("gene") %>%
  as.matrix()

##################### QuanTIseq ##############################

res_quan <-
  deconvolute_quantiseq(
    gene_expression_matrix = Mat,
    tumor = TRUE,
    arrays = TRUE,
    scale_mrna = TRUE
  )

##################### EPIC ##############################

res_EPIC <-
  deconvolute_epic(
    gene_expression_matrix = Mat,
    tumor = TRUE,
    scale_mrna = TRUE
  )

##################### MCP counter ##############################

res_MCP <-
  deconvolute_mcp_counter(
    gene_expression_matrix = Mat
  )
