# =============================================================================
# Script: 03_import_and_preview_data.R
# Description: Loads the NSCLC count matrix and metadata into R and performs
#              basic checks to confirm data integrity before analysis.
# Data: An integrated single-cell transcriptomic dataset for NSCLC
# Source: https://figshare.com/collections/An_integrated_single-cell_
#         transcriptomic_dataset_for1_non-small_cell_lung_cancer/6222221/3
# Input:  raw_data/RNA_rawcounts_matrix.rds
#         raw_data/metadata.csv
# Output: Objects 'cts' (count matrix) and 'metadata' loaded into R environment
# =============================================================================

library(Seurat)
library(tidyverse)

# --- Load data ---
cts <- readRDS("raw_data/RNA_rawcounts_matrix.rds")
metadata <- read.csv("raw_data/metadata.csv", stringsAsFactors = FALSE)

# --- Preview structure ---
cat("Dimensions of count matrix (genes x cells):\n")
print(dim(cts))

cat("\nFirst few cell barcodes:\n")
print(head(colnames(cts)))

cat("\nFirst few rows of metadata:\n")
print(head(metadata))

cat("\nNumber of cells in count matrix:", ncol(cts), "\n")
cat("Number of cells in metadata:", nrow(metadata), "\n")

# --- Verify cell count match ---
if (ncol(cts) == nrow(metadata)) {
  cat("✅ Cell numbers match between count matrix and metadata.\n")
} else {
  cat("⚠️ Cell numbers do NOT match! Please check your files.\n")
}