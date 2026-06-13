# =============================================================================
# Script: 04_create_obj_and_QC.R
# Description: Aligns metadata with count matrix, creates a Seurat object,
#              calculates QC metrics, and visualizes distributions to help
#              choose filtering thresholds.
# Input:  'cts' and 'metadata' objects from script 03
# Output: 'seurat_obj' Seurat object with QC metrics added
#          QC plots (violin plots, scatter plots, histograms)
# Note:   Expected dimensions before filtering: ~56,283 genes x 224,611 cells
# =============================================================================

library(Seurat)
library(ggplot2)

# --- Section 1: Align metadata and count matrix ---

# Set metadata row names to cell barcodes
rownames(metadata) <- metadata$X

# Keep only cells present in both count matrix and metadata
common_cells <- intersect(colnames(cts), metadata$X)
cts <- cts[, common_cells]
metadata <- metadata[common_cells, ]

# Verify alignment — this will throw an error if something is wrong
stopifnot(all(colnames(cts) == rownames(metadata)))
cat("✅ Metadata and count matrix are aligned.\n")
cat("Number of common cells:", length(common_cells), "\n")

# --- Section 2: Create Seurat object ---

seurat_obj <- CreateSeuratObject(
  counts = cts,
  meta.data = metadata,
  project = "NSCLC",
  min.cells = 3,      # Keep genes detected in at least 3 cells
  min.features = 200  # Keep cells with at least 200 detected genes
)

cat("Seurat object created:\n")
print(seurat_obj)

# --- Section 3: Calculate QC metrics ---

# Mitochondrial gene percentage (high % indicates dying/damaged cells)
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")

# View metadata in RStudio viewer (optional)
View(seurat_obj@meta.data)

# --- Section 4: Visualize QC metrics ---

# Violin plots — useful for smaller datasets
# For large datasets (~200K cells) histograms below are more informative
VlnPlot(seurat_obj, 
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        ncol = 3)

# Scatter plots to check relationships between QC metrics
FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  geom_smooth(method = "lm")

FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "percent.mt") +
  geom_smooth(method = "lm")

FeatureScatter(seurat_obj, feature1 = "nFeature_RNA", feature2 = "percent.mt") +
  geom_smooth(method = "lm")

# Histograms with density overlays — better for large datasets
meta <- seurat_obj@meta.data

# nFeature_RNA distribution
ggplot(meta, aes(x = nFeature_RNA)) +
  geom_histogram(bins = 100, fill = "skyblue", color = "black", alpha = 0.7) +
  geom_density(aes(y = after_stat(scaled) * max(after_stat(count))),
               color = "red", linewidth = 1) +
  theme_minimal() +
  labs(title = "nFeature_RNA Distribution", x = "nFeature_RNA", y = "Cell Count")

# nCount_RNA distribution
ggplot(meta, aes(x = nCount_RNA)) +
  geom_histogram(bins = 100, fill = "orange", color = "black", alpha = 0.7) +
  geom_density(aes(y = after_stat(scaled) * max(after_stat(count))),
               color = "red", linewidth = 1) +
  theme_minimal() +
  labs(title = "nCount_RNA Distribution", x = "nCount_RNA", y = "Cell Count")

# Mitochondrial % distribution
ggplot(meta, aes(x = percent.mt)) +
  geom_histogram(bins = 100, fill = "lightgreen", color = "black", alpha = 0.7) +
  geom_density(aes(y = after_stat(scaled) * max(after_stat(count))),
               color = "red", linewidth = 1) +
  theme_minimal() +
  labs(title = "Mitochondrial % Distribution", x = "percent.mt", y = "Cell Count")

# --- Section 5: Choose thresholds ---
# Based on the histograms above, thresholds used in script 05:
#   nFeature_RNA: 200 - 3000
#   nCount_RNA:   500 - 25000
#   percent.mt:   < 5%