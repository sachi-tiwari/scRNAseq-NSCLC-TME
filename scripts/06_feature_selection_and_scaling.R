# =============================================================================
# Script: 06_feature_selection_and_scaling.R
# Description: Identifies highly variable features (HVGs) for dimensionality
#              reduction, visualizes them, and scales the data for PCA.
# Input:  'seurat_obj_filtered' from script 05
# Output: 'seurat_obj_filtered' with variable features identified and data scaled
#          Variable feature plot
# Note:   2000 HVGs selected using VST method — standard for scRNA-seq analysis
#         Only variable features are scaled here (sufficient for clustering and
#         subpopulation analysis). Scaling all genes would be needed for
#         cell-cell communication analysis but is not part of this pipeline.
# =============================================================================

library(Seurat)
library(ggplot2)

# --- Section 1: Identify highly variable features ---

# VST (variance stabilizing transformation) method identifies genes that vary
# more than expected given their mean expression
seurat_obj_filtered <- FindVariableFeatures(
  seurat_obj_filtered,
  selection.method = "vst",
  nfeatures = 2000
)

cat("Number of variable features identified:", 
    length(VariableFeatures(seurat_obj_filtered)), "\n")

# --- Section 2: Visualize variable features ---

# Plot all features with top 10 most variable highlighted
plot1 <- VariableFeaturePlot(seurat_obj_filtered)

top10 <- head(VariableFeatures(seurat_obj_filtered), 10)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)

# Display both plots side by side
plot1 + plot2

# --- Section 3: Scale data ---

# Scaling centers and scales expression of each gene across cells
# Required before PCA — only scaling variable features here
seurat_obj_filtered <- ScaleData(
  seurat_obj_filtered,
  features = VariableFeatures(seurat_obj_filtered)
)

cat("Data scaling complete.\n")

# --- Save checkpoint ---
saveRDS(seurat_obj_filtered, 
        file = "results/seurat_obj_filtered_scaled.rds")
cat("Saved: results/seurat_obj_filtered_scaled.rds\n")