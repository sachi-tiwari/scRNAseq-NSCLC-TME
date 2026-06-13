# =============================================================================
# Script: 07_PCA_clustering_umap.R
# Description: Performs PCA for linear dimensionality reduction, determines
#              optimal number of PCs using elbow plot, clusters cells across
#              multiple resolutions, selects optimal resolution using clustree,
#              and runs UMAP for 2D visualization.
# Input:  'seurat_obj_filtered' from script 06
# Output: 'seurat_obj_filtered' with PCA, clusters, and UMAP embeddings
#          Saved RDS checkpoints in results/
#          Elbow plot, clustree plot, UMAP plot
# Parameters used:
#   PCA dims for clustering/UMAP: 1:15 (determined from elbow plot)
#   Final clustering resolution: 0.5 (determined from clustree analysis)
# =============================================================================

library(Seurat)
library(ggplot2)
library(clustree)

# --- Section 1: PCA ---

seurat_obj_filtered <- RunPCA(
  seurat_obj_filtered,
  features = VariableFeatures(seurat_obj_filtered)
)

# Save checkpoint after PCA — PCA is computationally expensive
# Load this checkpoint to resume from here without rerunning scripts 01-06
saveRDS(seurat_obj_filtered, file = "results/seurat_obj_filtered_afterPCA.rds")
# seurat_obj_filtered <- readRDS("results/seurat_obj_filtered_afterPCA.rds")

# Free memory by removing the unfiltered object (no longer needed)
rm(seurat_obj)
gc()

# --- Section 2: Inspect PCA results ---

cat("Number of PCs computed:", length(seurat_obj_filtered[["pca"]]@stdev), "\n")

# Top 5 features for first 5 PCs
print(seurat_obj_filtered[["pca"]], dims = 1:5, nfeatures = 5)

# Visualize PCA loadings
VizDimLoadings(seurat_obj_filtered, dims = 1:2, reduction = "pca")

# Plot cells in PCA space
DimPlot(seurat_obj_filtered, reduction = "pca")

# Elbow plot — used to determine how many PCs capture meaningful variation
# Look for the "elbow" where standard deviation levels off
# Based on this plot, dims 1:15 were selected
ElbowPlot(seurat_obj_filtered)

# --- Section 3: Clustering ---

# Build k-nearest neighbor graph using top 15 PCs
seurat_obj_filtered <- FindNeighbors(seurat_obj_filtered, dims = 1:15)

# Test multiple resolutions to find optimal clustering granularity
seurat_obj_filtered <- FindClusters(
  seurat_obj_filtered,
  resolution = c(0.1, 0.3, 0.5, 0.7, 1.0)
)

# Use clustree to visualize cluster stability across resolutions
# Each node = a cluster; edges show how clusters split/merge as resolution increases
# Choose a resolution where clusters are stable (minimal splitting/merging)
clustree(seurat_obj_filtered@meta.data, prefix = "RNA_snn_res.")

# Resolution 0.5 selected — stable clusters with biologically meaningful separation
Idents(seurat_obj_filtered) <- "RNA_snn_res.0.5"
cat("Active identity set to resolution 0.5\n")
cat("Number of clusters:", length(unique(Idents(seurat_obj_filtered))), "\n")

# --- Section 4: UMAP ---

# UMAP for 2D visualization of clusters
seurat_obj_filtered <- RunUMAP(seurat_obj_filtered, dims = 1:15)

# Visualize clusters on UMAP
DimPlot(seurat_obj_filtered, reduction = "umap", label = TRUE) +
  ggtitle("UMAP — Clusters at Resolution 0.5")

# Save checkpoint after UMAP
saveRDS(seurat_obj_filtered, file = "results/seurat_obj_filtered_afterUMAP.rds")
cat("Saved: results/seurat_obj_filtered_afterUMAP.rds\n")