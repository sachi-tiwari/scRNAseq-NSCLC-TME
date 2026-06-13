# =============================================================================
# Script: 11_subpopulation_analysis.R
# Description: Subsets T cells from the annotated dataset, re-clusters them
#              independently to identify T cell subpopulations, identifies
#              marker genes per subcluster, and cross-checks DEGs against
#              average expression and canonical T cell markers.
# Input:  'seurat_obj_filtered' with cell_type_cleaned from script 09
#          data/metadata/tcell_markers_refined.csv
# Output: 't_cells' — Seurat object containing only T cells, re-clustered
#         Marker genes, dot plots, UMAP, crosscheck table
# Parameters:
#   PCA dims: 1:15
#   Clustering resolution: 0.3 (14 clusters)
# =============================================================================

library(Seurat)
library(ggplot2)
library(dplyr)
library(tibble)
library(cowplot)

# --- Section 1: Subset T cells ---

t_cells <- subset(seurat_obj_filtered, subset = cell_type_cleaned == "T cells")
cat("Number of T cells:", ncol(t_cells), "\n")

# --- Section 2: Re-normalize and find variable features within T cells ---
# T cells are re-processed independently so that variable genes and PCs
# reflect variation within the T cell compartment, not the whole dataset

t_cells <- NormalizeData(t_cells)
t_cells <- FindVariableFeatures(t_cells, selection.method = "vst", nfeatures = 2000)

# --- Section 3: Scale and run PCA ---

t_cells <- ScaleData(t_cells)
t_cells <- RunPCA(t_cells, features = VariableFeatures(t_cells))

# Elbow plot to choose number of PCs
ElbowPlot(t_cells)

# --- Section 4: Clustering ---

t_cells <- FindNeighbors(t_cells, dims = 1:15)
t_cells <- FindClusters(t_cells, resolution = c(0.1, 0.3, 0.5, 0.7, 1.0))

# Resolution 0.3 selected — yields 14 biologically interpretable subclusters
DimPlot(t_cells, group.by = "RNA_snn_res.0.3", label = TRUE)
Idents(t_cells) <- "RNA_snn_res.0.3"

cat("Number of T cell subclusters:", length(unique(Idents(t_cells))), "\n")
print(table(t_cells$RNA_snn_res.0.3))

# --- Section 5: UMAP ---

t_cells <- RunUMAP(t_cells, dims = 1:15)
DimPlot(t_cells, reduction = "umap", label = TRUE, group.by = "RNA_snn_res.0.3") +
  ggtitle("T Cell Subpopulations in NSCLC")

# --- Section 6: Marker genes per subcluster ---

tcell_markers <- FindAllMarkers(
  t_cells,
  only.pos = TRUE,
  min.pct = 0.1,
  logfc.threshold = 0.1,
  test.use = "wilcox",
  min.diff.pct = 0.05
)

# Filter for statistically significant markers
tcell_markers_sig <- tcell_markers %>%
  filter(p_val_adj < 0.05) %>%
  arrange(cluster, desc(avg_log2FC))

# Top 10 markers per cluster by log2FC
top_markers <- tcell_markers %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC), .by_group = TRUE) %>%
  slice_head(n = 10) %>%
  ungroup()

print(top_markers)

# --- Section 7: Average expression per cluster ---

avg_expression <- AverageExpression(t_cells, assays = "RNA", slot = "data", 
                                    group.by = "RNA_snn_res.0.3")$RNA
avg_expression_df <- as.data.frame(avg_expression)

# For each cluster, get the top 10 genes by average expression
top_expressed_genes <- lapply(colnames(avg_expression_df), function(cluster) {
  sorted_genes <- avg_expression_df %>%
    select(all_of(cluster)) %>%
    rownames_to_column("gene") %>%
    arrange(desc(.data[[cluster]])) %>%
    slice_head(n = 10) %>%
    mutate(cluster = cluster)
  return(sorted_genes)
})

top_expressed_genes <- bind_rows(top_expressed_genes)
print(top_expressed_genes)

# --- Section 8: Canonical marker dot plot ---

marker_list <- read.csv("data/metadata/tcell_markers_refined.csv")
canonical_markers <- marker_list$gene
canonical_markers <- canonical_markers[canonical_markers %in% rownames(t_cells)]

p_dot <- DotPlot(t_cells, features = unique(canonical_markers), 
                 group.by = "RNA_snn_res.0.3") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8)) +
  ggtitle("Canonical Marker Expression Across T Cell Subclusters")
print(p_dot)

# --- Section 9: Cross-check DEGs, average expression, and canonical markers ---
# Note: AverageExpression() prefixes cluster numbers with "g" in column names
# (e.g., cluster "0" becomes column "g0")

crosscheck_results <- lapply(unique(t_cells$RNA_snn_res.0.3), function(cluster) {
  cluster <- as.character(cluster)
  cluster_col <- paste0("g", cluster)
  
  if (!(cluster_col %in% colnames(avg_expression_df))) {
    message(paste("Column", cluster_col, "not found. Skipping cluster", cluster))
    return(NULL)
  }
  # Get top 10 DEGs
  degs <- tcell_markers_sig %>%
    filter(cluster == !!cluster) %>%
    slice_head(n = 10) %>%
    pull(gene)
  # Get top 10 average expressed genes
  avg_expr <- avg_expression_df %>%
    select(all_of(cluster_col)) %>%
    rownames_to_column("gene") %>%
    arrange(desc(.data[[cluster_col]])) %>%
    slice_head(n = 10) %>%
    pull(gene)
  
  canonical <- marker_list$gene[marker_list$gene %in% rownames(t_cells)]
  
  # Find overlaps
  degs_avg_overlap <- intersect(degs, avg_expr)
  degs_canonical_overlap <- intersect(degs, canonical)
  avg_canonical_overlap <- intersect(avg_expr, canonical)
  all_overlap <- intersect(intersect(degs, avg_expr), canonical)
  
  
  
  return(data.frame(
    cluster = cluster,
    degs = paste(degs, collapse = ", "),
    avg_expr_genes = paste(avg_expr, collapse = ", "),
    degs_avg_overlap = paste(degs_avg_overlap, collapse = ", "),
    degs_canonical_overlap = paste(degs_canonical_overlap, collapse = ", "),
    avg_canonical_overlap = paste(avg_canonical_overlap, collapse = ", "),
    all_overlap = paste(all_overlap, collapse = ", ")
  ))
})


# Filter out NULL results and combine
crosscheck_results <- crosscheck_results[!sapply(crosscheck_results, is.null)]
crosscheck_df <- bind_rows(crosscheck_results)
print(crosscheck_df)

# --- Save checkpoint ---
saveRDS(t_cells, file = "results/t_cells_subclustered.rds")
write.csv(crosscheck_df, "results/tcell_crosscheck_table.csv", row.names = FALSE)
cat("Saved: results/t_cells_subclustered.rds\n")
cat("Saved: results/tcell_crosscheck_table.csv\n")