# =============================================================================
# Script: 12_subpopulation_scoring.R
# Description: Scores T cell subclusters against canonical marker gene sets
#              using AddModuleScore to assign biologically-informed subtype
#              labels. Validates scored annotations with UMAP and dot plot.
# Input:  't_cells' and 'marker_list' from script 11
#          marker_list must contain columns: gene, cell_type
# Output: 't_cells' with 'scored_cell_type' metadata column
#         figures/t_cell_umap_scored_nolabels.png
#         figures/t_cell_dot_validation_refined.png
# =============================================================================

library(Seurat)
library(ggplot2)
library(dplyr)
library(tibble)
library(cowplot)
library(ggrepel)

# --- Section 1: Verify marker list structure ---

cat("Checking marker_list structure:\n")
if (!all(c("gene", "cell_type") %in% colnames(marker_list))) {
  stop("marker_list must contain 'gene' and 'cell_type' columns.")
}
print(head(marker_list))

# --- Section 2: Create gene sets for scoring ---

gene_sets <- split(marker_list$gene, marker_list$cell_type)
gene_sets <- lapply(gene_sets, function(x) x[x %in% rownames(t_cells)])

# --- Section 3: Score cells for each T cell subtype ---

t_cells <- AddModuleScore(t_cells, features = gene_sets, name = "SubtypeScore", ctrl = 5)

# --- Section 4: Assign cell types based on highest score ---

score_columns <- paste0("SubtypeScore", seq_along(gene_sets))
t_cells$scored_cell_type <- apply(t_cells@meta.data[score_columns], 1, function(x) {
  if (all(is.na(x))) return("Unknown")
  names(gene_sets)[which.max(x)]
})

cat("scored_cell_type distribution:\n")
print(table(t_cells$scored_cell_type))

# --- Section 5: UMAP visualization ---

p_umap_scored <- DimPlot(t_cells, reduction = "umap", group.by = "scored_cell_type", label = FALSE) +
  ggtitle("T Cell Subpopulations with Scored Annotations") +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        legend.position = "right",
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8)) +
  guides(colour = guide_legend(override.aes = list(size = 3)))

ggsave("figures/t_cell_umap_scored_nolabels.png", plot = p_umap_scored, 
       width = 10, height = 8, dpi = 300)

# --- Section 6: Dot plot validation against canonical markers ---

canonical_markers <- marker_list$gene[marker_list$gene %in% rownames(t_cells)]
cat("Number of canonical markers in t_cells:", length(canonical_markers), "\n")
if (length(canonical_markers) == 0) {
  stop("No canonical markers found in t_cells. Check marker_list.")
}

p_dot_validation <- DotPlot(t_cells, features = head(canonical_markers, 20), group.by = "scored_cell_type") +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        axis.text.y = element_text(size = 8),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        legend.position = "right",
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8),
        panel.background = element_rect(fill = "white"),
        plot.background = element_rect(fill = "white")) +
  ggtitle("Canonical Marker Expression Across Scored T Cell Subpopulations") +
  guides(size = guide_legend(title = "Percent Expressed"), color = guide_legend(title = "Average Expression"))

ggsave("figures/t_cell_dot_validation_refined.png", plot = p_dot_validation, 
       width = 12, height = 8, dpi = 300)

cat("Saved: figures/t_cell_umap_scored_nolabels.png\n")
cat("Saved: figures/t_cell_dot_validation_refined.png\n")