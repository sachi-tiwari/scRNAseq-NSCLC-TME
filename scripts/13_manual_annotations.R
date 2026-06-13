# =============================================================================
# Script: 13_manual_annotations.R
# Description: Assigns final manual cell type labels to T cell subclusters
#              based on marker gene analysis (scripts 11-12) and biological
#              interpretation. Maps cluster numbers (resolution 0.3) to
#              14 T cell subtype labels.
# Input:  't_cells' from script 12
# Output: 't_cells' with 'manual_annotation' metadata column
#         figures/t_cell_umap_manual.png
#         results/t_cells_final_annotated.rds
# =============================================================================

library(Seurat)
library(ggplot2)
library(dplyr)
library(plyr)      # For mapvalues
library(ggrepel)   # For non-overlapping labels

# --- Section 1: Map cluster numbers to cell type labels ---

cluster_ids <- as.character(0:13)
cell_types <- c(
  "Naive CD4+ T cells",                 # 0
  "Naive CD8+ T cells",                 # 1
  "Central memory CD4+ T cells",        # 2
  "Effector memory CD4+ T cells",       # 3
  "Effector CD8+ T cells",               # 4
  "Memory CD8+ T cells",                 # 5
  "Exhausted CD8+ T cells",              # 6
  "Progenitor exhausted CD8+ T cells",   # 7
  "Terminally exhausted CD8+ T cells",   # 8
  "IFN-γ+ Th1 cells",                   # 9
  "IL-17A+ Th17 cells",                  # 10
  "T follicular helper cells",           # 11
  "FOXP3+ Tumor-infiltrating Tregs",     # 12
  "Cycling T cells"                      # 13
)

t_cells$manual_annotation <- plyr::mapvalues(
  t_cells$RNA_snn_res.0.3,
  from = cluster_ids,
  to = cell_types
)

cat("manual_annotation distribution:\n")
print(table(t_cells$manual_annotation))

# --- Section 2: UMAP visualization with manual annotations ---

p_umap_manual <- DimPlot(t_cells, group.by = "manual_annotation", 
                         label = TRUE, label.size = 4, repel = TRUE) +
  ggtitle("Manual Annotation of T Cell Subclusters") +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        legend.position = "right",
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8)) +
  guides(colour = guide_legend(override.aes = list(size = 3)))

ggsave("figures/t_cell_umap_manual.png", plot = p_umap_manual, 
       width = 10, height = 8, dpi = 300)

cat("Saved: figures/t_cell_umap_manual.png\n")

# --- Save final annotated object ---

saveRDS(t_cells, file = "results/t_cells_final_annotated.rds")
cat("Saved: results/t_cells_final_annotated.rds\n")