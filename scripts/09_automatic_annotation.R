# =============================================================================
# Script: 09_automatic_annotation.R
# Description: Performs automated cell type annotation using SingleR with
#              two reference datasets (HPCA and Monaco Immune Data).
#              Raw SingleR labels are then harmonized into cleaned cell type
#              categories for downstream analysis.
# Input:  'seurat_obj_filtered' with UMAP from script 07
# Output: 'seurat_obj_filtered' with two new metadata columns:
#           - SingleR_multi_label: raw SingleR labels
#           - cell_type_cleaned: harmonized cell type labels
#         results/cell_type_counts.csv
#         results/reference_counts.csv
#         results/cleaned_cell_type_counts.csv
#         results/annotated_seurat_object.rds
# =============================================================================

library(Seurat)
library(SingleR)
library(celldex)
library(ggplot2)
library(scCustomize)

# --- Section 1: Prepare data for SingleR ---

# Extract normalized expression matrix from RNA assay
norm_exp_mat <- Seurat::GetAssayData(seurat_obj_filtered, 
                                     assay = "RNA", 
                                     slot = "data")

# --- Section 2: Load reference datasets ---

# HPCA: broad cell type reference (30 cell types)
ref_hpca <- celldex::HumanPrimaryCellAtlasData()

# Monaco: immune-focused reference (29 immune cell types)
ref_monaco <- celldex::MonacoImmuneData()

# --- Section 3: Run SingleR with both references ---

# SingleR assigns each cell the best matching label from either reference
multi_res <- SingleR(
  test = norm_exp_mat,
  ref = list(HPCA = ref_hpca, Monaco = ref_monaco),
  labels = list(ref_hpca$label.main, ref_monaco$label.main)
)

# Add results to Seurat metadata
seurat_obj_filtered$SingleR_multi_label <- multi_res$labels
seurat_obj_filtered$SingleR_best_ref <- multi_res$reference

cat("SingleR annotation complete.\n")
cat("Cell type distribution (raw labels):\n")
print(table(seurat_obj_filtered$SingleR_multi_label))
print(table(seurat_obj_filtered$SingleR_best_ref))

# Export raw counts
write.csv(as.data.frame(table(seurat_obj_filtered$SingleR_multi_label)),
          "results/cell_type_counts.csv")
write.csv(as.data.frame(table(seurat_obj_filtered$SingleR_best_ref)),
          "results/reference_counts.csv")

# Visualize raw SingleR annotation on UMAP
DimPlot(seurat_obj_filtered, 
        group.by = "SingleR_multi_label", 
        label = TRUE, repel = TRUE) +
  ggtitle("SingleR Cell Type Annotation (Raw Labels)")

# --- Section 4: Harmonize labels into cleaned cell type categories ---

# Map diverse SingleR labels into unified biological categories
label_map <- c(
  # T cells
  "T cells" = "T cells", "T_cells" = "T cells",
  "CD4+ T cells" = "T cells", "CD8+ T cells" = "T cells",
  "Pre-B_cell_CD34-" = "T cells", "Pro-B_cell_CD34+" = "T cells",
  # B cells
  "B cells" = "B cells", "B_cell" = "B cells",
  # Monocytes
  "Monocyte" = "Monocytes", "Monocytes" = "Monocytes",
  # NK cells
  "NK cells" = "NK cells", "NK_cell" = "NK cells",
  # Dendritic cells
  "Dendritic cells" = "Dendritic cells", "DC" = "Dendritic cells",
  # Macrophages
  "Macrophage" = "Macrophages",
  # Neutrophils
  "Neutrophils" = "Neutrophils",
  # Basophils
  "Basophils" = "Basophils",
  # Epithelial cells
  "Epithelial_cells" = "Epithelial cells",
  # Endothelial cells
  "Endothelial_cells" = "Endothelial cells",
  # Fibroblasts
  "Fibroblasts" = "Fibroblasts",
  # Smooth muscle cells
  "Smooth_muscle_cells" = "Smooth muscle cells",
  # Platelets
  "Platelets" = "Platelets",
  # Progenitors
  "Progenitors" = "Progenitors", "CMP" = "Progenitors",
  "GMP" = "Progenitors", "Myelocyte" = "Progenitors",
  "Pro-Myelocyte" = "Progenitors", "MSC" = "Progenitors",
  "HSC_CD34+" = "Progenitors", "HSC_-G-CSF" = "Progenitors",
  "BM" = "Progenitors", "BM & Prog." = "Progenitors"
)

# Apply mapping
orig_labels <- as.character(seurat_obj_filtered$SingleR_multi_label)
cleaned_labels <- label_map[orig_labels]
cleaned_labels[is.na(cleaned_labels)] <- "Other"
names(cleaned_labels) <- NULL
seurat_obj_filtered$cell_type_cleaned <- cleaned_labels

cat("Cleaned cell type distribution:\n")
print(table(seurat_obj_filtered$cell_type_cleaned))

# --- Section 5: Visualize cleaned annotations ---

p_celltypes <- DimPlot(seurat_obj_filtered, group.by = "cell_type_cleaned", label = FALSE) +
  theme(panel.background = element_rect(fill = "white"),
        plot.background = element_rect(fill = "white"))

p_celltypes_labeled <- LabelClusters(plot = p_celltypes, id = "cell_type_cleaned") +
  ggtitle("Cleaned Cell Type Annotation")

ggsave("figures/umap_celltypes_cleaned.png", plot = p_celltypes_labeled, width = 10, height = 8, dpi = 300, bg = "white")

# Proportion bar plot

p_proportion <- Proportion_Plot(
  seurat_object = seurat_obj_filtered,
  group.by = "cell_type_cleaned",
  plot_type = "bar",
  plot_scale = "percent"
) + theme(panel.background = element_rect(fill = "white"),
          plot.background = element_rect(fill = "white"))
ggsave("figures/celltype_proportions.png", plot = p_proportion, width = 10, height = 6, dpi = 300, bg = "white")

# --- Section 6: Export results ---

celltype_counts <- as.data.frame(table(seurat_obj_filtered$cell_type_cleaned))
colnames(celltype_counts) <- c("CellType", "Count")
celltype_counts$Proportion <- celltype_counts$Count / sum(celltype_counts$Count)

write.csv(celltype_counts, "results/cleaned_cell_type_counts.csv", row.names = FALSE)
cat("Saved: results/cleaned_cell_type_counts.csv\n")

# --- Save annotated object ---
saveRDS(seurat_obj_filtered, "results/annotated_seurat_object.rds")
cat("Saved: results/annotated_seurat_object.rds\n")
