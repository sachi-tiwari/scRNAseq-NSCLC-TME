# =============================================================================
# Script: 05_filter_and_normalize.R
# Description: Filters low quality cells based on QC thresholds determined
#              in script 04, then normalizes the filtered data using
#              log normalization.
# Input:  'seurat_obj' from script 04
# Output: 'seurat_obj_filtered' — filtered and normalized Seurat object
#          Histogram plots showing QC distributions after filtering
# Filtering thresholds applied:
#   nFeature_RNA: 200 - 3000  (removes empty droplets and doublets)
#   nCount_RNA:   500 - 25000 (removes low quality and multiplet cells)
#   percent.mt:   < 5%        (removes dying/damaged cells)
# =============================================================================

library(Seurat)
library(ggplot2)

# --- Section 1: Filter low quality cells ---

seurat_obj_filtered <- subset(seurat_obj,
                              subset = nFeature_RNA > 200 & nFeature_RNA < 3000 &
                                nCount_RNA > 500 & nCount_RNA < 25000 &
                                percent.mt < 5)

# Compare cell numbers before and after filtering
cat("Before filtering:\n")
print(seurat_obj)
cat("\nAfter filtering:\n")
print(seurat_obj_filtered)

# --- Section 2: Visualize QC metrics after filtering ---

meta <- seurat_obj_filtered@meta.data

# nFeature_RNA distribution after filtering
ggplot(meta, aes(x = nFeature_RNA)) +
  geom_histogram(bins = 100, fill = "skyblue", color = "black", alpha = 0.7) +
  geom_density(aes(y = after_stat(scaled) * max(after_stat(count))),
               color = "red", linewidth = 1) +
  theme_minimal() +
  labs(title = "nFeature_RNA Distribution (Post-filtering)",
       x = "nFeature_RNA", y = "Cell Count")

# nCount_RNA distribution after filtering
ggplot(meta, aes(x = nCount_RNA)) +
  geom_histogram(bins = 100, fill = "orange", color = "black", alpha = 0.7) +
  geom_density(aes(y = after_stat(scaled) * max(after_stat(count))),
               color = "red", linewidth = 1) +
  theme_minimal() +
  labs(title = "nCount_RNA Distribution (Post-filtering)",
       x = "nCount_RNA", y = "Cell Count")

# Mitochondrial % distribution after filtering
ggplot(meta, aes(x = percent.mt)) +
  geom_histogram(bins = 100, fill = "lightgreen", color = "black", alpha = 0.7) +
  geom_density(aes(y = after_stat(scaled) * max(after_stat(count))),
               color = "red", linewidth = 1) +
  theme_minimal() +
  labs(title = "Mitochondrial % Distribution (Post-filtering)",
       x = "percent.mt", y = "Cell Count")

# --- Section 3: Normalize data ---

# Log normalization: normalizes each cell to total count of 10,000
# then log-transforms. This is the standard method for scRNA-seq data.
# Using default settings: normalization.method = "LogNormalize", scale.factor = 10000

seurat_obj_filtered <- NormalizeData(seurat_obj_filtered)

cat("Normalization complete.\n")

# --- Save checkpoint ---
saveRDS(seurat_obj_filtered, file = "results/seurat_obj_filtered_normalized.rds")
cat("Saved: results/seurat_obj_filtered_normalized.rds\n")