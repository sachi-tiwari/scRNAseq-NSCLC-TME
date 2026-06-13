# =============================================================================
# Script: 08_marker_gene.R
# Description: Identifies marker genes for each cluster using FindAllMarkers
#              with Wilcoxon rank-sum test. Uses the presto package backend
#              for fast computation on large datasets (~200K cells).
# Input:  'seurat_obj_filtered' with clusters from script 07
# Output: results/all_markers_res0.5.rds — full marker gene results
#         results/top_markers.csv — top 5 markers per cluster
# Parameters:
#   only.pos = TRUE        (positive markers only)
#   min.pct = 0.25         (gene expressed in at least 25% of cells)
#   logfc.threshold = 0.25 (minimum log2 fold change)
# Note: presto package is used automatically by Seurat when installed.
#       It dramatically speeds up FindAllMarkers for large datasets.
#       Install via script 01 if not already installed.
# =============================================================================

library(Seurat)
library(dplyr)
library(presto)

# --- Section 1: Find marker genes for all clusters ---

all_markers <- FindAllMarkers(
  seurat_obj_filtered,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

cat("Marker gene identification complete.\n")
cat("Total marker genes found:", nrow(all_markers), "\n")

# Save full results
saveRDS(all_markers, file = "results/all_markers_res0.5.rds")
cat("Saved: results/all_markers_res0.5.rds\n")

# --- Section 2: Extract and save top 5 markers per cluster ---

top_markers <- all_markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC)

print(top_markers)

write.csv(top_markers, "results/top_markers.csv", row.names = TRUE)
cat("Saved: results/top_markers.csv\n")