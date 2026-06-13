# Pipeline Overview & Methods Notes

This document describes the analytical decisions made at each step of the scRNA-seq pipeline.

## Data Source

The dataset was downloaded from Figshare (public repository):
- Count matrix RDS: https://figshare.com/ndownloader/files/39537094
- Metadata CSV: https://figshare.com/ndownloader/files/39537250

Initial dimensions: ~224,611 cells x 56,283 genes

## Step 1 - Quality Control (Script 04)

QC metrics calculated:
- nFeature_RNA: number of unique genes detected per cell
- nCount_RNA: total UMI counts per cell
- percent.mt: percentage of reads mapping to mitochondrial genes

Filtering thresholds applied:
- nFeature_RNA: 200 - 3000
- nCount_RNA: 500 - 25000
- percent.mt: < 5%

These thresholds remove likely empty droplets, doublets, and dying/damaged cells.


## Step 2 - Normalization (Script 05)

NormalizeData() with LogNormalize method, scale factor 10000 (Seurat defaults). Normalizes each cell's total counts to 10,000 and log-transforms the result.

## Step 3 - Feature Selection & Scaling (Script 06)

2000 highly variable genes identified using the VST method. ScaleData() applied to these variable features for PCA.

## Step 4 - PCA, Clustering, UMAP (Script 07)

- PCA run on 2000 variable features
- ElbowPlot used to select 15 PCs
- Multiple resolutions tested (0.1, 0.3, 0.5, 0.7, 1.0)
- clustree used to visualize cluster stability across resolutions - resolution 0.5 selected
- UMAP run with dims 1:15 for visualization

## Step 5 - Marker Gene Identification (Script 08)

FindAllMarkers() run with:
- only.pos = TRUE
- min.pct = 0.25
- logfc.threshold = 0.25

presto package used as backend for fast Wilcoxon test (essential for ~200K cells). Top 5 markers per cluster exported to top_markers.csv.

## Step 6 - Automated Annotation (Script 09)

SingleR used with two reference datasets:
1. HumanPrimaryCellAtlasData (HPCA) - broad cell type reference
2. MonacoImmuneData - immune-focused reference

Raw labels harmonized into cleaned categories using a manual mapping dictionary (cell_type_cleaned).

## Step 7 - Annotation Diagnostics (Script 10)

- Pruned labels: SingleR's conservative annotation, ambiguous cells set to NA
- Delta distribution: confidence margin between top and second-best annotation scores
- Marker support scoring: checks what fraction of cells in each annotated cluster express expected canonical markers (canonical_marker.csv, manually curated from literature)
  - fraction_expressing > 0.3 and avg_expression > 0.7 = "Strong" marker support

## Step 8 - T Cell Subpopulation Analysis (Script 11)

T cells subsetted and re-analyzed independently:
- Re-normalized and re-ran variable feature selection within T cells only
- PCA + UMAP (dims 1:15)
- Resolution 0.3 selected - 14 clusters
- Marker genes identified per subcluster
- Cross-check: DEGs vs average expression vs canonical marker list overlap assessed

## Step 9 - Scoring & Manual Annotation (Scripts 12-13)

AddModuleScore used to score each cell against canonical marker gene sets (tcell_markers_refined.csv), assigning a scored_cell_type based on highest score per cell.

Final manual annotation (script 13) assigns biologically interpretable labels to each of the 14 clusters based on marker gene analysis and cross-checking from scripts 11-12:

- Exhaustion trajectory: Naive -> Progenitor exhausted -> Exhausted -> Terminally exhausted CD8+
- Helper subsets: Th1 (IFN-gamma+), Th17 (IL-17A+), Tfh
- Immunosuppressive: FOXP3+ Tumor-infiltrating Tregs
- Proliferating: Cycling T cells (MKI67+, TOP2A+)

## Notes for Reproducibility

- RDS checkpoints saved at several steps (after PCA, after UMAP, after annotation) to allow resuming without re-running the full pipeline
- .RData environment saves and large .rds files are gitignored due to size
- canonical_marker.csv (broad cell types) is manually curated from established immunology literature, not auto-generated
