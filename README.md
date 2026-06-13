# scRNA-seq Analysis of the NSCLC Tumor Microenvironment

A complete single-cell RNA-sequencing pipeline for characterizing immune cell populations and T cell exhaustion states in non-small cell lung cancer (NSCLC).

## Project Overview

This project performs scRNA-seq analysis of the NSCLC tumor microenvironment (TME), focusing on:

- Profiling major immune and non-immune cell populations
- Deep characterization of T cell subpopulations, including activation and exhaustion states
- Identifying marker genes for each annotated cell type
- Applying both automated (SingleR) and manual annotation strategies

## Repository Structure

```
scRNAseq-NSCLC-TME/

scripts/          - R analysis scripts (numbered in pipeline order)

data/

metadata/        - Cell metadata and canonical marker gene lists

raw_data/        - Raw count matrix (not tracked, see download script)

results/          - Output tables (marker genes, cell type counts)

figures/          - Output plots (UMAPs, dotplots, feature plots)

docs/             - Detailed methods documentation

README.md

```
## Dataset

- Source: Publicly available NSCLC scRNA-seq dataset hosted on Figshare
  - Count matrix: https://figshare.com/ndownloader/files/39537094
  - Metadata: https://figshare.com/ndownloader/files/39537250
- Scale: ~224,611 cells x 56,283 genes (before QC filtering)

Raw data files are not included in this repository due to size. Run `scripts/02_download_data.R` to download them automatically.

## Pipeline Summary

1. **Quality Control & Filtering** - Removed low-quality cells based on gene counts, UMI counts, and mitochondrial percentage
2. **Normalization** - Log normalization (LogNormalize, scale factor 10,000)
3. **Feature Selection & Scaling** - 2,000 highly variable genes selected using VST method
4. **PCA, Clustering, UMAP** - 15 PCs, resolution 0.5, 2D visualization
5. **Marker Gene Identification** - FindAllMarkers with Wilcoxon test
6. **Cell Type Annotation** - SingleR with HPCA and Monaco Immune references
7. **Annotation Diagnostics** - Confidence assessment using pruned labels and marker validation
8. **T Cell Subpopulation Analysis** - Re-clustered T cells (resolution 0.3, 14 subclusters)
9. **Scoring & Manual Annotation** - AddModuleScore-based scoring and final manual labeling

## T Cell Subpopulations Identified

| Cluster | Cell Type |
|---------|-----------|
| 0 | Naive CD4+ T cells |
| 1 | Naive CD8+ T cells |
| 2 | Central memory CD4+ T cells |
| 3 | Effector memory CD4+ T cells |
| 4 | Effector CD8+ T cells |
| 5 | Memory CD8+ T cells |
| 6 | Exhausted CD8+ T cells |
| 7 | Progenitor exhausted CD8+ T cells |
| 8 | Terminally exhausted CD8+ T cells |
| 9 | IFN-gamma+ Th1 cells |
| 10 | IL-17A+ Th17 cells |
| 11 | T follicular helper cells |
| 12 | FOXP3+ Tumor-infiltrating Tregs |
| 13 | Cycling T cells |

## Requirements

R >= 4.1.0

Key packages: Seurat (>=4.0), SingleR, celldex, presto, clustree, scCustomize, ggplot2, dplyr, tibble, cowplot, ggrepel, plyr

Install all packages by running `scripts/01_install_and_load_packages.R`.

## How to Reproduce

```r
# Step 1 - Set up directories
source("scripts/00_create_directories.R")

# Step 2 - Install packages
source("scripts/01_install_and_load_packages.R")

# Step 3 - Download data
source("scripts/02_download_data.R")

# Step 4-13 - Run scripts in numbered order
# Each script builds on objects saved by the previous one
```

## Author

Sachi Tiwari
M.Sc. Bioinformatics, Babasaheb Bhimrao Ambedkar University, Lucknow
Email: sachitiwari722@gmail.com
LinkedIn: linkedin.com/in/sachi-tiwari
GitHub: github.com/sachi-tiwari

## License

This project is open-source under the MIT License.
