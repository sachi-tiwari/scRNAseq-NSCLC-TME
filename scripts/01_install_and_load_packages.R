# =============================================================================
# Script: 01_install_and_load_packages.R
# Description: Installs and loads all R packages required for the NSCLC
#              scRNA-seq analysis pipeline.
# Note: Run this script once before starting the analysis.
#       Installation lines can be commented out after first run.
# Input: None
# Output: None (packages loaded into R environment)
# =============================================================================

# --- Install CRAN packages ---
install.packages(c("Seurat", "tidyverse", "patchwork", "Matrix", "future",
                   "ggplot2", "dplyr", "tibble", "cowplot", "ggrepel",
                   "clustree", "scCustomize"))

# --- Install Bioconductor packages ---
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("SingleR", "celldex"))

# --- Install presto (fast Wilcoxon test for large datasets) ---
if (!requireNamespace("devtools", quietly = TRUE))
  install.packages("devtools")

devtools::install_github("immunogenomics/presto")

# --- Load all packages ---
library(Seurat)
library(tidyverse)
library(patchwork)
library(Matrix)
library(future)
library(ggplot2)
library(dplyr)
library(tibble)
library(cowplot)
library(ggrepel)
library(clustree)
library(scCustomize)
library(SingleR)
library(celldex)
library(presto)