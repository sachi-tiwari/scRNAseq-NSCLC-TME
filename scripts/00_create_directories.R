# =============================================================================
# Script: 00_create_directories.R
# Description: Creates the project folder structure for the NSCLC scRNA-seq
#              analysis pipeline. Run this script first before anything else.
# Input: None
# Output: Creates folders — raw_data, data/metadata, results, figures, docs
# =============================================================================

# Set working directory to the project root
# Change this path to where you cloned/downloaded this repository
project_dir <- getwd()  # Uses current directory by default
setwd(project_dir)
cat("Working directory set to:", project_dir, "\n")

# Create subdirectories if they don't already exist
dir.create("raw_data", showWarnings = FALSE)
dir.create("data/metadata", showWarnings = FALSE, recursive = TRUE)
dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)
dir.create("docs", showWarnings = FALSE)

cat("Project directories created (or already exist):\n")
cat(list.files(project_dir), sep = "\n")