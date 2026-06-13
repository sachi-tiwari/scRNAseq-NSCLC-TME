# =============================================================================
# Script: 02_download_data.R
# Description: Downloads the raw NSCLC scRNA-seq count matrix and metadata
#              from Figshare (public repository).
# Note: The count matrix is a large file (~1GB). If the download times out,
#       increase the timeout limit or download manually.
#       Count matrix: https://figshare.com/ndownloader/files/39537094
#       Metadata:     https://figshare.com/ndownloader/files/39537250
# Input: None
# Output: raw_data/RNA_rawcounts_matrix.rds
#         raw_data/metadata.csv
# =============================================================================

# Increase timeout limit for large file download (in seconds)
options(timeout = 300)

# Download raw count matrix
download.file(
  url = "https://figshare.com/ndownloader/files/39537094",
  destfile = "raw_data/RNA_rawcounts_matrix.rds",
  mode = "wb"
)

# Download metadata
download.file(
  url = "https://figshare.com/ndownloader/files/39537250",
  destfile = "raw_data/metadata.csv",
  mode = "wb"
)

cat("Downloads complete:\n")
cat("  raw_data/RNA_rawcounts_matrix.rds\n")
cat("  raw_data/metadata.csv\n")
cat("If any download failed, use the URLs in the header above.\n")