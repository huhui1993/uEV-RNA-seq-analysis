#!/usr/bin/env Rscript
# =============================================================================
# Install Required R Packages for Urine Exosome Transcriptome Analysis
# Run this script first before any analysis scripts
# Usage: Rscript install_packages.R
# =============================================================================

cat("=== Installing required R packages ===\n\n")

# Vector of CRAN packages
cran_packages <- c(
  # Data manipulation
  "reshape", "dplyr", "tidyr",
  # Visualization
  "ggplot2", "ggpubr", "ggsci", "pheatmap", "RColorBrewer",
  "scales", "ComplexHeatmap",
  # Statistics
  "genefilter",
  # PCA
  "FactoMineR", "factoextra",
  # Time series
  "Mfuzz", "Biobase",
  # Venn diagrams
  "venn", "VennDiagram",
  # Miscellaneous
  "stringr", "corrplot"
)

# Vector of Bioconductor packages
bioc_packages <- c(
  "DESeq2",          # Differential expression
  "limma",           # Linear models + avereps
  "clusterProfiler", # GO/KEGG enrichment
  "org.Hs.eg.db",    # Human genome annotation
  "GSVA",            # Gene set variation analysis
  "GSEABase",        # Gene set enrichment
  "edgeR",           # RNA-seq analysis (CPM)
  "sva",             # Surrogate variable analysis
  "scran",           # Single-cell analysis utilities
  "scater"           # Single-cell QC
)

# Install CRAN packages
cat("1. Installing CRAN packages...\n")
install_if_missing <- function(pkg) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
    cat(sprintf("   Installed: %s\n", pkg))
  } else {
    cat(sprintf("   Already installed: %s\n", pkg))
  }
}
invisible(lapply(cran_packages, install_if_missing))

# Install Bioconductor packages
cat("\n2. Installing Bioconductor packages...\n")
if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}
invisible(lapply(bioc_packages, function(pkg) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    BiocManager::install(pkg, ask = FALSE)
    cat(sprintf("   Installed: %s\n", pkg))
  } else {
    cat(sprintf("   Already installed: %s\n", pkg))
  }
}))

# Install msigdbr (MSigDB gene sets)
cat("\n3. Installing msigdbr...\n")
if (!require("msigdbr", quietly = TRUE)) {
  install.packages("msigdbr", repos = "https://cloud.r-project.org")
}

# Check FactoInvestigate (may need additional dependencies)
cat("\n4. Checking FactoInvestigate...\n")
if (!require("FactoInvestigate", quietly = TRUE)) {
  cat("   FactoInvestigate not installed (optional, used for PCA report generation)\n")
}

cat("\n=== All packages installed successfully! ===\n")
