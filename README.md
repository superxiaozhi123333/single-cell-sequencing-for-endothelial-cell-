## Contents

- [Overview](#overview)
- [Repo Contents](#repo-contents)
- [System Requirements](#system-requirements)
- [Installation Guide](#installation-guide)
- [Demo](#demo)
- [License](./LICENSE)
- [Issues](https://github.com/superxiaozhi123333/single-cell-sequencing-for-endothelial-cell-)

# Overview
Single-cell RNA sequencing (scRNA-seq) generates high-dimensional data with thousands of genes and cells, posing challenges in analysis due to noise, batch effects, and dimensionality. This pipeline (scAnalyze) provides an end-to-end workflow for processing scRNA-seq data, including quality control, normalization, integration, clustering, and differential expression analysis. It leverages state-of-the-art tools (e.g., Seurat) to reduce dimensionality, correct batch effects, and identify cell-type markers, enabling robust biological interpretation in both simulated and real datasets.
# Repo Contents
R/: Core R scripts for data preprocessing, integration, clustering, and visualization (e.g., preprocess.R, integration.R, DE_analysis.R).
data/: Simulated demo dataset (demo_counts.csv.gz, demo_metadata.csv) and example real datasets (e.g., GSE131907 subset).
docs/: Comprehensive documentation, including step-by-step tutorials and parameter explanations.
tests/: Unit tests for key functions (e.g., quality control checks, integration validity) using testthat.
vignettes/: Interactive HTML vignettes for R session help (e.g., scAnalyze_basics.html, advanced_analysis.html).

# System Requirements

## Hardware Requirements
A standard desktop/laptop is sufficient for small datasets (≤5,000 cells). For large datasets (≥100,000 cells), we recommend:
RAM: 32+ GB (for memory-intensive steps like integration)
CPU: 8+ cores, 3.0+ GHz/core (for parallelized operations)
Storage: 100+ GB (for raw data and intermediate outputs)
GPU (optional): NVIDIA GPU with 8+ GB VRAM (accelerates UMAP/TSNE by ~2x)
## Software Requirements

### OS Requirements
Tested on:
Linux: Ubuntu 22.04 LTS
macOS: Sonoma 14.0+
Windows: 11 (64-bit, with WSL2 or Rtools43)
Dependencies

R ≥ v4.3.0 (recommended v4.4.0)
R Packages (versions locked for reproducibility):

Seurat (v4.4.0),
dplyr (v1.2.0),
patchwork (v1.1.2),
BiocManager (v1.30.20),
cowplot (v1.1.1),
ggplot2 (v3.5.2),
stringr (v1.7.0),
paletteer (v1.6.0),
MySeuratWrappers (v0.1.2),
SCpubr (v0.2.6),
viridis (v0.6.3),
ComplexHeatmap (v2.16.2),
limma (v3.54.1),
ggrepel (v0.9.3),
clusterProfiler (v3.16.1),
UCell (v1.1.0),
Cell-ID (v1.21),
Harmony (v1.0)

# Installation Guide

## Stable Release

Install the latest stable version from CRAN:

install.packages("scAnalyze")  
Development Version (GitHub)
Pre-Install Dependencies
First, install required R packages:

install.packages(c("Seurat", "dplyr", "patchwork", "cowplot", "ggplot2", "stringr", "ggrepel"))  
BiocManager::install(c("ComplexHeatmap", "limma", "SCpubr", "viridis", "MySeuratWrappers"))  

Install scAnalyze
Install the development version with vignettes:

require(devtools)  
install_github("yourusername/scAnalyze", build_vignettes = TRUE, force = TRUE)  
require(scAnalyze)  
vignette("scAnalyze_basics", package = "scAnalyze")  # View introductory vignette  

Installation Time
Stable release: ~5 minutes (download + dependencies).
Development version with vignettes: ~15 minutes (on a 16 GB RAM, 8-core machine).

Demo
Run the Demo Pipeline
Step 1: Download Demo Data
The simulated dataset includes 500 cells (250 Normal + 250 Tumor) and 2,000 genes. Download it from:
bash
wget https://example.com/scAnalyze_demo/demo_counts.csv.gz  
wget https://example.com/scAnalyze_demo/demo_metadata.csv  

Step 2: Run Analysis
In R:

# Load package and data  
library(scAnalyze)  
counts <- read.csv("demo_counts.csv.gz", row.names = 1)  
metadata <- read.csv("demo_metadata.csv")  

# Run full pipeline (preprocessing → integration → clustering → DE analysis)  
results <- scAnalyze_pipeline(counts, metadata)  

# View outputs  
print(results$qc_report)  # Quality control summary  
plot(results$umap_plot)   # UMAP clustering plot  
head(results$DE_genes)    # Top differential genes  
Expected Outputs
QC Report: Table of filtered cells, mitochondrial gene stats.
Visualizations: UMAP/TSNE plots (umap.png, tsne.png), volcano plot (volcano.png).
Results Files: Clustering assignments (clusters.csv), differential genes (DE_genes.csv).

Run Time
~20 minutes on a 16 GB RAM, 8-core machine (without GPU acceleration).

MIT License

Copyright (c) 2024 SingleronBio

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.