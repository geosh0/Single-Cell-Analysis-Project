# Modular Single-Cell RNA-Seq Analysis Pipeline
## 🧬 Project Overview
This repository hosts a portfolio of Single-Cell RNA-Seq (scRNA-seq) analyses performed on seven distinct public datasets.
The primary goal of this project was not just to analyze biological data, but to engineer a scalable, reproducible software architecture. Instead of repeating code across notebook files, this project refactors common logic into a shared utility package, adhering to DRY (Don't Repeat Yourself) principles. This separation of concerns allows for dataset-specific nuances (like file parsing) to be handled locally, while mathematical operations (normalization, dimensionality reduction) remain universal.

## 🏗️ Software Architecture
The project is organized into two distinct layers: Local Adapters and Shared Logic.

1. The Shared Core (shared_utils/)
This folder contains the "product" logic—code that is agnostic to the data source and can be reused across any single-cell project.
* sc_processor.py: The mathematical engine. It handles Quality Control (QC), Log-Normalization, Highly Variable Gene (HVG) selection, Scaling, and PCA.
* sc_clustering.py: The evaluation engine. It runs multiple clustering algorithms (KMeans, Spectral, Agglomerative) and benchmarks them against biological ground truths (e.g., Cell Type, Timepoint) using Adjusted Mutual Information (AMI) and Adjusted Rand Index (ARI).

2. The Local Implementations (GSE_XXXXX/)
Each dataset has its own dedicated folder containing:
* sc_loader.py: A custom ETL script designed to handle the specific messy reality of that dataset (e.g., parsing raw Series Matrix text files, aggregating SRA runs, or pivoting FPKM tracking files).
* notebook.ipynb: The narrative interface. It imports the local loader to get the data, then calls the shared core to process it.


## 🚀 Key Features 
* Data Ingestion: Handled 7 different file formats, including raw text counts, compressed .gz archives, and complex SRA metadata tables.
* Automated QC: Standardized filtering pipelines based on library size and gene detection rates.
* Multi-Target Benchmarking: The clustering engine can test algorithms against multiple biological variables simultaneously (e.g., comparing if clusters align better with Genotype vs. Developmental Stage).
* Reproducibility: All random states and parameter configurations are centralized, ensuring results can be replicated.

## 📊 Methodology
For every dataset, the following standardized pipeline is applied:
1. ETL: Data is loaded and aligned with SRA Metadata using the local sc_loader.
2. Filtering: Low-quality cells and genes are removed based on dynamic thresholds (e.g., min genes > 2000).
3. Normalization: Library size normalization followed by Log1p transformation.
4. Feature Selection: Selection of the top 2,000-3,000 Highly Variable Genes (HVGs).
5. Dimensionality Reduction: Z-Score scaling followed by PCA (Principal Component Analysis).
6. Unsupervised Learning: Clustering is performed over a range of k values and scored against metadata labels using AMI/ARI metrics.

## 🔗 Data Sources
All data used in this project is publicly available via the NCBI Gene Expression Omnibus (GEO).
>  Data Repository:  [GEO-NCBI](https://www.ncbi.nlm.nih.gov/geo/)

---
## ⚠️ NOTE

This project was developed by a graduate student as a part of the bachelors thesis to cluster single cell datasets. While the mathematical pipelines utilize standard libraries (Scikit-Learn, Pandas, NumPy) and follow best practices for single-cell analysis, this code is intended for educational and demonstration purposes rather than production-grade clinical use. Feedback and contributions are welcome.
