# Shared Utilities: The Core Engine

This package contains the universal logic used across all seven datasets. By centralizing the mathematical operations and benchmarking steps, this project adheres to **DRY (Don't Repeat Yourself)** software engineering principles.

## 🛠️ Modules

### 1. `sc_processor.py` (The Mathematical Core)
This module handles the transformation of raw gene expression matrices into analyzable principal components. It ensures that every dataset, regardless of source format, undergoes the same rigorous statistical treatment.

**Key Functions:**
*   `filter_tpm_matrix`: Removes low-quality cells (based on gene count) and low-expression genes.
*   `log_transform`: Applies `log1p` normalization to stabilize variance.
*   `select_highly_variable_genes`: Identifies the top $N$ genes driving biological heterogeneity.
*   `scale_data`: Performs Z-score scaling (StandardScaler) for PCA readiness.
*   `run_pca_pipeline`: Executes dimensionality reduction and returns variance ratios.

### 2. `sc_clustering.py` (The Benchmarking Engine)
This module provides a framework for objectively evaluating unsupervised learning results. Instead of arbitrarily choosing a clustering algorithm, this tool runs multiple algorithms and scores them against known biological ground truths.

**Key Features:**
*   **Multi-Target Benchmarking:** Can test clusters against a list of metadata columns (e.g., `['TimePoint', 'Genotype']`) simultaneously to determine the primary driver of variance.
*   **Algorithms:** K-Means, Agglomerative (Hierarchical), and Spectral Clustering.
*   **Metrics:** Uses **Adjusted Mutual Information (AMI)** and **Adjusted Rand Index (ARI)** to quantify cluster accuracy.
