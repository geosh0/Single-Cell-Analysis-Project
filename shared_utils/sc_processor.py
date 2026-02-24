import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from typing import Tuple, Optional, List
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA

# ==========================================
# CORE PROCESSING FUNCTIONS
# ==========================================

def filter_tpm_matrix(
    df: pd.DataFrame, 
    min_tpm: float = 1.0, 
    min_genes_per_sample: int = 500, 
    min_samples_per_gene: int = 3
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Filters a TPM expression matrix based on detection thresholds.
    
    Args:
        df: Genes x Samples DataFrame.
        min_tpm: Threshold to consider a gene 'detected'.
        min_genes_per_sample: Minimum detected genes required to keep a sample.
        min_samples_per_gene: Minimum samples required to keep a gene.

    Returns:
        Tuple containing:
        1. Filtered DataFrame.
        2. QC Metrics DataFrame (metadata for the kept samples).
    """
    if df.empty:
        print("Warning: Input DataFrame is empty.")
        return pd.DataFrame(), pd.DataFrame()

    print(f"--- QC Filtering ---")
    print(f"Initial shape: {df.shape}")

    # 1. Sample QC
    # Calculate how many genes exceed the TPM threshold per sample
    genes_detected = (df > min_tpm).sum(axis=0)
    samples_mask = genes_detected >= min_genes_per_sample
    df_filtered_samples = df.loc[:, samples_mask]
    
    print(f"Samples kept: {df_filtered_samples.shape[1]} (Threshold: {min_genes_per_sample} genes > {min_tpm} TPM)")

    if df_filtered_samples.empty:
        print("No samples passed QC filtering.")
        return pd.DataFrame(), pd.DataFrame()

    # 2. Gene QC
    # Calculate how many samples have this gene expressed (using the filtered samples)
    samples_per_gene = (df_filtered_samples > min_tpm).sum(axis=1)
    genes_mask = samples_per_gene >= min_samples_per_gene
    df_final = df_filtered_samples.loc[genes_mask, :]
    
    print(f"Genes kept: {df_final.shape[0]} (Threshold: detected in {min_samples_per_gene} samples)")

    # 3. Generate Metrics for downstream analysis
    qc_metrics = pd.DataFrame(index=df_final.columns)
    qc_metrics['n_genes_detected'] = genes_detected[df_final.columns]
    qc_metrics['total_tpm'] = df_final.sum(axis=0)

    return df_final, qc_metrics

def log_transform(df: pd.DataFrame, method: str = 'log1p') -> pd.DataFrame:
    """
    Applies log transformation to the expression matrix.
    """
    if df.empty:
        return df

    if method == 'log2':
        return np.log2(df + 1)
    else:
        # Default to natural log (ln(x+1))
        return np.log1p(df)

def select_highly_variable_genes(
    df: pd.DataFrame, 
    n_top_genes: int = 2000
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Selects highly variable genes (HVGs) based on dispersion (variance/mean).
    Assumes df is log-transformed.
    """
    if df.empty:
        return pd.DataFrame(), pd.DataFrame()

    # Calculate statistics
    mean_expr = df.mean(axis=1)
    var_expr = df.var(axis=1)
    
    # Calculate dispersion (Var / Mean)
    # Adding epsilon to avoid division by zero
    epsilon = 1e-9
    dispersion = var_expr / (mean_expr + epsilon)
    dispersion = dispersion.fillna(0)
    
    # Select top N genes
    n_top_genes = min(n_top_genes, len(dispersion))
    hvg_indices = dispersion.nlargest(n_top_genes).index
    
    df_hvg = df.loc[hvg_indices]
    
    # Pack metrics for plotting
    metrics = pd.DataFrame({
        'mean': mean_expr,
        'variance': var_expr,
        'dispersion': dispersion,
        'is_hvg': dispersion.index.isin(hvg_indices)
    })
    
    print(f"HVG Selection: {df_hvg.shape[0]} genes selected.")
    return df_hvg, metrics

def scale_data(df: pd.DataFrame) -> pd.DataFrame:
    """
    Scales data to zero mean and unit variance (Z-score).
    Operates Gene-wise (rows are genes, columns are samples).
    """
    if df.empty:
        return pd.DataFrame()

    # StandardScaler expects (Samples x Features). 
    # Our df is (Genes x Samples). 
    # So we Transpose -> Scale -> Transpose back.
    scaler = StandardScaler()
    
    scaled_values = scaler.fit_transform(df.T).T
    
    df_scaled = pd.DataFrame(
        scaled_values, 
        index=df.index, 
        columns=df.columns
    )
    return df_scaled

def run_pca_pipeline(
    df_scaled: pd.DataFrame, 
    n_components: int = 50
) -> Tuple[pd.DataFrame, np.ndarray, PCA]:
    """
    Runs PCA on scaled data.
    Returns: PC Coordinates, Explained Variance Ratio, PCA Object
    """
    if df_scaled.empty:
        return pd.DataFrame(), np.array([]), None
        
    # PCA expects (Samples x Features)
    X = df_scaled.T 
    
    # Handle n_components validation
    n_samples, n_features = X.shape
    actual_n = min(n_components, n_samples, n_features)
    
    pca = PCA(n_components=actual_n, random_state=42)
    X_pca = pca.fit_transform(X)
    
    # Create DataFrame for results
    pc_cols = [f'PC{i+1}' for i in range(actual_n)]
    df_pca = pd.DataFrame(X_pca, index=X.index, columns=pc_cols)
    
    return df_pca, pca.explained_variance_ratio_, pca

# ==========================================
# PLOTTING FUNCTIONS
# ==========================================

def plot_expression_distribution(df: pd.DataFrame, title: str = "Mean Log-Expression Distribution"):
    """Visualizes the distribution of mean expression per gene."""
    if df.empty: return
    mean_expr = df.mean(axis=1)
    plt.figure(figsize=(10, 6))
    sns.histplot(mean_expr, bins=50, kde=True, color='teal')
    plt.title(title)
    plt.xlabel('Mean Expression value')
    plt.ylabel('Count of Genes')
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    plt.show()

def plot_hvg_dispersion(metrics: pd.DataFrame):
    """Plots Mean vs Variance, highlighting HVGs."""
    if metrics.empty: return
    plt.figure(figsize=(8, 6))
    # Plot all genes
    sns.scatterplot(
        data=metrics[~metrics['is_hvg']], x='mean', y='variance', 
        color='gray', alpha=0.3, s=10, label='Other Genes'
    )
    # Plot HVGs
    sns.scatterplot(
        data=metrics[metrics['is_hvg']], x='mean', y='variance', 
        color='red', alpha=0.6, s=15, label='HVGs'
    )
    plt.xscale('log'); plt.yscale('log')
    plt.title('Gene Dispersion (Mean vs Variance)')
    plt.xlabel('Mean Expression'); plt.ylabel('Variance')
    plt.legend(); plt.show()

def plot_pca_results(df_pca: pd.DataFrame, variance_ratio: np.ndarray):
    """Plots Scree Plot and PC1 vs PC2."""
    if df_pca.empty: return
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))
    
    # 1. Scree Plot
    n_pcs = len(variance_ratio)
    axes[0].plot(range(1, n_pcs+1), variance_ratio, 'o-', color='steelblue')
    axes[0].set_title('Scree Plot (Variance Explained)')
    axes[0].set_xlabel('Principal Component')
    axes[0].set_ylabel('Variance Ratio')
    
    # 2. PC1 vs PC2
    if n_pcs >= 2:
        sns.scatterplot(
            data=df_pca, x='PC1', y='PC2', 
            s=100, color='purple', ax=axes[1]
        )
        for idx, row in df_pca.iterrows():
            axes[1].text(row['PC1']+.02, row['PC2']+.02, idx, fontsize=8)
            
        axes[1].set_xlabel(f'PC1 ({variance_ratio[0]*100:.1f}%)')
        axes[1].set_ylabel(f'PC2 ({variance_ratio[1]*100:.1f}%)')
        axes[1].set_title('PCA: PC1 vs PC2')
        axes[1].grid(True, linestyle='--', alpha=0.5)
        
    plt.tight_layout(); plt.show()