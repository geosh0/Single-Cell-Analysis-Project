import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from typing import List, Tuple, Optional, Union

def clean_and_validate_df(df: pd.DataFrame, index_name: str = 'GENE') -> pd.DataFrame:
    """
    Cleans a gene expression DataFrame by ensuring numeric types, removing artifacts,
    handling duplicates, and filling missing values.

    Args:
        df (pd.DataFrame): Raw gene expression dataframe.
        index_name (str): The name of the index column to check for artifact headers.

    Returns:
        pd.DataFrame: A cleaned, float-type dataframe with unique indices.
    """
    # Create a copy to avoid SettingWithCopy warnings
    df_clean = df.copy()

    # 1. Remove rows that might be leftover headers (e.g., index matches the column name)
    # Using specific check to avoid errors if index is not unique yet
    if index_name in df_clean.index:
        print(f" > Dropping artifact row matching index name: '{index_name}'")
        df_clean = df_clean.drop(index_name, axis=0)

    # 2. Convert to float
    # We use apply(pd.to_numeric) to handle entire columns at once safely
    try:
        df_clean = df_clean.astype(float)
    except ValueError:
        print(" > Warning: Could not strictly convert to float. Coercing non-numeric values to NaN.")
        df_clean = df_clean.apply(pd.to_numeric, errors='coerce')

    # 3. Handle duplicates (Summing fragments of same gene)
    if df_clean.index.has_duplicates:
        n_dupes = df_clean.index.duplicated().sum()
        print(f" > Handling {n_dupes} duplicate gene entries by summing.")
        df_clean = df_clean.groupby(df_clean.index).sum()

    # 4. Fill NaNs
    if df_clean.isnull().values.any():
        print(" > Filling NaN values with 0.")
        df_clean = df_clean.fillna(0)
    
    return df_clean

def merge_expression_tables(dfs: List[pd.DataFrame]) -> pd.DataFrame:
    """
    Merges a list of expression tables along the columns (samples).
    
    Args:
        dfs (List[pd.DataFrame]): A list of dataframes to merge.
        
    Returns:
        pd.DataFrame: A single merged and cleaned dataframe.
    """
    if not dfs:
        raise ValueError("Input list of DataFrames is empty.")
    
    print(f"Merging {len(dfs)} dataframes...")
    
    # Concatenate along columns (axis=1)
    # join='outer' ensures we keep all genes, filling missing ones with NaN (which acts as 0 later)
    combined = pd.concat(dfs, axis=1, join='outer')
    
    # Run standard cleaning on the combined result
    combined_clean = clean_and_validate_df(combined)
    
    print(f"Merge successful. Final Shape: {combined_clean.shape}")
    return combined_clean

def plot_qc_metrics(df: pd.DataFrame, 
                    detection_threshold: float = 0.0, 
                    figsize: Tuple[int, int] = (18, 5)) -> Tuple[pd.Series, pd.Series]:
    """
    Generates standard QC plots for Single Cell expression data.
    
    Args:
        df (pd.DataFrame): Genes x Samples DataFrame (TPM or Counts).
        detection_threshold (float): Cutoff to consider a gene 'detected' (usually 0 or 1).
        figsize (Tuple[int, int]): Tuple determining plot size.
        
    Returns:
        Tuple[pd.Series, pd.Series]: (Library Sizes, Detected Gene Counts)
    """
    # Set a clean style for the plots
    sns.set_theme(style="whitegrid")
    
    fig, axes = plt.subplots(1, 3, figsize=figsize)
    
    # 1. Library Size (Total TPM/Counts per sample)
    lib_size = df.sum(axis=0)
    sns.barplot(x=lib_size.index, y=lib_size.values, ax=axes[0], color='#4c72b0')
    axes[0].set_title('Total Expression per Sample (Library Size)')
    axes[0].set_ylabel('Total TPM/Counts')
    axes[0].tick_params(axis='x', rotation=90, labelsize=8)

    # 2. Detected Genes per Sample
    n_genes = (df > detection_threshold).sum(axis=0)
    sns.barplot(x=n_genes.index, y=n_genes.values, ax=axes[1], color='#55a868')
    axes[1].set_title(f'Detected Genes per Sample (>{detection_threshold})')
    axes[1].set_ylabel('Count of Genes')
    axes[1].tick_params(axis='x', rotation=90, labelsize=8)

    # 3. Gene Detection Frequency (How many samples have a specific gene?)
    gene_freq = (df > detection_threshold).sum(axis=1)
    sns.histplot(gene_freq, bins=30, kde=False, ax=axes[2], color='#8172b3')
    axes[2].set_title('Gene Detection Frequency')
    axes[2].set_xlabel('Number of Samples Gene is Detected In')
    axes[2].set_ylabel('Number of Genes')
    
    plt.tight_layout()
    plt.show()

    return lib_size, n_genes