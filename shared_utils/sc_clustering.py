import pandas as pd
import numpy as np
from typing import List, Dict, Union, Optional
from sklearn.cluster import KMeans, AgglomerativeClustering, SpectralClustering
from sklearn.metrics import adjusted_rand_score, adjusted_mutual_info_score

def run_clustering_benchmark(
    pca_df: pd.DataFrame, 
    cell_metadata: pd.DataFrame, 
    n_clusters_range: List[int],
    target_cols: List[str]
) -> pd.DataFrame:
    """
    Benchmarks Baseline Clustering Models (KMeans, HClust, Spectral) against 
    ONE OR MORE ground truth columns.
    
    Args:
        pca_df: PCA coordinates (Samples x PCs).
        cell_metadata: Metadata containing ground truth labels.
        n_clusters_range: List of k values to test.
        target_cols: List of column names in metadata to test against (e.g., ['Stage', 'Cross']).
        
    Returns:
        pd.DataFrame: Leaderboard of AMI/ARI scores for each method/k/target.
    """
    print("\n--- Starting Multi-Target Clustering Benchmark ---")
    
    all_results = [] 
    
    # --- 1. Alignment & Validation ---
    # Intersection prevents errors if indices are slightly different
    common_index = pca_df.index.intersection(cell_metadata.index)
    
    if len(common_index) < len(pca_df):
        print(f"Warning: Dropping {len(pca_df) - len(common_index)} samples not found in metadata.")
        
    # Align dataframes
    pca_aligned = pca_df.loc[common_index]
    meta_aligned = cell_metadata.loc[common_index]
    
    # Validate Targets
    valid_targets = [col for col in target_cols if col in meta_aligned.columns]
    if not valid_targets:
        print(f"Error: None of the target columns {target_cols} found in metadata.")
        return pd.DataFrame()
    
    print(f"Benchmarking against: {valid_targets}")

    # --- 2. Define Models (User-Specific Configuration) ---
    models = {
        'KMeans': lambda k: KMeans(
            n_clusters=k, n_init=10, init='random', 
            algorithm='lloyd', max_iter=10, random_state=42
        ),
        'HClust': lambda k: AgglomerativeClustering(
            n_clusters=k, linkage='single', metric='euclidean'
        ), 
        'Spectral': lambda k: SpectralClustering(
            n_clusters=k, affinity='nearest_neighbors', 
            assign_labels='kmeans', random_state=42
        )
    }
    
    # --- 3. Execution Loop ---
    for k in n_clusters_range:
        print(f"\nProcessing k = {k}...", end=" ")
        
        for method_name, model_func in models.items():
            try:
                # Instantiate & Fit
                model = model_func(k)
                X = pca_aligned.values
                
                if hasattr(model, 'fit_predict'):
                    labels_pred = model.fit_predict(X)
                else:
                    model.fit(X)
                    labels_pred = model.labels_
                
                # --- Score against ALL targets ---
                row = {'Method': method_name, 'k': k}
                
                for target in valid_targets:
                    y_true = meta_aligned[target].astype(str).values
                    
                    ami = adjusted_mutual_info_score(y_true, labels_pred)
                    ari = adjusted_rand_score(y_true, labels_pred)
                    
                    # Store with dynamic column names
                    row[f'AMI_{target}'] = round(ami, 3)
                    row[f'ARI_{target}'] = round(ari, 3)
                
                all_results.append(row)
                
            except Exception as e:
                print(f"[{method_name} Failed: {e}]")

    print("\nDone.")

    # --- 4. Format Output ---
    results_df = pd.DataFrame(all_results)
    
    if not results_df.empty:
        # Sort by the AMI of the first target provided (primary objective)
        primary_sort = f'AMI_{valid_targets[0]}'
        if primary_sort in results_df.columns:
            results_df = results_df.sort_values(by=primary_sort, ascending=False)
        
    return results_df