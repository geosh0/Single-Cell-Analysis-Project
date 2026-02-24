import pandas as pd
import os
import glob
from typing import Tuple

def load_dataset_5(
    expression_path: str, 
    metadata_path: str
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Loads .fpkm_tracking files, pivots them into a matrix, 
    and filters for Single-Cell samples only.
    """
    
    # --- 1. Load SRA Metadata ---
    try:
        sra_df = pd.read_csv(metadata_path)
        # Ensure ID column is string for merging
        if 'GEO_Accession (exp)' in sra_df.columns:
            sra_df['GEO_Accession (exp)'] = sra_df['GEO_Accession (exp)'].astype(str)
    except Exception as e:
        raise ValueError(f"Could not load metadata: {e}")

    # --- 2. Load Expression Files ---
    all_data = []
    files = glob.glob(os.path.join(expression_path, '*.fpkm_tracking'))
    print(f"Found {len(files)} files. Loading...")
    
    for filename in files:
        try:
            df = pd.read_csv(filename, sep="\t")
            
            # Extract Sample ID (GSM...) from filename
            sample_id = os.path.basename(filename).split('_')[0]
            
            # Keep only necessary columns to save memory before merging
            # We assume columns are 'tracking_id' (gene) and 'FPKM'
            if 'tracking_id' in df.columns and 'FPKM' in df.columns:
                subset = df[['tracking_id', 'FPKM']].copy()
                subset['sample_Id'] = sample_id
                all_data.append(subset)
                
        except Exception as e:
            print(f"Error loading {filename}: {e}")

    if not all_data:
        raise ValueError("No expression data loaded.")

    # Combine into one long dataframe
    long_df = pd.concat(all_data, ignore_index=True)
    
    # --- 3. Merge Metadata ---
    # We merge *before* pivoting to filter for Single Cells efficiently
    merged_df = pd.merge(
        long_df,
        sra_df[['GEO_Accession (exp)', 'cell_type', 'technology_used']],
        left_on='sample_Id',
        right_on='GEO_Accession (exp)',
        how='left'
    )
    
    # --- 4. Filter for Single Cell Only ---
    # Your specific logic: keep only 'single-cell RNA-seq'
    if 'technology_used' in merged_df.columns:
        sc_subset = merged_df[merged_df['technology_used'] == 'single-cell RNA-seq']
        print(f"Filtering: Kept {sc_subset['sample_Id'].nunique()} Single-Cell samples (dropped Bulk).")
    else:
        print("Warning: 'technology_used' column not found. Keeping all data.")
        sc_subset = merged_df

    # --- 5. Pivot to Matrix (Genes x Samples) ---
    print("Pivoting to create Expression Matrix...")
    expression_matrix = sc_subset.pivot_table(
        index='tracking_id', 
        columns='sample_Id', 
        values='FPKM', 
        fill_value=0
    )
    
    # --- 6. Prepare Final Metadata ---
    # Create a metadata table indexed by sample_Id
    metadata_final = sc_subset[['sample_Id', 'cell_type', 'technology_used']].drop_duplicates().set_index('sample_Id')
    
    # Align
    metadata_final = metadata_final.loc[expression_matrix.columns]
    
    print(f"Final Matrix Shape: {expression_matrix.shape}")
    return expression_matrix, metadata_final