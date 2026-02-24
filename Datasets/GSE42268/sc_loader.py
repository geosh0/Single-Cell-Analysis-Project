import pandas as pd
import os
import glob
from typing import Tuple, List, Dict

def load_gse42268_data(
    data_folder: str, 
    metadata_path: str
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Loads individual TXT expression files, pivots them into a matrix,
    and merges sample information with SRA metadata.

    Args:
        data_folder (str): Path to folder containing .txt files.
        metadata_path (str): Path to SRA Run Table CSV.

    Returns:
        Tuple[pd.DataFrame, pd.DataFrame]: 
            1. Expression Matrix (Index=Genes, Columns=Sample_IDs)
            2. Metadata DataFrame (Index=Sample_IDs)
    """
    # --- 1. Load Files & Extract Info ---
    files = glob.glob(os.path.join(data_folder, '*.txt'))
    print(f"Found {len(files)} files to load in: {data_folder}")
    
    if not files:
        raise ValueError("No .txt files found. Check your DATA_PATH.")

    expression_dict = {}
    sample_info_list = []
    
    for filename in files:
        try:
            # Parse filename for ID and Group (Your specific logic)
            base = os.path.basename(filename)
            # Example filename: "GSM1035628_ES_1.txt" (Hypothetical)
            sample_id = base.split('_')[0]
            
            # Logic to extract 'group' from filename
            group_raw = base.replace(sample_id + '_', '').split('.')[0]
            group = group_raw.split('_')[0] if '_' in group_raw else group_raw
            
            # Read Data
            # Assumes 3 columns: id, gene, fpkm
            df_temp = pd.read_csv(filename, sep="\t", header=0)
            
            # We assume column 1 is Gene Symbol and column 2 is FPKM
            # (Adjust indices if your txt files are different)
            df_temp.columns = ['id', 'gene_symbol', 'fpkm']
            
            # Handle duplicate genes within one file (summing FPKM)
            if df_temp['gene_symbol'].duplicated().any():
                df_temp = df_temp.groupby('gene_symbol')['fpkm'].sum().reset_index()
            
            # Store expression vector
            expression_dict[sample_id] = df_temp.set_index('gene_symbol')['fpkm']
            
            # Store metadata derived from filename
            sample_info_list.append({
                'Sample_ID': sample_id,
                'Filename_Group': group
            })
            
        except Exception as e:
            print(f"Error loading {filename}: {e}")

    # --- 2. Create Expression Matrix ---
    # Combine all vectors into one dataframe (Genes x Samples)
    print("Constructing expression matrix...")
    expression_matrix = pd.DataFrame(expression_dict)
    
    # Fill NaNs with 0 (genes present in some samples but not others)
    expression_matrix = expression_matrix.fillna(0)
    
    print(f"Expression Matrix Shape: {expression_matrix.shape}")

    # --- 3. Process Metadata ---
    # Create dataframe from filename info
    local_meta = pd.DataFrame(sample_info_list)
    
    # Load SRA Metadata
    if os.path.exists(metadata_path):
        sra_df = pd.read_csv(metadata_path)
        
        # Ensure join keys are strings
        sra_df['GEO_Accession (exp)'] = sra_df['GEO_Accession (exp)'].astype(str)
        local_meta['Sample_ID'] = local_meta['Sample_ID'].astype(str)
        
        # Merge
        # We merge 'left' to keep all samples we actually loaded data for
        metadata_merged = pd.merge(
            local_meta,
            sra_df,
            left_on='Sample_ID',
            right_on='GEO_Accession (exp)',
            how='left'
        )
        
        # Set Index to match Expression Matrix columns
        metadata_merged.set_index('Sample_ID', inplace=True)
        
        # Clean up
        if 'GEO_Accession (exp)' in metadata_merged.columns:
            metadata_merged.drop(columns=['GEO_Accession (exp)'], inplace=True)
            
        print(f"Metadata merged. Shape: {metadata_merged.shape}")
        return expression_matrix, metadata_merged
    
    else:
        print("Warning: Metadata path not found. Returning filename info only.")
        local_meta.set_index('Sample_ID', inplace=True)
        return expression_matrix, local_meta