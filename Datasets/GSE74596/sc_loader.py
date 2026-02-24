import pandas as pd
import os
import glob
import re
import numpy as np
from typing import Tuple

def _process_sra_metadata(sra_path: str) -> pd.DataFrame:
    """Internal helper to parse and aggregate SRA metadata."""
    print(f"Loading SRA metadata from: {sra_path}")
    try:
        sra_df = pd.read_csv(sra_path)
    except Exception as e:
        print(f"Error loading SRA: {e}")
        return pd.DataFrame()

    # Extract NKT Subset (NKT0, NKT1, etc.)
    # We use source_name as the primary source
    if 'source_name' in sra_df.columns:
        sra_df['NKT_Subset'] = sra_df['source_name'].str.extract(r'(NKT[0-9]+)')
        sra_df['NKT_Subset'] = sra_df['NKT_Subset'].fillna('Unknown')
    else:
        print("Warning: 'source_name' column missing. Cannot parse NKT subsets.")
        sra_df['NKT_Subset'] = 'Unknown'

    # Aggregation Logic (Handling multiple runs per cell)
    # We group by GSM ID ('Sample Name')
    print("Aggregating SRA metadata per cell (GSM ID)...")
    
    agg_funcs = {
        'NKT_Subset': 'first',      # Should be constant for the sample
        'BioSample': 'first',
        'tissue': 'first' if 'tissue' in sra_df.columns else lambda x: None,
        'strain': 'first' if 'strain' in sra_df.columns else lambda x: None
    }
    
    # Only use columns that actually exist
    valid_funcs = {k: v for k, v in agg_funcs.items() if k in sra_df.columns}
    
    sra_agg = sra_df.groupby('Sample Name').agg(valid_funcs)
    sra_agg.index.name = 'GSM_ID'
    
    print(f"Metadata Aggregated. Found {sra_agg.shape[0]} unique cells.")
    return sra_agg

def _load_counts(data_dir: str) -> pd.DataFrame:
    """Internal helper to load individual count files."""
    print(f"Scanning for count files in: {data_dir}")
    
    # Pattern specific to GSE74596
    files = glob.glob(os.path.join(data_dir, 'GSM*_MappedReads_Annotations.txt*'))
    print(f"Found {len(files)} files.")
    
    if not files:
        return pd.DataFrame()

    series_list = []
    
    for fpath in files:
        fname = os.path.basename(fpath)
        # Extract GSM ID (GSM12345_...)
        match = re.match(r'(GSM[0-9]+)_', fname)
        if not match: continue
        gsm_id = match.group(1)
        
        try:
            # Load (2 columns: Gene, Count)
            # No header usually implies column 0=Gene, 1=Count
            df = pd.read_csv(fpath, sep='\t', header=None, names=['Gene', 'Count'])
            
            # Numeric validation
            df['Count'] = pd.to_numeric(df['Count'], errors='coerce').fillna(0)
            
            # Remove duplicates (take first)
            if df['Gene'].duplicated().any():
                df = df.drop_duplicates(subset=['Gene'], keep='first')
            
            # Create Series
            s = df.set_index('Gene')['Count']
            s.name = gsm_id
            series_list.append(s)
            
        except Exception as e:
            print(f"Error loading {fname}: {e}")

    if not series_list:
        return pd.DataFrame()

    print("Concatenating matrix...")
    # Outer join to include all genes found in any file
    # This creates Genes x Samples
    matrix = pd.concat(series_list, axis=1, join='outer').fillna(0)
    
    return matrix

def load_gse74596_data(
    data_dir: str, 
    sra_path: str
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Master loader for GSE74596.
    Returns: (Expression Matrix [Genes x Samples], Metadata)
    """
    # 1. Load Metadata
    meta_df = _process_sra_metadata(sra_path)
    
    # 2. Load Counts
    raw_counts = _load_counts(data_dir)
    
    if raw_counts.empty:
        raise ValueError("No count data loaded.")
        
    # 3. Align
    common = raw_counts.columns.intersection(meta_df.index)
    
    if len(common) == 0:
        print("Error: No intersection between File GSM IDs and SRA GSM IDs.")
        return raw_counts, meta_df # Return raw to debug
        
    print(f"Aligning: Keeping {len(common)} cells matched in both.")
    
    raw_counts = raw_counts[common]
    meta_df = meta_df.loc[common]
    
    print(f"Final Matrix Shape: {raw_counts.shape}")
    return raw_counts, meta_df