import pandas as pd
import os
import re
import numpy as np
from typing import Tuple

# --- HARDCODED CONTROLS (Specific to GSE52583) ---
BULK_CONTROL_GSMS = ["GSM1271882", "GSM1271944"]
NO_CELL_CONTROL_GSM = ["GSM1271883"]

def _load_raw_fpkm_files(data_dir: str) -> pd.DataFrame:
    """Internal helper to load individual FPKM tracking files."""
    if not os.path.isdir(data_dir):
        raise FileNotFoundError(f"Directory not found: {data_dir}")

    files = [f for f in os.listdir(data_dir) if f.endswith('.fpkm_tracking.gz')]
    print(f"Found {len(files)} FPKM files to process.")
    
    gsm_pattern = re.compile(r'(GSM\d+)')
    series_list = []
    
    for filename in files:
        # Extract GSM ID
        match = gsm_pattern.search(filename)
        if not match: continue
        gsm_id = match.group(1)
        
        file_path = os.path.join(data_dir, filename)
        
        try:
            # Load Data
            df = pd.read_csv(file_path, sep='\t', compression='gzip', header=0)
            
            # --- User's Column Logic ---
            gene_col, fpkm_col = None, None
            
            # Find Gene Column
            if 'tracking_id' in df.columns: gene_col = 'tracking_id'
            elif len(df.columns) > 0: gene_col = df.columns[0]
            
            # Find FPKM Column
            if 'FPKM' in df.columns: fpkm_col = 'FPKM'
            elif len(df.columns) > 9:
                # Heuristic check for Column 9 (Cufflinks standard)
                col_9 = df.columns[9]
                # Check if it looks numeric
                if pd.api.types.is_numeric_dtype(df[col_9]) or \
                   df[col_9].astype(str).str.contains('e', case=False).any():
                    fpkm_col = col_9
            
            if not gene_col or not fpkm_col:
                print(f"Skipping {gsm_id}: Could not identify columns.")
                continue
                
            # Clean and Store
            df[fpkm_col] = pd.to_numeric(df[fpkm_col], errors='coerce')
            df = df.dropna(subset=[gene_col, fpkm_col])
            
            # Handle duplicates (averaging)
            if df.duplicated(subset=[gene_col]).any():
                s = df.groupby(gene_col)[fpkm_col].mean()
            else:
                s = df.set_index(gene_col)[fpkm_col]
            
            s.name = gsm_id
            series_list.append(s)
            
        except Exception as e:
            print(f"Error loading {filename}: {e}")
            
    if not series_list:
        return pd.DataFrame()
        
    print(f"Concatenating {len(series_list)} samples...")
    combined_df = pd.concat(series_list, axis=1, join='outer')
    
    # Handle duplicates in the index (Gene IDs) if they exist across files
    if combined_df.index.has_duplicates:
        combined_df = combined_df.groupby(level=0).mean()
        
    return combined_df.fillna(0)

def load_gse52583_data(
    expression_dir: str, 
    sra_path: str,
    filter_single_cells: bool = True
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Master loader for GSE52583.
    Returns (Expression Matrix, Metadata).
    If filter_single_cells=True, removes Bulk/No-Cell controls.
    """
    # 1. Load Expression
    raw_df = _load_raw_fpkm_files(expression_dir)
    if raw_df.empty:
        raise ValueError("No expression data loaded.")
        
    # 2. Load SRA
    try:
        sra_df = pd.read_csv(sra_path)
    except Exception as e:
        print(f"Error loading SRA: {e}")
        sra_df = pd.DataFrame()
        
    # 3. Build Metadata
    metadata_rows = []
    
    for gsm_id in raw_df.columns:
        # Defaults
        age, genotype, instrument = "Unknown", "Unknown", "Unknown"
        
        # Get SRA info
        if not sra_df.empty:
            match = sra_df[sra_df['Sample Name'] == gsm_id]
            if not match.empty:
                age = str(match['AGE'].iloc[0])
                genotype = str(match['genotype'].iloc[0])
                instrument = str(match['Instrument'].iloc[0])

        # Determine Sample Type (Control Logic)
        if gsm_id in BULK_CONTROL_GSMS:
            s_type = "Bulk_200cell"
        elif gsm_id in NO_CELL_CONTROL_GSM:
            s_type = "No_Cell_Control"
        else:
            s_type = "Single_Cell"
            
        # Clean Labels
        clean_age = age.replace(' ', '_').replace('.', '_')
        clean_genotype = "SftpcGFP" if "Sftpc-Cre" in genotype else ("WT" if "wild type" in genotype.lower() else "Other")
        
        combined_label = f"{clean_age}_{clean_genotype}"
        
        metadata_rows.append({
            'GSM_ID': gsm_id,
            'Sample_Type': s_type,
            'Age': clean_age,
            'Genotype': clean_genotype,
            'Combined_Label': combined_label
        })
        
    meta_df = pd.DataFrame(metadata_rows).set_index('GSM_ID')
    
    # 4. Filter (Optional but recommended)
    if filter_single_cells:
        valid_samples = meta_df[meta_df['Sample_Type'] == 'Single_Cell'].index
        # Intersect with loaded data
        valid_samples = valid_samples.intersection(raw_df.columns)
        
        print(f"Filtering: Keeping {len(valid_samples)} Single Cells (Dropped {len(raw_df.columns) - len(valid_samples)} controls)")
        
        raw_df = raw_df[valid_samples]
        meta_df = meta_df.loc[valid_samples]
        
    return raw_df, meta_df