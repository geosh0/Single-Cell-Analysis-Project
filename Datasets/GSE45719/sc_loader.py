import pandas as pd
import os
import gzip
import re
import numpy as np
from typing import Tuple

def load_gse45719_data(
    expression_dir: str, 
    sra_path: str
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Loads GSE45719 data (RPKM single-cell data).
    1. Reads individual .txt.gz files.
    2. Merges them into a Gene x Sample matrix.
    3. Parses SRA metadata to create biological labels.

    Returns:
        Tuple[pd.DataFrame, pd.DataFrame]: (Expression Matrix, Metadata)
    """
    
    # --- PART 1: LOAD EXPRESSION FILES ---
    if not os.path.exists(expression_dir):
        raise FileNotFoundError(f"Directory not found: {expression_dir}")

    # Files usually look like: GSM1111111_expression.txt.gz
    files = [f for f in os.listdir(expression_dir) if f.endswith('expression.txt.gz')]
    print(f"Found {len(files)} expression files in {expression_dir}")

    gsm_pattern = re.compile(r'(GSM\d+)')
    all_series = []
    
    for filename in files:
        # Extract GSM ID
        match = gsm_pattern.search(filename)
        if not match: continue
        gsm_id = match.group(1)
        
        file_path = os.path.join(expression_dir, filename)
        
        try:
            # Smart loading: Check header for comments
            # GSE45719 files often have a header row, sometimes comments
            with gzip.open(file_path, 'rt') as f:
                header_line = f.readline().strip()
                # If header starts with #, remove it to get column names
                if header_line.startswith('#'):
                    header_line = header_line[1:]
                col_names = header_line.split('\t')

            # Read CSV using the parsed header
            df_temp = pd.read_csv(
                file_path, 
                sep='\t', 
                compression='gzip', 
                header=0, 
                names=col_names, 
                skiprows=1 # Skip the header line we just read manually
            )
            
            # Standardize columns (Gene_symbol, RPKM)
            # Adjust these keys if the specific txt file differs
            if 'Gene_symbol' in df_temp.columns and 'RPKM' in df_temp.columns:
                # Create a Series: Index=Gene, Value=RPKM, Name=SampleID
                series = df_temp.set_index('Gene_symbol')['RPKM']
                series.name = gsm_id
                
                # Handle duplicate genes in a single file (isoforms) by summing
                if series.index.duplicated().any():
                    series = series.groupby(level=0).sum()
                    
                all_series.append(series)
            
        except Exception as e:
            print(f"Warning: Failed to load {filename}. Error: {e}")

    if not all_series:
        raise ValueError("No valid data loaded.")

    print(f"Concatenating {len(all_series)} samples into matrix...")
    # Concatenate columns (axis=1). 'outer' join keeps all genes found in any file.
    expression_df = pd.concat(all_series, axis=1, join='outer')
    
    # Fill missing values (gene detected in sample A but not B = 0)
    expression_df = expression_df.fillna(0)
    print(f"Expression Matrix Shape: {expression_df.shape}")

    # --- PART 2: PARSE METADATA ---
    print("Parsing SRA Metadata...")
    try:
        sra_df = pd.read_csv(sra_path)
    except Exception as e:
        print(f"Error loading SRA CSV: {e}")
        return expression_df, pd.DataFrame(index=expression_df.columns)

    metadata_rows = []
    
    # Only create metadata for samples we actually loaded
    for gsm_id in expression_df.columns:
        # Find row in SRA
        sra_match = sra_df[sra_df['Sample Name'] == gsm_id]
        
        if sra_match.empty:
            metadata_rows.append({'GSM_ID': gsm_id, 'Combined_Label': 'Unknown'})
            continue
            
        row = sra_match.iloc[0]
        
        # Extract fields safely
        dev_stage = str(row.get('Developmental_Stage', 'Unknown'))
        source = str(row.get('source_name', 'Unknown'))
        strain = str(row.get('strain', 'Unknown'))
        
        # Logic: Create a clean "Stage" label
        # If source_name looks like a stage (e.g. "Zygot"), use it
        stage_label = source if source != 'Unknown' and len(source) < 20 else dev_stage
        stage_label = stage_label.replace(' ', '_').replace('(', '').replace(')', '')

        # Logic: Create a clean "Cross" label
        if "CAST" in strain and "C57BL" in strain:
            cross_label = "Hybrid_Cross"
        elif "C57BL" in strain:
            cross_label = "B6_Pure"
        elif "CAST" in strain:
            cross_label = "CAST_Pure"
        else:
            cross_label = "Other"

        combined = f"{stage_label}_{cross_label}"
        
        metadata_rows.append({
            'GSM_ID': gsm_id,
            'Stage': stage_label,
            'Cross': cross_label,
            'Combined_Label': combined
        })
        
    metadata_df = pd.DataFrame(metadata_rows).set_index('GSM_ID')
    
    # Ensure alignment
    metadata_df = metadata_df.reindex(expression_df.columns)
    
    print(f"Metadata prepared. Shape: {metadata_df.shape}")
    return expression_df, metadata_df