import pandas as pd
import os
import glob
import re
from typing import Tuple

def parse_series_matrix(series_matrix_path: str) -> pd.DataFrame:
    """
    Parses the GSE65528_series_matrix.txt file to extract sample metadata.
    (Refactored from your provided code to be a clean helper function)
    """
    samples_metadata = []
    current_sample_info = {}
    temp_characteristics = {} 
    temp_descriptions = {}
    gsm_ids_list = []

    try:
        with open(series_matrix_path, 'r', encoding='utf-8') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("!series_matrix_table"): continue

                parts = line.split("\t")
                header = parts[0]
                values = [val.strip('"') for val in parts[1:]]

                if header == "!Sample_geo_accession":
                    gsm_ids_list = values
                    for i, gsm_id in enumerate(gsm_ids_list):
                        samples_metadata.append({'GSM_ID': gsm_id, 'col_idx': i})
                    continue

                if header == "!Sample_characteristics_ch1":
                    for i, entry in enumerate(values):
                        if i not in temp_characteristics: temp_characteristics[i] = {}
                        if ":" in entry:
                            key, value = entry.split(":", 1)
                            key_clean = key.strip().lower().replace(" ", "_").replace("(", "").replace(")", "").replace("?", "")
                            temp_characteristics[i][key_clean] = value.strip()
                    continue
                
                if header == "!Sample_description":
                    for i, desc_entry in enumerate(values):
                        if i not in temp_descriptions: temp_descriptions[i] = {}
                        
                        if "Column " in desc_entry and " of processed data file" in desc_entry:
                            temp_descriptions[i]['plate_column'] = desc_entry.replace("Column ", "").replace(" of processed data file", "")
                        elif "_table.txt" in desc_entry: 
                            filename_only = os.path.basename(desc_entry)
                            if "_p1_" in filename_only: temp_descriptions[i]['plate_id'] = 'p1'
                            elif "_p2_" in filename_only: temp_descriptions[i]['plate_id'] = 'p2'
                    continue

        # Combine Data
        for sample_dict in samples_metadata:
            idx = sample_dict['col_idx']
            if idx in temp_characteristics: sample_dict.update(temp_characteristics[idx])
            if idx in temp_descriptions: sample_dict.update(temp_descriptions[idx])
        
        df = pd.DataFrame(samples_metadata)
        
        # Create Merge Key: "p1_1", "p2_12", etc.
        if 'plate_id' in df.columns and 'plate_column' in df.columns:
            df['cell_id_merge'] = df['plate_id'] + "_" + df['plate_column']
        else:
            print("Warning: Could not create cell_id_merge from Series Matrix.")

        # Derive Infection Status (Your logic)
        def get_status(row):
            phrodo = str(row.get('phrodo_positive', '')).lower() == 'yes'
            gfp = str(row.get('gfp_positive', '')).lower() == 'yes'
            if phrodo and gfp: return "Live_Bacteria"
            elif phrodo: return "Dead_Bacteria"
            else: return "No_Bacteria"

        df['Infection_Status'] = df.apply(get_status, axis=1)
        
        # Clean up columns
        if 'time_after_salmonella_exposure_h' in df.columns:
            df['TimePoint'] = pd.to_numeric(df['time_after_salmonella_exposure_h'], errors='coerce').fillna(0).astype(int)
        
        return df.set_index('cell_id_merge')

    except Exception as e:
        print(f"Error parsing Series Matrix: {e}")
        return pd.DataFrame()

def load_gse65528_data(
    data_dir: str, 
    series_matrix_path: str
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Loads Counts files, Pivots them to Matrix, and Merges with parsed Metadata.
    Returns: (Expression Matrix, Metadata)
    """
    
    # --- 1. Load & Melt Count Files ---
    count_files = glob.glob(os.path.join(data_dir, '*TPM_table.txt.gz'))
    print(f"Found {len(count_files)} count files.")
    
    melted_chunks = []
    
    for filename in count_files:
        try:
            # Extract Plate ID (p1 or p2)
            if "_p1_" in filename: plate_id = "p1"
            elif "_p2_" in filename: plate_id = "p2"
            else: continue
            
            # Read Matrix (Genes x 96 Wells)
            # Columns are just numbers "1", "2"... we need to prefix them with plate_id
            df = pd.read_csv(filename, sep="\t", index_col=0, compression='gzip')
            
            # Rename columns to match metadata logic: "p1_1", "p1_2"
            # Note: The original file has columns like "1", "2". 
            df.columns = [f"{plate_id}_{col}" for col in df.columns]
            
            # Melt to long format temporarily to stack p1 and p2
            # We reset index to keep Gene IDs
            df_melted = df.reset_index().melt(id_vars=df.index.name, var_name='cell_id_merge', value_name='TPM')
            
            # Standardize Gene Column Name
            gene_col = df_melted.columns[0]
            df_melted.rename(columns={gene_col: 'Gene'}, inplace=True)
            
            melted_chunks.append(df_melted)
            
        except Exception as e:
            print(f"Error loading {filename}: {e}")
            
    if not melted_chunks:
        raise ValueError("No count data loaded.")
        
    # Combine p1 and p2 data
    full_long_df = pd.concat(melted_chunks, ignore_index=True)
    
    # --- 2. Pivot to Matrix (Genes x Cells) ---
    print("Pivoting to Expression Matrix...")
    expression_matrix = full_long_df.pivot_table(
        index='Gene', 
        columns='cell_id_merge', 
        values='TPM', 
        fill_value=0
    )
    print(f"Matrix Shape: {expression_matrix.shape}")

    # --- 3. Parse Metadata ---
    print("Parsing Metadata...")
    meta_df = parse_series_matrix(series_matrix_path)
    
    # Filter Metadata to match available cells in Matrix
    common_cells = expression_matrix.columns.intersection(meta_df.index)
    meta_df = meta_df.loc[common_cells]
    expression_matrix = expression_matrix[common_cells]
    
    # Create Combined Label for Clustering
    # Combine TimePoint and Infection Status
    meta_df['Combined_Label'] = meta_df['TimePoint'].astype(str) + "h_" + meta_df['Infection_Status']
    
    print(f"Final Aligned Data: {expression_matrix.shape} matrix, {meta_df.shape} metadata.")
    
    return expression_matrix, meta_df