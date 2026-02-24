import pandas as pd
import numpy as np
from typing import Optional, Dict

def _parse_condition_row(row: pd.Series) -> Dict[str, str]:
    """
    Internal helper to parse a single row of SRA metadata into a descriptive label.
    Logic is specific to the GSE41265 dataset annotations.
    """
    # 1. Extract values safely and ensure they are strings
    cell_count_val = str(row.get('cell_count', 'Unknown')).strip()
    protocol_val = str(row.get('protocol', '')).strip()
    source_name_val = str(row.get('source_name', '')).strip()
    treatment_val = str(row.get('treatment', '')).strip()

    # 2. Determine Cell Type (Part 1)
    # Using 'in' allows for partial matches
    if "1 cell" in cell_count_val:
        label_part1 = "SC" # Single Cell
    elif "10,000 cells" in cell_count_val:
        label_part1 = "Bulk_10k"
    else:
        label_part1 = "Unknown_Cells"

    # 3. Determine Protocol (Part 2)
    if "molecular barcodes (MB)" in protocol_val:
        label_part2 = "MB_Protocol"
    elif protocol_val == "" or protocol_val == "nan": 
        # Empty implies standard SMARTer for this specific dataset
        label_part2 = "Std_Protocol"
    else:
        label_part2 = "Other_Protocol"

    # 4. Determine Stimulation (Part 3)
    # Complex logic based on treatment and source name columns
    if "BMDC (4h LPS stim)" in source_name_val:
        label_part3 = "LPS_4h"
    elif treatment_val == "LPS-stimulation" and (source_name_val == "" or source_name_val == "nan"):
        # Heuristics for unspecified duration based on other attributes
        if label_part1 == "SC" and label_part2 == "Std_Protocol":
             label_part3 = "LPS_4h" # Assumption based on dataset paper
        elif label_part1 == "Bulk_10k":
             label_part3 = "LPS_4h"
        elif label_part1 == "SC" and label_part2 == "MB_Protocol":
             label_part3 = "LPS_4h"
        else:
             label_part3 = "LPS_UnknownDetail"
    else: 
        label_part3 = "Unstimulated" # Default fallback

    # 5. Return the constructed dictionary
    return {
        'Combined_Label': f"{label_part1}_{label_part2}_{label_part3}",
        'Cell_Count_Detail': cell_count_val,
        'Protocol_Detail': label_part2,
        'Stimulation_Detail': label_part3
    }

def prepare_sample_metadata(
    sra_run_table_df: pd.DataFrame, 
    qc_metrics_df: Optional[pd.DataFrame] = None, 
    pca_df: Optional[pd.DataFrame] = None
) -> pd.DataFrame:
    """
    Prepares a comprehensive metadata table for samples.
    Merges SRA technical details with the analytical index (QC or PCA).
    
    Args:
        sra_run_table_df: Raw metadata from SRA (SraRunTable.csv).
        qc_metrics_df: DataFrame containing QC metrics.
        pca_df: DataFrame containing PCA coordinates (defines final valid samples).
        
    Returns:
        pd.DataFrame: Aligned metadata with 'Combined_Label'.
    """
    
    # --- Input Validation ---
    if sra_run_table_df.empty:
        print("Error: SRA run table is empty.")
        return pd.DataFrame()

    print(f"--- Preparing Sample Metadata ---")
    
    # --- 1. Consolidate SRA info per Sample Name (GSM ID) ---
    # We rename 'Sample Name' to 'GSM_ID' to match our expression matrix columns
    sra_clean = sra_run_table_df.copy()
    if 'Sample Name' in sra_clean.columns:
        sra_clean = sra_clean.rename(columns={'Sample Name': 'GSM_ID'})
    
    # We set GSM_ID as index to facilitate merging
    # Grouping by ID ensures we don't have duplicates (taking the first entry found)
    sample_attributes = sra_clean.groupby('GSM_ID').first()

    # --- 2. Generate Labels (Iterate and Parse) ---
    conditions = []
    for gsm_id, row in sample_attributes.iterrows():
        parsed_data = _parse_condition_row(row)
        # Add the ID back so we can set it as index later
        parsed_data['GSM_ID'] = gsm_id 
        conditions.append(parsed_data)
        
    # Create the conditions dataframe
    conditions_df = pd.DataFrame(conditions).set_index('GSM_ID')

    # --- 3. Determine the "Master Index" ---
    # We want metadata only for the cells that exist in our analysis (PCA or QC)
    target_index = None
    
    if pca_df is not None and not pca_df.empty:
        target_index = pca_df.index
        print(f"Aligning metadata to PCA index ({len(target_index)} samples)")
    elif qc_metrics_df is not None and not qc_metrics_df.empty:
        target_index = qc_metrics_df.index
        print(f"Aligning metadata to QC index ({len(target_index)} samples)")
    else:
        print("Warning: No QC or PCA data provided. Returning all SRA samples.")
        return conditions_df

    # --- 4. Merge ---
    # We use reindex to strictly keep only our target samples
    # This automatically fills missing samples with NaN
    final_metadata = conditions_df.reindex(target_index)

    # Fill NaNs for any samples in our data but missing from SRA
    final_metadata['Combined_Label'] = final_metadata['Combined_Label'].fillna('Unknown_Condition')
    
    # If QC metrics were provided, append them to the metadata for convenience
    if qc_metrics_df is not None:
        final_metadata = final_metadata.join(qc_metrics_df, how='left')

    print(f"Metadata prepared. Final Shape: {final_metadata.shape}")
    return final_metadata