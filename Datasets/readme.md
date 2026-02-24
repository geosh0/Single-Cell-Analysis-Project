# 📊 Dataset Catalog & Engineering Challenges

Each of the seven datasets in this portfolio presented unique data engineering challenges. While the downstream analysis (via `shared_utils`) was standardized, the **Extract-Transform-Load (ETL)** process required custom logic for each study.

| Dataset | Type | The Engineering Challenge (ETL) |
| :--- | :--- | :--- |
| **GSE41265** | TPM | **Merging Data Streams:** Required merging a main TPM table with a separate Molecular Barcode (UMB) table while handling specific column renaming logic. |
| **GSE42268** | FPKM | **Filename Parsing:** Biological group information was embedded in the raw filenames rather than a metadata table, requiring regex extraction during loading. |
| **GSE45719** | RPKM | **Developmental Staging:** Involved parsing complex developmental stages (zygote to blastocyst) and cross-referencing them with genotype information. |
| **GSE52583** | FPKM | **Control Filtering:** The raw data included "Bulk" controls and "No-Cell" technical controls mixed with single cells. The loader had to identify and filter these out based on GSM IDs. |
| **GSE60361** | FPKM | **Format Pivoting:** Data was provided in "Long" format `.fpkm_tracking` files. The loader had to pivot these into a Wide Matrix (Genes x Cells) and filter based on the `technology_used` metadata column. |
| **GSE65528** | TPM | **Complex Parsing:** Metadata was locked inside a raw text `Series_Matrix.txt` file (GEO format). I wrote a custom parser to extract infection status and timepoints, then merged this with melted count tables. |
| **GSE74596** | Raw Counts | **SRA Aggregation:** A single biological cell was often split across multiple sequencing runs (SRR IDs). The loader implemented aggregation logic to sum reads and metadata across multiple runs to reconstruct the single-cell profile. |
