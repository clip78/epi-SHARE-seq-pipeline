# Example SHARE-seq Input Data

This directory contains example input FASTQ files for the SHARE-seq pipeline. These files are synthetic and designed to demonstrate the expected file structure and barcode placement.

## Files

*   **`rna_R1.fastq.gz`**: Synthetic scRNA-seq Read 1 (cDNA).
*   **`rna_R2.fastq.gz`**: Synthetic scRNA-seq Read 2 (UMI + Barcodes).
*   **`atac_R1.fastq.gz`**: Synthetic scATAC-seq Read 1 (Genomic).
*   **`atac_R2.fastq.gz`**: Synthetic scATAC-seq Read 2 (Genomic + Barcodes).
*   **`whitelist.txt`**: A subset whitelist containing the valid barcode triplet used in the synthetic data.
*   **`generate_example_data.py`**: Python script used to generate these files.

## File Specifications

### scRNA-seq
*   **Read 1**: Contains the biological sequence (transcript).
*   **Read 2**: Contains the UMI and Cell Barcodes.
    *   **Structure**: `[UMI (10bp)] + [Spacer/Linker] + [Barcodes (last 99bp)]`
    *   The pipeline extracts the **first 10bp** as the UMI.
    *   The pipeline extracts the **last 99bp** to find cell barcodes.

### scATAC-seq
*   **Read 1**: Contains the forward biological read (genomic).
*   **Read 2**: Contains the reverse biological read (genomic) followed by cell barcodes.
    *   **Structure**: `[Genomic Sequence] + [Barcodes (last 99bp)]`
    *   The pipeline strips the **last 99bp** to isolate the genomic sequence.
    *   The **last 99bp** are processed to find cell barcodes.

## Barcode Mapping (Last 99bp)

The SHARE-seq pipeline expects the cell barcodes to be located at specific positions within the last 99bp of the read (Read 2 for 2-file input, or Barcode Read for 3-file input).

Assuming a 0-indexed 99bp sequence:

1.  **Barcode 1 (R1)**: 8bp length.
    *   Located around index 15.
    *   Exact extraction window in code: `[14:24]`.
    *   Expected position: Indices `15` to `22`.

2.  **Barcode 2 (R2)**: 8bp length.
    *   Located around index 53.
    *   Exact extraction window in code: `[52:62]`.
    *   Expected position: Indices `53` to `60`.

3.  **Barcode 3 (R3)**: 8bp length.
    *   Located around index 91.
    *   Exact extraction window in code: `[90:99]`.
    *   Expected position: Indices `91` to `98`.

The `whitelist.txt` file provided contains the barcode triplet `AAAAAAAA` + `CCCCCCCC` + `GGGGGGGG` (24bp) which matches the valid reads in the example files.
