# Example SHARE-seq Input Data

This directory contains example input FASTQ files for the SHARE-seq pipeline. These files are synthetic and designed to demonstrate the expected file structure and barcode placement for SHARE-seq libraries.

## Files

For each library type (scRNA-seq and scATAC-seq), there are 4 FASTQ files:

### scRNA-seq
*   **`rna_R1.fastq.gz`**: Read 1 - Biological Sequence (cDNA).
*   **`rna_R2.fastq.gz`**: Read 2 - UMI (Unique Molecular Identifier).
*   **`rna_I1.fastq.gz`**: Index 1 - Cell Barcodes (99bp).
*   **`rna_I2.fastq.gz`**: Index 2 - Library/Sample Barcode (PKR).

### scATAC-seq
*   **`atac_R1.fastq.gz`**: Read 1 - Biological Sequence (Genomic Forward).
*   **`atac_R2.fastq.gz`**: Read 2 - Biological Sequence (Genomic Reverse).
*   **`atac_I1.fastq.gz`**: Index 1 - Cell Barcodes (99bp).
*   **`atac_I2.fastq.gz`**: Index 2 - Library/Sample Barcode (PKR).

## File Specifications

### scRNA-seq Details
*   **Read 1 (R1)**: Contains the transcript sequence.
*   **Read 2 (R2)**: Contains the 10bp UMI at the beginning of the read.
    *   **Structure**: `[UMI (10bp)]`
*   **Index 1 (I1)**: Contains the cell barcodes (last 99bp).
    *   **Structure**: `[Spacer/Linker] + [Barcodes (last 99bp)]`
*   **Index 2 (I2)**: Contains the 8bp sample index (PKR).

### scATAC-seq Details
*   **Read 1 (R1)**: Contains the forward genomic read.
*   **Read 2 (R2)**: Contains the reverse genomic read.
*   **Index 1 (I1)**: Contains the cell barcodes (last 99bp).
    *   **Structure**: `[Spacer/Linker] + [Barcodes (last 99bp)]`
*   **Index 2 (I2)**: Contains the 8bp sample index (PKR).

## Barcode Mapping (Last 99bp of I1)

The SHARE-seq pipeline expects the cell barcodes to be located at specific positions within the last 99bp of the Cell Barcode Read (Index 1).

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

## Other Files
*   **`whitelist.txt`**: A subset whitelist containing the valid barcode triplet used in the synthetic data.
*   **`generate_example_data.py`**: Python script used to generate these files.
