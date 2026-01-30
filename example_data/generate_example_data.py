import gzip
import os

def write_fastq(filename, sequences, names):
    with gzip.open(filename, 'wt') as f:
        for i, (seq, name) in enumerate(zip(sequences, names)):
            # Create a quality string of the same length as the sequence
            qual = 'F' * len(seq)
            f.write(f"@{name}\n{seq}\n+\n{qual}\n")

# Setup
os.makedirs("example_data", exist_ok=True)

# Whitelist
r1_bc = "AAAAAAAA"
r2_bc = "CCCCCCCC"
r3_bc = "GGGGGGGG"
whitelist_line = r1_bc + r2_bc + r3_bc
with open("example_data/whitelist.txt", "w") as f:
    f.write(whitelist_line + "\n")

# Construct 99bp barcode sequence
# Indices: 0-based.
# R1 at 15. R2 at 53. R3 at 91.
base_99 = ["T"] * 99
for i, char in enumerate(r1_bc):
    base_99[15 + i] = char
for i, char in enumerate(r2_bc):
    base_99[53 + i] = char
for i, char in enumerate(r3_bc):
    base_99[91 + i] = char

valid_barcode_seq = "".join(base_99)

# Invalid barcode sequence (mutated)
invalid_barcode_seq = valid_barcode_seq.replace("A", "G").replace("C", "T").replace("G", "A")

# Common Names
names = ["READ1", "READ2"]

# --- RNA ---
# R1: cDNA (Biological)
rna_r1_seqs = ["ACGTACGTACGTACGTACGT", "TGCATGCATGCATGCATGCA"]

# R2: UMI (10bp)
# Note: The pipeline expects the UMI to be in R2.
rna_r2_seqs = [
    "NNNNNNNNNN", # UMI for Read 1
    "NNNNNNNNNN"  # UMI for Read 2
]

# Barcode Read (I2)
rna_bc_seqs = [
    valid_barcode_seq,
    invalid_barcode_seq
]

write_fastq("example_data/rna_R1.fastq.gz", rna_r1_seqs, names)
write_fastq("example_data/rna_R2.fastq.gz", rna_r2_seqs, names)
write_fastq("example_data/rna_barcode.fastq.gz", rna_bc_seqs, names)


# --- ATAC ---
# R1: Genomic Read 1
atac_r1_seqs = ["GCTAGCTAGCTAGCTAGCTA", "CGATCGATCGATCGATCGAT"]

# R2: Genomic Read 2
atac_r2_seqs = ["ATATATATATATATATATAT", "GCGCGCGCGCGCGCGCGCGC"]

# Barcode Read (I2)
atac_bc_seqs = [
    valid_barcode_seq,
    invalid_barcode_seq
]

write_fastq("example_data/atac_R1.fastq.gz", atac_r1_seqs, names)
write_fastq("example_data/atac_R2.fastq.gz", atac_r2_seqs, names)
write_fastq("example_data/atac_barcode.fastq.gz", atac_bc_seqs, names)

print("Files generated.")
