import gzip
import os

def write_fastq(filename, sequences, names):
    with gzip.open(filename, 'wt') as f:
        for i, (seq, name) in enumerate(zip(sequences, names)):
            f.write(f"@{name}\n{seq}\n+\n{'F'*len(seq)}\n")

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

# RNA
rna_r1_seqs = ["ACGTACGTACGTACGTACGT", "TGCATGCATGCATGCATGCA"]
rna_names = ["READ1", "READ2"]
rna_r2_seqs = [
    "NNNNNNNNNN" + valid_barcode_seq, # UMI + Valid
    "NNNNNNNNNN" + invalid_barcode_seq # UMI + Invalid
]

write_fastq("example_data/rna_R1.fastq.gz", rna_r1_seqs, rna_names)
write_fastq("example_data/rna_R2.fastq.gz", rna_r2_seqs, rna_names)

# ATAC
atac_r1_seqs = ["GCTAGCTAGCTAGCTAGCTA", "CGATCGATCGATCGATCGAT"]
atac_names = ["READ1", "READ2"]
atac_r2_seqs = [
    "ATATATATATATATATATAT" + valid_barcode_seq, # Genomic + Valid
    "ATATATATATATATATATAT" + invalid_barcode_seq # Genomic + Invalid
]

write_fastq("example_data/atac_R1.fastq.gz", atac_r1_seqs, atac_names)
write_fastq("example_data/atac_R2.fastq.gz", atac_r2_seqs, atac_names)

print("Files generated.")
