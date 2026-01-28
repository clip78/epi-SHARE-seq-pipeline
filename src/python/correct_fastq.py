#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Correct fastq
"""

import argparse
import xopen
from collections import deque
import multiprocessing

def parse_arguments():
    parser = argparse.ArgumentParser(description="Perform barcode error correction on read 2 FASTQ file; write corrected barcodes into read names of both read 1 and read 2 FASTQ files; generate QC statistics file.")
    parser.add_argument("input_read1_fastq_file", help="Filename for uncorrected input read 1 FASTQ file")
    parser.add_argument("input_read2_fastq_file", help="Filename for uncorrected input read 2 FASTQ file")
    parser.add_argument("output_read1_fastq_file", help="Filename for corrected output read 1 FASTQ file")
    parser.add_argument("output_read2_fastq_file", help="Filename for corrected output read 2 FASTQ file")
    parser.add_argument("whitelist_file", help="Filename for whitelisted combinations of R1R2R3 barcodes, one per line")
    parser.add_argument("sample_type", choices=["ATAC", "RNA"], help="Sample modality")
    parser.add_argument("prefix", help="Prefix for naming output QC txt file")
    parser.add_argument("pkr", nargs="?", help="PKR name")
    parser.add_argument("--barcode_fastq", help="Filename for separate barcode FASTQ file", default=None)
    
    return parser.parse_args()

def get_barcodes(whitelist_file):
    """
    Read barcode whitelist file, split into R1, R2, and R3 barcodes
    """
    r1_barcodes, r2_barcodes, r3_barcodes = set(), set(), set()
    with open(whitelist_file) as f:
        for line in f:
            r1_barcodes.add(line[:8])
            r2_barcodes.add(line[8:16])
            r3_barcodes.add(line[16:24])
    
    return r1_barcodes, r2_barcodes, r3_barcodes

def check_putative_barcode(barcode_str, barcode_set, quality_str):
    """
    Procedure: check exact match of barcode, then 1 mismatch, then 1bp left/right shift
    """

    # Helper function to find match (exact or 1-mismatch)
    def find_match(putative_barcode, barcode_set):
        # 1. Exact match
        if putative_barcode in barcode_set:
            return putative_barcode

        # 2. Mismatch search
        corrected_barcode = None
        # Generate neighbors of the *read* barcode
        for i, base in enumerate(putative_barcode):
            for x in 'ACGTN':
                if base != x:
                    neighbor = putative_barcode[:i] + x + putative_barcode[i + 1:]
                    if neighbor in barcode_set:
                        if corrected_barcode is not None and corrected_barcode != neighbor:
                            # Ambiguous match
                            return None
                        corrected_barcode = neighbor
        return corrected_barcode

    # Check exact location first
    value = find_match(barcode_str[1:9], barcode_set)
    quality = quality_str[1:9]
    if value is None:
        # Check 1bp shift left
        value = find_match(barcode_str[:8], barcode_set)
        quality = quality_str[:8]
        if value is None:
            # check 1bp shift right
            # round 3 is shorter so add "N" for those
            if len(barcode_str) < 10: 
                value = find_match(barcode_str[2:]+"N", barcode_set)
                quality = quality_str[2:]+"F"
            else:
                value = find_match(barcode_str[2:], barcode_set)
                quality = quality_str[2:]
                    
    return value, quality

# Global variables for worker processes
g_r1_barcodes = None
g_r2_barcodes = None
g_r3_barcodes = None
g_sample_type = None
g_pkr = None
g_has_barcode_file = None

def worker_init(r1_barcodes, r2_barcodes, r3_barcodes, sample_type, pkr, has_barcode_file):
    global g_r1_barcodes, g_r2_barcodes, g_r3_barcodes
    global g_sample_type, g_pkr, g_has_barcode_file
    g_r1_barcodes = r1_barcodes
    g_r2_barcodes = r2_barcodes
    g_r3_barcodes = r3_barcodes
    g_sample_type = sample_type
    g_pkr = pkr
    g_has_barcode_file = has_barcode_file

def process_chunk(chunk):
    """
    Process a chunk of reads.
    chunk: list of tuples ((name1, seq1, qual1), (name2, seq2, qual2), (name_bc, seq_bc, qual_bc) or None)
    Returns: (list_of_r1_out_strings, list_of_r2_out_strings, stats_dictionary)
    """
    out_r1 = []
    out_r2 = []
    stats = {"match": 0, "mismatch": 0, "poly_g": 0, "read2_poly_g": 0}

    for record in chunk:
        r1_rec, r2_rec, bc_rec = record
        name1, sequence1, quality1 = r1_rec
        name2, sequence2, quality2 = r2_rec

        sequence_bc = bc_rec[1] if bc_rec else None
        quality_bc = bc_rec[2] if bc_rec else None

        # Logic from original process_fastqs
        if g_has_barcode_file:
            read_2_barcode_sequence = sequence_bc
            read_2_barcode_quality = quality_bc
            if len(read_2_barcode_sequence) > 99:
                 read_2_barcode_sequence = read_2_barcode_sequence[-99:]
                 read_2_barcode_quality = read_2_barcode_quality[-99:]
        else:
            read_2_barcode_sequence = sequence2[-99:]
            read_2_barcode_quality = quality2[-99:]

        if len(read_2_barcode_sequence) < 99:
            stats["mismatch"] += 1
            continue

        r1_str, r2_str, r3_str = read_2_barcode_sequence[14:24], read_2_barcode_sequence[52:62], read_2_barcode_sequence[90:99]
        q1_str, q2_str, q3_str = read_2_barcode_quality[14:24], read_2_barcode_quality[52:62], read_2_barcode_quality[90:99]

        r1 = r2 = r3 = None
        r1, q1 = check_putative_barcode(r1_str, g_r1_barcodes, q1_str)
        r2, q2 = check_putative_barcode(r2_str, g_r2_barcodes, q2_str)
        r3, q3 = check_putative_barcode(r3_str, g_r3_barcodes, q3_str)

        if sequence2[:10] == "G"*10:
            stats["read2_poly_g"] += 1

        elif r1 and r2 and r3:
            stats["match"] += 1

            if g_sample_type == "RNA":
                corrected_header = name1.split(" ")[0] + "_" + ",".join(filter(None, [r1, r2, r3, g_pkr])) + "_" + sequence2[:10]
                corrected_read1 = f"{corrected_header}\n{sequence1}\n+\n{quality1}\n"
                out_r1.append(corrected_read1)

                corrected_sequence2 = r1 + r2 + r3 + sequence2[:10]
                corrected_quality2 = q1 + q2 + q3 + quality2[:10]
                corrected_read2 = f"{corrected_header}\n{corrected_sequence2}\n+\n{corrected_quality2}\n"
                out_r2.append(corrected_read2)

            elif g_sample_type == "ATAC":
                corrected_header = name1.split(" ")[0] + "_" + ",".join(filter(None, [r1, r2, r3, g_pkr]))
                corrected_read1 = f"{corrected_header}\n{sequence1}\n+\n{quality1}\n"
                out_r1.append(corrected_read1)

                if g_has_barcode_file:
                    sequence2_out = sequence2
                    quality2_out = quality2
                else:
                    sequence2_out = sequence2[:-99]
                    quality2_out = quality2[:-99]

                corrected_read2 = f"{corrected_header}\n{sequence2_out}\n+\n{quality2_out}\n"
                out_r2.append(corrected_read2)

        elif "G"*8 in r1_str and "G"*8 in r2_str and "G"*8 in r3_str:
            stats["poly_g"] += 1

        else:
            stats["mismatch"] += 1

    return out_r1, out_r2, stats

def chunk_generator(iterators, read1_fh, read2_fh, barcode_fh, chunk_size=1000):
    """
    Yields chunks of reads.
    """
    while True:
        chunk = []
        try:
            for _ in range(chunk_size):
                reads = next(iterators)

                r1_name = reads[0].strip()
                r2_name = reads[1].strip()
                bc_name = reads[2].strip() if barcode_fh else None

                r1_seq = next(read1_fh).strip()
                r2_seq = next(read2_fh).strip()
                bc_seq = next(barcode_fh).strip() if barcode_fh else None

                # skip +
                next(read1_fh)
                next(read2_fh)
                if barcode_fh: next(barcode_fh)

                r1_qual = next(read1_fh).strip()
                r2_qual = next(read2_fh).strip()
                bc_qual = next(barcode_fh).strip() if barcode_fh else None

                chunk.append(((r1_name, r1_seq, r1_qual), (r2_name, r2_seq, r2_qual), (bc_name, bc_seq, bc_qual) if barcode_fh else None))
        except StopIteration:
            if chunk:
                yield chunk
            return

        yield chunk

def process_fastqs(input_read1_fastq_file, input_read2_fastq_file,
                  output_read1_fastq_file, output_read2_fastq_file,
                  r1_barcodes, r2_barcodes, r3_barcodes,
                  sample_type, pkr, prefix, barcode_fastq_file=None):

    # QC counters
    total_stats = {"match": 0, "mismatch": 0, "poly_g": 0, "read2_poly_g": 0}

    read1_out_writer = xopen.xopen(output_read1_fastq_file, mode = 'w')
    read2_out_writer = xopen.xopen(output_read2_fastq_file, mode ='w')

    # Determine file openers
    read1_fh = xopen.xopen(input_read1_fastq_file, mode="r", threads=8)
    read2_fh = xopen.xopen(input_read2_fastq_file, mode="r", threads=8)
    barcode_fh = xopen.xopen(barcode_fastq_file, mode="r", threads=8) if barcode_fastq_file else None

    try:
        if barcode_fh:
             iterators = zip(read1_fh, read2_fh, barcode_fh)
        else:
             iterators = zip(read1_fh, read2_fh)

        # Initialize pool
        # Determine number of processes. Use typical default or reasonable number.
        # Since I/O is in main, maybe CPU count is fine.
        pool = multiprocessing.Pool(initializer=worker_init,
                                    initargs=(r1_barcodes, r2_barcodes, r3_barcodes, sample_type, pkr, bool(barcode_fh)))

        try:
            generator = chunk_generator(iterators, read1_fh, read2_fh, barcode_fh, chunk_size=1000)
            
            # Using imap for ordered results
            for r1_out, r2_out, stats in pool.imap(process_chunk, generator):
                # Write output
                if r1_out:
                    read1_out_writer.write("".join(r1_out))
                if r2_out:
                    read2_out_writer.write("".join(r2_out))
                
                # Aggregate stats
                for k, v in stats.items():
                    total_stats[k] += v

        finally:
            pool.close()
            pool.join()

    finally:
        read1_fh.close()
        read2_fh.close()
        if barcode_fh:
            barcode_fh.close()
        read1_out_writer.close()
        read2_out_writer.close()

    
    # write QC stats
    with open(f"{prefix}_barcode_qc.txt", "w") as f:
        fields = ["library", "match", "mismatch", "poly_G_barcode", "poly_G_in_first_10bp"]
        f.write("\t".join(fields) + "\n")
        f.write("%s\t%s\t%s\t%s\t%s" % (prefix, total_stats["match"], total_stats["mismatch"], total_stats["poly_g"], total_stats["read2_poly_g"]))

def main():
    args = parse_arguments()
    input_read1_fastq_file = getattr(args, "input_read1_fastq_file")
    input_read2_fastq_file = getattr(args, "input_read2_fastq_file")
    output_read1_fastq_file = getattr(args, "output_read1_fastq_file")
    output_read2_fastq_file = getattr(args, "output_read2_fastq_file")
    whitelist_file = getattr(args, "whitelist_file")
    sample_type = getattr(args, "sample_type")
    prefix = getattr(args, "prefix")
    pkr = getattr(args, "pkr")
    barcode_fastq_file = getattr(args, "barcode_fastq")
    
    # read whitelist, get lists of barcodes
    (r1_barcodes, r2_barcodes, r3_barcodes) = get_barcodes(whitelist_file)

    # write corrected FASTQs and QC stats
    process_fastqs(input_read1_fastq_file, input_read2_fastq_file,
                   output_read1_fastq_file, output_read2_fastq_file,
                   r1_barcodes, r2_barcodes, r3_barcodes,
                   sample_type, pkr, prefix, barcode_fastq_file)

if __name__ == "__main__":
    main()
