#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Make uncorrected R1 and R2 FASTQs from unmapped BAM file.
"""

import argparse
import pysam

def parse_arguments():
    parser = argparse.ArgumentParser(description="Make R1 and R2 FASTQs from unmapped BAM file")
    parser.add_argument("bam_file", help="Filename for unmapped BAM file")
    parser.add_argument("pkr", help="PKR name")
    parser.add_argument("prefix", help="Prefix for naming output FASTQ files")
    parser.add_argument("r1_barcode_set_file", help="File containing biosample splits in R1 barcodes, one split per lane")
    parser.add_argument("--r2_barcode_file", help= "File containing R2 barcodes, one line")
    parser.add_argument("--r3_barcode_file", help= "File containing R3 barcodes, one line")
    
    return parser.parse_args()

# DNA base complements
TRANS_TABLE = str.maketrans("ATCGN", "TAGCN")

def reverse_complement(sequence):
    """
    Return reverse complement of DNA sequence.
    """
    return sequence.translate(TRANS_TABLE)[::-1]

class BufferedFastqWriter:
    """
    Buffered writer for FASTQ files.
    """
    def __init__(self, filename, buffer_size=65536):
        self.file = open(filename, "w")
        self.buffer = []
        self.current_size = 0
        self.buffer_size = buffer_size

    def write(self, name, sequence, quality):
        entry = f"@{name}\n{sequence}\n+\n{quality}\n"
        self.buffer.append(entry)
        self.current_size += len(entry)
        if self.current_size >= self.buffer_size:
            self.flush()

    def flush(self):
        if self.buffer:
            self.file.write("".join(self.buffer))
            self.buffer = []
            self.current_size = 0

    def close(self):
        self.flush()
        self.file.close()

def check_putative_barcode(barcode_str, barcode_matching_dict):
    """
    Procedure: check exact match of barcode, then 1 mismatch, then 1bp left/right shift
    """
    # check exact location first
    exact = True
    corrected_barcode = barcode_matching_dict.get(barcode_str[1:9]) 
    if corrected_barcode is None:
        exact = False
        # check 1bp shift left
        corrected_barcode = barcode_matching_dict.get(barcode_str[:8]) 
        if corrected_barcode is None:
            exact = False
            # check 1bp shift right
            corrected_barcode = barcode_matching_dict.get(barcode_str[2:])
    return corrected_barcode, exact

def create_barcode_matching_dict(barcode_list):
    """
    Create dictionary mapping barcodes to themselves (exact matches), as well as
    to their 1bp mismatch possibilities
    """
    barcode_matching_dict  = dict() # {potential match: original barcode}
    for barcode in barcode_list:
        barcode_matching_dict[barcode] = barcode # exact match
        for i, base in enumerate(barcode):
            for x in 'ACGTN':
                if base != x:
                    # add mismatch possibilities at pos i
                    mismatch = barcode[:i] + x + barcode[i+1:]
                    barcode_matching_dict[mismatch] = barcode
    return barcode_matching_dict

def create_barcode_subset_dict(file_path):
    """
    Create dictionary mapping barcodes to their R1 barcode subsets
    """    
    with open(file_path) as f:
        barcode_subset_dict = dict() # {barcode: subset}
        for line in f:
            line = line.strip().split()
            subset = line[0] # first entry of line is subset name  
            for barcode in line[1:]:
                barcode_subset_dict[barcode] = subset
    return barcode_subset_dict

def write_fastqs(bam_file, read_1_writers, read_2_writers, r1_barcode_subset_dict, r1_barcode_matches, prefix):
    """
    Get reads from BAM file and error-correct R1 barcode to determine which R1 barcode subset's
    FASTQ to write the read to. Barcodes written into FASTQs are not error-corrected. 
    
    read_1_writers and read_2_writers are dictionaries mapping R1 barcode subset names
    to corresponding BufferedFastqWriter objects.
    """ 
    bam = pysam.Samfile(bam_file, "rb", check_sq=False)
    query_name = read_1 = None
    exact_match = nonexact_match = nonmatch = poly_g_barcode = 0

    for read in bam:
        if read.is_read1:
            # save and continue processing
            read_1 = read
            query_name = read_1.query_name
            
        else:
            # check that query names are the same
            if read.query_name == query_name:
                read_2 = read
                # get tags
                barcode_tag = read.get_tag("RX")
                quality_tag = read.get_tag("QX")
                
                # get 10bp sequence containing R1 barcode (additional 1bp padding used for checking shifts)
                r1_barcode_window = barcode_tag[14:24] 
                # get error-corrected R1 barcode
                r1_barcode, exact = check_putative_barcode(r1_barcode_window, r1_barcode_matches)
                 
                # write reads to appropriate FASTQ files by checking which R1 barcode subset the corrected R1 barcode belongs to
                if r1_barcode in r1_barcode_subset_dict:
                    subset = r1_barcode_subset_dict[r1_barcode]

                    # Process Read 1
                    r1_seq = read_1.query_sequence
                    r1_qual = read_1.qual
                    if read_1.is_reverse:
                        r1_seq = reverse_complement(r1_seq)
                        r1_qual = r1_qual[::-1]
                    read_1_writers[subset].write(read_1.query_name, r1_seq, r1_qual)

                    # Process Read 2 (construct sequence and quality)
                    r2_seq = read_2.query_sequence + barcode_tag
                    r2_qual = read_2.qual + quality_tag

                    if read_2.is_reverse:
                        r2_seq = reverse_complement(r2_seq)
                        r2_qual = r2_qual[::-1]

                    read_2_writers[subset].write(read_2.query_name, r2_seq, r2_qual)
                    
                # increment QC counter
                if r1_barcode:
                    if exact:
                        exact_match += 1
                    else:
                        nonexact_match += 1
                elif "G"*8 in r1_barcode_window:
                    poly_g_barcode += 1
                else:
                    nonmatch += 1

    # Close all writers
    for writer in read_1_writers.values():
        writer.close()
    for writer in read_2_writers.values():
        writer.close()

    # write QC stats
    with open(f"{prefix}_R1_barcode_qc.txt", "w") as f:
        f.write("%s\t%s\t%s\t%s\t%s\n" % (prefix, exact_match, nonexact_match, nonmatch, poly_g_barcode))

def main():
    args = parse_arguments()
    bam_file = getattr(args, "bam_file")
    pkr = getattr(args, "pkr")
    prefix = getattr(args, "prefix")
    r1_barcode_set_file = getattr(args, "r1_barcode_set_file")
    r2_barcode_file = getattr(args, "r2_barcode_file")
    r3_barcode_file = getattr(args, "r3_barcode_file")
    
    # create dictionary mapping R1 barcodes to their R1 barcode subsets
    r1_barcode_subset_dict = create_barcode_subset_dict(r1_barcode_set_file)
    # create dictionaries of R1 barcode exact matches and 1bp mismatches
    r1_barcode_matches = create_barcode_matching_dict(r1_barcode_subset_dict.keys())
    
    # read R2 and R3 barcode files; if not passed in, use R1 barcodes
    if r2_barcode_file:
        with open(r2_barcode_file) as f:
            r2_barcodes = [barcode for barcode in f.read().rstrip().split()]
    else:
        r2_barcodes = r1_barcode_subset_dict.keys()
        
    if r3_barcode_file:
        with open(r3_barcode_file) as f:
            r3_barcodes = [barcode for barcode in f.read().rstrip().split()]
    else:
        r3_barcodes = r1_barcode_subset_dict.keys()
    
    # reads are written into one R1 FASTQ and one R2 FASTQ per R1 barcode subset;
    # create dictionaries of file pointers associated with each R1 barcode subset,
    # make whitelist for each R1 barcode subset
    read_1_writers = dict()
    read_2_writers = dict()
    for subset in set(r1_barcode_subset_dict.values()):
        fp = BufferedFastqWriter(prefix + "_" + pkr + "_" + subset + "_R1.fastq")
        read_1_writers[subset] = fp
        fp = BufferedFastqWriter(prefix + "_" + pkr + "_" + subset + "_R2.fastq")
        read_2_writers[subset] = fp
        # get possible combinations of R1R2R3, write to whitelist
        r1_barcodes = [k for k,v in r1_barcode_subset_dict.items() if v == subset]
        whitelist_barcodes = [r1+r2+r3 for r1 in r1_barcodes for r2 in r2_barcodes for r3 in r3_barcodes]
        with open(prefix + "_" + pkr + "_" + subset + "_whitelist.txt", "w") as f:
            f.write("\n".join(whitelist_barcodes))
    
    # write reads to FASTQs
    write_fastqs(bam_file, read_1_writers, read_2_writers, r1_barcode_subset_dict, r1_barcode_matches, prefix)
            
if __name__ == "__main__":
    main()
