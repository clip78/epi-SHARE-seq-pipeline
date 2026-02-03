version 1.0

import "subwf-preprocess.wdl" as subwf

workflow demultiplex_fastq {
    input {
        Array[File] read1
        Array[File] read2
        Array[File] index1
        Array[File] index2
        String rna_primers
        String atac_primers
        File metaCsv
        String dockerImage = "us.gcr.io/buenrostro-share-seq/share_task_preprocess"
    }

    # Process each lane
    scatter (i in range(length(read1))) {
        # 1. Convert FASTQs (R1, R2, I1, I2) to an unmapped BAM with tags.
        # This simulates the output of BasecallsToBams (which starts from BCLs).
        call FastqToBam {
            input:
                r1 = read1[i],
                r2 = read2[i],
                i1 = index1[i],
                i2 = index2[i],
                lane = i + 1,
                dockerImage = dockerImage
        }

        # 2. Look up metadata (PKR, library, etc.) based on the lane/barcode info.
        call subwf.BamLookUp {
            input:
                bam = FastqToBam.out_bam,
                metaCsv = metaCsv,
                bucket = ""
        }

        # 3. Demultiplex/Split the BAM into raw FASTQs per library.
        call subwf.BamToRawFastq {
            input:
                bam = FastqToBam.out_bam,
                pkrId = BamLookUp.pkrId,
                library = BamLookUp.library,
                sampleType = BamLookUp.sampleType,
                genome = BamLookUp.genome,
                notes = BamLookUp.notes,
                R1barcodeSet = BamLookUp.R1barcodeSet,
                dockerImage = dockerImage
        }
    }
}

task FastqToBam {
    input {
        File r1
        File r2
        File i1
        File i2
        Int lane
        String dockerImage
    }

    command <<<
        set -e
        # Use python to merge the FASTQs into a SAM file with RX/QX tags
        python3 <<CODE
import gzip
import sys

def open_file(f):
    if f.endswith('.gz'):
        return gzip.open(f, 'rt')
    else:
        return open(f, 'r')

r1_path = "~{r1}"
r2_path = "~{r2}"
i1_path = "~{i1}"
i2_path = "~{i2}"
lane = "~{lane}"
output_sam = "lane" + lane + "_L" + lane + ".sam"

with open_file(r1_path) as f1, open_file(r2_path) as f2, open_file(i1_path) as i1, open_file(i2_path) as i2, open(output_sam, 'w') as out:
    # Write SAM header
    out.write("@HD\tVN:1.6\tSO:unsorted\n")
    out.write(f"@RG\tID:Lane{lane}\tSM:Sample\tLB:Library\tPL:ILLUMINA\n")

    while True:
        # Read 4 lines from each FASTQ
        n1 = f1.readline().strip()
        s1 = f1.readline().strip()
        p1 = f1.readline().strip()
        q1 = f1.readline().strip()

        n2 = f2.readline().strip()
        s2 = f2.readline().strip()
        p2 = f2.readline().strip()
        q2 = f2.readline().strip()

        ni1 = i1.readline().strip()
        si1 = i1.readline().strip()
        pi1 = i1.readline().strip()
        qi1 = i1.readline().strip()

        ni2 = i2.readline().strip()
        si2 = i2.readline().strip()
        pi2 = i2.readline().strip()
        qi2 = i2.readline().strip()

        if not n1: break

        # Construct RX/QX tags (Concatenate I1 and I2)
        rx = si1 + si2
        qx = qi1 + qi2

        # Write Read 1 record
        read_name = n1.split()[0][1:]
        out.write(f"{read_name}\t77\t*\t0\t0\t*\t*\t0\t0\t{s1}\t{q1}\tRG:Z:Lane{lane}\tRX:Z:{rx}\tQX:Z:{qx}\n")

        # Write Read 2 record
        out.write(f"{read_name}\t141\t*\t0\t0\t*\t*\t0\t0\t{s2}\t{q2}\tRG:Z:Lane{lane}\tRX:Z:{rx}\tQX:Z:{qx}\n")

CODE

        # Convert SAM to BAM
        samtools view -bS "lane~{lane}_L~{lane}.sam" > "lane~{lane}_L~{lane}.bam"
        rm "lane~{lane}_L~{lane}.sam"
    >>>

    output {
        File out_bam = "lane~{lane}_L~{lane}.bam"
    }

    runtime {
        docker: dockerImage
        memory: "4 GB"
    }
}
