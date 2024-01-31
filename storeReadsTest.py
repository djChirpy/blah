import pysam
import sys
from Bio import pairwise2

counter = 0

with pysam.AlignmentFile(sys.argv[1], 'rb') as bam:
    iter = bam.fetch("chr20_48531169_rs3819809_T")
    barcodes = set()
    reads1 = {}
    reads2 = {}
    reads3 = {}
    for read in iter:

        if read.is_read1 and read.is_proper_pair and read.is_paired:
            name = read.query_name
            #print(name)
            reads1[name] = read
        elif read.is_read2 and read.is_proper_pair and read.is_paired:
            name = read.query_name
            #print(name)
            reads2[name] = read
            
    print(str(len(reads1)))
    print(str(len(reads2)))

    for name, read in reads1.items():
        reads3[name] = [read, reads2[name]]

    print(str(len(reads3)))

    for name, reads in reads3.items():
        print(str(len(reads3)))
        print(str(reads[0]))
        print(str(reads[1]))

        

        #alignments = pairwise2.align.globalxx(
