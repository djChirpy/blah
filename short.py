#!/usr/bin/env python

# Derives barcode/oligo associations from MPRA sequencing FASTQs
# Usage: python barcodAssociation.py fwdReadFastQ rvsReadFastQ
# ReferenceOligoFile outputPrefix


import sys
import gzip
import re
from Bio import SeqIO, pairwise2



def revComplement(seq):
    '''Returns the reverse complement of a DNA sequence passed to it.
    Accepts a string as input, and returns a string. Valid bases for
    input are "ATCGW"'''
    compSeq = ''
    codes = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'W': 'W', 'N': 'N'}
    for letter in seq:
        compSeq += codes[letter]
    return(compSeq[::-1])


def buildRef():
    '''Returns a dictionary of ordered oligo variant inserts and names'''
    ref = open(sys.argv[3])
    variants = {}
    for line in ref:
        fields = re.split('\s+', line.rstrip())
        name = str(fields[0]) + '_' + str(fields[1]) + \
            '_' + str(fields[2]) + '_' + str(fields[9])
        var = str(fields[10])
        name = name.replace(',', '^')
        variants[var] = name
    ref.close()
    return variants

oligos = buildRef()
unmatched = {}
barcodes = {}
chunkSize = int(sys.argv[5])

with gzip.open(sys.argv[1], 'rt') as fwdReads, gzip.open(sys.argv[2], 'rt') as rvReads, open(sys.argv[4] + '_matched.out', 'w') as matchedOutFile, open(sys.argv[4] + '_unmatched.out', 'w') as unmatchedOutFile, open(sys.argv[4] + '_barcodes.out', 'w') as barcodeOutFile:

    fwd = SeqIO.parse(fwdReads, 'fastq')
    rv = SeqIO.parse(rvReads, 'fastq')
    eof = False
    good = 0
    bad = 0

    # write headers
    matchedOutFile.write('Name\tCoverage\tBarcodes\n')
    unmatchedOutFile.write('Sequence\tCoverage\tBarcodes\n')
    barcodeOutFile.write(
        'Barcode\tCoverage\tNumberMatchingOligos\tNameOfMatchingOligo\tMatchingOligoSeq\n')

    #for i in range (2):
    while not eof:
        for i in range(chunkSize):
            # read from fwd and reverse fastqs until end
            if eof:
                continue
            try:
                read1 = next(fwd).seq
                read2 = revComplement(next(rv).seq)
            except:
                eof = True
                continue

            # barcode corresponds to the first 19 bases of the fwd read
            barcode = str(read1[:19])
            # trimming barcode and static sequence from fwd read
            read1 = read1[57:]
            # trimming static sequence from rv read
            read2 = read2[:len(read2) - 16]
            # align trimmed reads
            alignments = pairwise2.align.globalxx(
                read1,
                read2,
                one_alignment_only=True)

            consensus = ''

            # Combine aligned reads into single consensus sequence.
            for a in zip(str(alignments[0][0]), str(alignments[0][1])):
                if a[0] == a[1]:
                    consensus = consensus + str(a[0])
                elif a[0] == '-':
                    consensus = consensus + str(a[1])
                elif a[1] == '-':
                    consensus = consensus + str(a[0])
                else:
                    consensus = consensus + 'N'

            # sequence ends up being the reverse complement of oligs as ordered
            consensus = revComplement(consensus)
            
            # Check for seq in reference list
            try:
                seqName = oligos[consensus]
                good += 1
            except KeyError:
                seqName = consensus
                bad += 1
            try:
                barcodes[barcode][0] += 1
                try:
                    idx = barcodes[barcode].index(seqName)
                    barcodes[barcode][idx + 1] += 1
                except ValueError:
                    barcodes[barcode][1] += 1
                    barcodes[barcode].append(seqName)
                    barcodes[barcode].append(1)
            
            except KeyError:
                barcodes[barcode] = [1,1,seqName,1]        
                    
                
        for k, v in barcodes.items():
            barcodeOutFile.write(k + '\t')
            for variant in v:
                barcodeOutFile.write(str(variant) + '\t')
            barcodeOutFile.write('\n')
        oligos = buildRef()
        unmatched = {}
        barcodes = {}

    print('Good: ' + str(good))
    print('Bad: ' + str(bad))
