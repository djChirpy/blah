import sys#, argparse 
import gzip
from Bio import SeqIO, pairwise2
from Bio.pairwise2 import format_alignment
from Bio.Align.Generic import Alignment
from Bio.Align import AlignInfo
from Bio.Alphabet import IUPAC, Gapped


#parser = argparse.ArgumentParser()
#parser.add_argument()

def revComplement(seq):
    '''Returns the reverse complement of a DNA sequence passed to it. 
    Accepts a string as input, and returns a string. Valid bases for 
    input are "ATCGW"'''
    compSeq = ''
    codes = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'W': 'W'}
    for letter in seq:
        compSeq += codes[letter]
    return(compSeq[::-1])



with gzip.open(sys.argv[1], 'rt') as fwdReads, open(sys.argv[2], 'rt') as rvReads:
    fwd = SeqIO.parse(fwdReads, 'fastq')
    rv = SeqIO.parse(rvReads, 'fastq')
    
    read1 = next(fwd).seq
    read2 = revComplement(next(rv).seq)
    
    print(read1)
    print(read2)
    
    alignments = pairwise2.align.globalxx(read1,read2,one_alignment_only=True)
    alignment = Alignment(Gapped(IUPAC.unambiguous_dna, '-'))
    alignment.add_sequence('Alpha',str(alignments[0][0]))
    alignment.add_sequence('Beta',str(alignments[0][1]))
    summaryAlignment = AlignInfo.SummaryInfo(alignment)

    print(summaryAlignment.dumb_consensus())
    
    print(alignment[0][0])
    print(alignment[0][1])
