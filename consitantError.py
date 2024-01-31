import sys
import re


oligoSummary={}

def summarizeBarcode(oligos, barcode):
    if len(oligos) == 1:
        for oligo, reads in oligos.items():
            
            try:
                oligoSummary[oligo][barcode] = [reads[0],reads[1]]
            except KeyError:
                oligoSummary[oligo] = {barcode: [reads[0],reads[1]]}
                
with open(sys.argv[1], 'r') as inFile, open(sys.argv[1][:-3] + 'error_summary.out', 'w') as outFile::
    
    barcode = ''
    oligos = {}
    
    for line in inFile:
        fields = re.split('\s+', line)
        oligo = fields[0]
        
        if oligo[:3] == 'chr':
            summarizeBarcode(oligos, barcode)
            oligos = {}
            barcode = oligo
            continue
        else:
            
            perfectReads = fields[2]
            imperfectReads = fields[4]
            readCounts = [perfectReads, imperfectReads]
            
        
            oligos[oligo] = readCounts
        
    for oligo, barcodes in oligoSummary.items():
        numBarcodes = 0
        perfect = 0
        imperfect = 0      
        #print(oligo)
        for barcode, reads in barcodes.items():
            numBarcodes += 1
            uniquFile.write(oligo + '\t' + barcode + '\t' + str(reads[0]) + '\t' + str(reads[1] + '\n'))
            perfect += int(reads[0])
            imperfect += int(reads[1])
        
        summary.write(oligo + '\t' + str(numBarcodes) + '\t' + str(perfect) + '\t' + str(imperfect) + '\n')
