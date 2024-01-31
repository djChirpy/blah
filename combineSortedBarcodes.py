import sys
import re

def isNumber(num):
    try:
        float(num)
        return True
    except ValueError:
        return False
        

def revComplement(seq):
    '''Returns the reverse complement of a DNA sequence passed to it.
    Accepts a string as input, and returns a string. Valid bases for
    input are "ATCGW"'''
    compSeq = ''
    codes = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'W': 'W', 'N': 'N'}
    for letter in seq:
        compSeq += codes[letter]
    return(compSeq[::-1])

def writeBarcode(barcode, oligo):
    barcodeReads = 0
    for k, v, in oligos.items():
        barcodeReads += v
    outPut.write(barcode + '\t' + str(len(oligos)) + '\t' + str(barcodeReads))
    for k, v in oligos.items():
        outPut.write('\t' + k + '\t' + str(v))
    outPut.write('\n')


with open(sys.argv[1], 'rt') as reads, open(sys.argv[2], 'w') as outPut:
    #fields = re.split('\t', reads.readline().rstrip())
    #barcode = fields[0]
    #seq = fields[1]
    barcodeCount = 0
    barcodeOligos = 0
    storedBarcode = 'first'
    oligos = {}
    
    testCounter = 0
    for line in reads:
        try:       
            fields = re.split('\s+', line)
            barcode = fields[0].strip()
            barcodeCount += int(fields[1])
            barcodeOligos += int(fields[2])
        except ValueError:
            print(line)
            continue
        
        if barcode != storedBarcode:
            writeBarcode(storedBarcode, oligos)
            storedBarcode = barcode
            oligos = {}
            barcodeCount = 0
            barcodeOligos = 0
    
        
        for field in fields[3:]:
            if isNumber(field):
                count = int(field)
                if seq in oligos:
                    oligos[seq] += count    
                else:
                    oligos[seq] = count
                seq = ''
            else:
                seq = field.strip()

