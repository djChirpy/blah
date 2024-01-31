# Pulls individual SNPs from input data files corresponding to sets of
# input regions. Writes two files; a list of SNPs with ref/alt, etc, and
# second, a bed-formated list.


import sys
import csv

with open(sys.argv[1]) as assoc, open(sys.argv[2]) as regions, open(sys.argv[3] + ".out", 'w') as outFile:

    assocReader = csv.DictReader(
        assoc,
        fieldnames=(
            'CHR',
            'SNP',
            'BP',
            'A1',
            'F_A',
            'F_U',
            'A2',
            'CHISQ',
            'P',
            'OR'), delimiter = '\t')
    
    regionReader = csv.DictReader(
        regions,
        fieldnames=(
            'CHR',
            'StartPos',
            'EndPos',
            'name',
            'V5',
            'Strand'), delimiter = '\t')
    
    snpWriter = csv.DictWriter(
        outFile,
        fieldnames=(
            'CHR',
            'POS',
            'ID',
            'REF',
            'ALT',
            'CHISQ',
            'P',
            'OR'),
        delimiter='\t')
   

    snpWriter.writeheader()

    assocHeader = next(assocReader)
    print(assocHeader)
    
    snp = next(assocReader)

    for region in regionReader:
        print(str(region) + ' ' + snp['CHR'] + ' ' + snp['BP'])
        while int(snp['CHR']) < int(region['CHR'][3:]):
            snp = next(assocReader)
        while int(snp['CHR']) == int(region['CHR'][3:]) and int(snp['BP']) < int(region['StartPos']):
            snp = next(assocReader)
        while int(snp['CHR']) == int(region['CHR'][3:]) and int(snp['BP']) >= int(region['StartPos']) and int(snp['BP']) <= int(region['EndPos']):
            outRow = {
                'CHR': snp['CHR'],
                'POS': snp['BP'],
                'ID': snp['SNP'],
                'REF': snp['A1'],
                'ALT': snp['A2'],
                'CHISQ': snp['CHISQ'],
                'P': snp['P'],
                'OR': snp['OR']
                }
            
            snpWriter.writerow(outRow)
           
            snp = next(assocReader)
      
