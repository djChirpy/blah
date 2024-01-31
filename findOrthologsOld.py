# usage python findOrthologs.py -i TWIST.annotated.bed -o test -de /seq/vgb/swofford/ref/canFam3EnsemblGeneAnnoUCSC021819.gz -or /seq/vgb/swofford/ref/HumanDogOrthologsHumanInclusive021220.txt -gd /seq/vgb/swofford/ref/Homo_sapiens.GRCh37.87.chr.gtf.gz

import re
import csv
import argparse
from csv import DictReader
import gzip
import os



parser = argparse.ArgumentParser()
parser.add_argument('-i', '--inFile') #TWIST human design file
parser.add_argument('-o', '--outFile')
parser.add_argument('-or', '--refFile') #ensembl ortholog annotation from biomart
parser.add_argument('-gd', '--geneAnno') #ensembl gene annotation hg19
parser.add_argument('-d', '--dogEns')

args = parser.parse_args()

'''Object type for storing gene information from dog Ensembl UCSC track'''
class gene:
    def __init__(self, geneID, TSS, TES, cdsStart, cdsEnd, exonCount, 
                 exonStarts, exonEnds, exonIntervals, strand, transcriptID, 
                 cdsStarts, cdsEnds, cdsIntervals, UTR5, UTR3, chrom):
        self.geneID = geneID
        self.TSS = TSS
        self.TES = TES
        self.cdsStart = cdsStart
        self.cdsEnd = cdsEnd
        self.exonCount = exonCount
        self.exonStarts = exonStarts
        self.exonEnds = exonEnds
        self.exonIntervals = exonIntervals
        self.strand = strand
        self.transcriptID = transcriptID
        self.cdsStarts = cdsStarts
        self.cdsEnds = cdsEnds
        self.cdsIntervals = cdsIntervals
        self.UTR5 = UTR5
        self.UTR3 = UTR3
        self.chrom = chrom
        

''' Object type for storing human/dog ortholog information. '''
class ortholog:
    def __init__(self, dogGeneName, humanGeneName, dogGeneID, geneID, 
                 dogChr, humanChr, dogStart, dogEnd, humanStart, 
                 humanEnd, conservationScore, orthologConfidence, orthologType):
        self.dogGeneName = dogGeneName
        self.humanGeneName = humanGeneName
        self.dogGeneID = dogGeneID
        self.geneID = geneID
        self.dogChr = dogChr
        self.humanChr = humanChr
        self.dogStart = dogStart
        self.dogEnd = dogEnd
        self.humanStart = humanStart
        self.humanEnd = humanEnd
        self.conservationScore = conservationScore
        self.orthologConfidence = orthologConfidence
        self.orthologType = orthologType
        
'''Builds a dictionary of gene names/ids and their associated ortholog objects
to facilitate lookups.
'''
def buildRef(refFile):
    with open(refFile, 'r') as refFile: 
        ref = {}
        refReader =  csv.DictReader(refFile, delimiter='\t')
        for line in refReader:
            dogGeneName = line['Dog gene name']
            humanGeneName = line['Gene name']
            dogGeneID = line['Dog gene stable ID']
            geneID = line['Gene stable ID']
            dogChr = line['Dog chromosome/scaffold name']
            humanChr = line['Chromosome/scaffold name']
            dogStart = line['Dog chromosome/scaffold start (bp)']
            dogEnd = line['Dog chromosome/scaffold end (bp)']
            humanStart = line['Gene start (bp)']
            humanEnd = line['Gene end (bp)']
            conservationScore = line['Dog Gene-order conservation score']
            orthologConfidence = line['Dog orthology confidence [0 low, 1 high]']
            orthologType = line['Dog homology type']
            

            try:
                ref[geneID].append(ortholog(dogGeneName=dogGeneName, humanGeneName=humanGeneName, dogGeneID=dogGeneID, 
                                       geneID=geneID, dogChr=dogChr, humanChr=humanChr, dogStart=dogStart, 
                                       dogEnd=dogEnd, humanStart=humanStart, humanEnd=humanEnd, 
                                       conservationScore=conservationScore, orthologConfidence=orthologConfidence, 
                                       orthologType=orthologType))
            except KeyError:
                ref[geneID] = [ortholog(dogGeneName=dogGeneName, humanGeneName=humanGeneName, dogGeneID=dogGeneID, 
                                       geneID=geneID, dogChr=dogChr, humanChr=humanChr, dogStart=dogStart, 
                                       dogEnd=dogEnd, humanStart=humanStart, humanEnd=humanEnd, 
                                       conservationScore=conservationScore, orthologConfidence=orthologConfidence, 
                                       orthologType=orthologType)]
            try:      
                ref[humanGeneName].append(ortholog(dogGeneName=dogGeneName, humanGeneName=humanGeneName, dogGeneID=dogGeneID, 
                                       geneID=geneID, dogChr=dogChr, humanChr=humanChr, dogStart=dogStart, 
                                       dogEnd=dogEnd, humanStart=humanStart, humanEnd=humanEnd, 
                                       conservationScore=conservationScore, orthologConfidence=orthologConfidence, 
                                       orthologType=orthologType))
            except KeyError:
                ref[humanGeneName] = [ortholog(dogGeneName=dogGeneName, humanGeneName=humanGeneName, dogGeneID=dogGeneID, 
                                       geneID=geneID, dogChr=dogChr, humanChr=humanChr, dogStart=dogStart, 
                                       dogEnd=dogEnd, humanStart=humanStart, humanEnd=humanEnd, 
                                       conservationScore=conservationScore, orthologConfidence=orthologConfidence, 
                                       orthologType=orthologType)]
                
            
        print('Finished building Reference. ' + str(int(len(ref.keys())/2)) + ' orthologs recorded.')
        return ref

'''Builds a lookup table for searching hg19 for gene ids associated with gene names found in TWIST 
but not in hg38
'''
def buildLookup(gtf):
    hg19Lookup = {}
    with gzip.open(gtf, 'rt') as gtfFile:
        #skip header lines
        for line in gtfFile:
            if line[:1] != '#':
                break
            
        gtfReader = DictReader(gtfFile, delimiter='\t', fieldnames = ['seqname', 'source', 
                                                                      'feature', 'start', 
                                                                      'end', 'score',
                                                                      'strand','frame', 
                                                                      'attribute'])
        for line in gtfReader:
            biotype = line['attribute'].split('biotype')[1].split(';')[0][2:-1]
            
            gene_name = line['attribute'].split('gene_name')[1].split(';')[0][2:-1]#.split("\"")[0]
            gene_id = line['attribute'].split('gene_id')[1].split(';')[0][2:-1]#.split("\"")[0]
            hg19Lookup[gene_name] = gene_id
            
    print('Finished building hg19 lookup. ' + str(len(hg19Lookup.keys())) + ' gene/id pairs.')
    return(hg19Lookup)    
        
def buildDogGeneReference(dogEns):
    dogGenes = {}
    with gzip.open(dogEns, 'rt') as dogRef:
        reader = csv.DictReader(dogRef, delimiter = '\t')
        for line in reader:
            geneID = line['name2'].split('.')[0]
            TSS = int(line['txStart'])
            TES = int(line['txEnd'])
            cdsStart = int(line['cdsStart'])
            cdsEnd = int(line['cdsEnd'])
            exonCount = int(line['exonCount'])
            chrom = line['chrom']
            exonStarts = [int(x) for x in line['exonStarts'].split(',')[:-1]]
            exonEnds = [int(x) for x in line['exonEnds'].split(',')[:-1]]
            exonIntervals = []
            for start, end in zip(exonStarts, exonEnds):
                exonIntervals.append([start, end])
            strand = line['strand']
            transcriptID = line['name'].split('.')[0]
            cdsStarts = [int(x) for x in line['exonStarts'].split(',')[:-1]]
            cdsStarts[0] = cdsStart
            cdsEnds = [int(x) for x in line['exonEnds'].split(',')[:-1]]
            cdsEnds[exonCount - 1] = cdsEnd
            cdsIntervals = []
            for start, end in zip(cdsStarts, cdsEnds):
                cdsIntervals.append([start, end])
            
                
            '''    
            for range in tempInterval:
                print(range)
            print('break')
            for range in tempInterval:
                print(range)
            print('break2')
                
            
            for start, end in cdsIntervals:
                print(str(start) + '-' + str(end), end = ',')
            print('')
            for start, end in exonIntervals:
                print(str(start) + '-' + str(end), end = ',')  
            print('\n')
            '''
            if strand == '+':
                UTR5 = (TSS, cdsStart)
                UTR3 = (cdsEnd, TES)
            else:
                UTR5 = (cdsEnd, TES)
                UTR3 = (TSS, cdsStart)
            #print(UTR5)
            #print(UTR3)
            dogGenes[geneID] = gene(geneID = geneID, TSS = TSS, TES = TES, 
                                    cdsStart = cdsStart, cdsEnd = cdsEnd, 
                                    exonCount = exonCount, 
                                    exonStarts = exonStarts, 
                                    exonEnds = exonEnds, 
                                    exonIntervals = exonIntervals, 
                                    strand = strand, 
                                    transcriptID = transcriptID, 
                                    cdsStarts = cdsStarts, 
                                    cdsEnds = cdsEnds, 
                                    cdsIntervals = cdsIntervals, 
                                    UTR5 = UTR5, UTR3 = UTR3, chrom = chrom)
    return(dogGenes)
            
                           
print('Building dog gene details')
dogEns = buildDogGeneReference(args.dogEns)
print(str(len(dogEns)) + ' genes recorded from ' + args.dogEns)
print('Building ortholog reference from:' + args.refFile)
reference = buildRef(args.refFile)
print('Building hg19 gene/id lookup from:' + args.geneAnno)
hg19Lookup = buildLookup(args.geneAnno)

tempID = reference['COL11A1'][0].dogGeneID
'''for entry in reference['COL11A1']:
    print(entry)
    for interval in dogEns[entry].cdsIntervals:
        print(dogEns[entry].chr)
'''
print(dogEns[tempID].UTR5)
for entry in dogEns[tempID].cdsIntervals:
    print(dogEns[tempID].chrom + '\t' + str(entry[0]) + '\t' + str(entry[1]))

    
with open(args.inFile, 'r') as inFile, \
    open(args.outFile + 'primary_one2one.tmpout', 'w') as primary_one2one, \
    open(args.outFile + 'primary_otherOrtho.tmpout', 'w') as primary_otherOrtho, \
    open(args.outFile + 'primary_noOrtho.tmpout', 'w') as primary_noOrtho, \
    open(args.outFile + 'secondary_one2one.tmpout', 'w') as secondary_one2one, \
    open(args.outFile + 'secondary_otherOrtho.tmpout', 'w') as secondary_otherOrtho, \
    open(args.outFile + 'secondary_noOrtho.tmpout', 'w') as secondary_noOrtho, \
    open(args.outFile + 'unknownGenes.tmpout', 'w') as unknown:
    
    reader = csv.DictReader(inFile, delimiter = '\t', fieldnames=['Chr','start','end','description'])
    fieldnames = ['Chr', 'start', 'end', 'name', 'confidence(o/1)', 'geneName', 'geneID', 'dogGeneName', 'dogGeneID', 'seqType']
    fieldnamesBase = ['Chr', 'start', 'end', 'name', 'geneName']
    pone2one = csv.DictWriter(primary_one2one, fieldnames = fieldnames, delimiter = '\t')
    #pone2one.writeheader()
    #primary_one2one.write('Chr\tstart\tend\tname\tconfidence(0/1)\tgeneName\tgeneID\tdogGeneName\tdogGeneID\tseqType\n')
    pone2many = csv.DictWriter(primary_otherOrtho, fieldnames = fieldnames, delimiter = '\t')
    #pone2many.writeheader()
    #primary_otherOrtho.write('Chr\tstart\tend\tname\tconfidence(0/1)\tgeneName\tgeneID\tdogGeneName\tdogGeneID\tseqType\n')
    pnoOrtho = csv.DictWriter(primary_noOrtho, fieldnames = fieldnamesBase, delimiter = '\t')
    #pnoOrtho.writeheader()
    #primary_noOrtho.write('Chr\tstart\tend\tname\tgeneName\n')
    sone2one = csv.DictWriter(secondary_one2one, fieldnames = fieldnames, delimiter = '\t')
    #sone2one.writeheader()
    #secondary_one2one.write('Chr\tstart\tend\tname\tconfidence(0/1)\tgeneName\tgeneID\tdogGeneName\tdogGeneID\tseqType\n')
    sone2many = csv.DictWriter(secondary_otherOrtho, fieldnames = fieldnames, delimiter = '\t')
    #sone2many.writeheader()
    #secondary_otherOrtho.write('Chr\tstart\tend\tname\tconfidence(0/1)\tgeneName\tgeneID\tdogGeneName\tdogGeneID\tseqType\n')
    snoOrtho = csv.DictWriter(secondary_noOrtho, fieldnames = fieldnamesBase, delimiter = '\t')
    #snoOrtho.writeheader()
    
    unk = csv.DictWriter(unknown, fieldnames = ['Chr', 'start', 'end', 'name', 'twistName', 'hg19Name'], delimiter = '\t')
    #unk.writeheader()
    
    
    #secondary_one2one.write('Chr\tstart\tend\tname\tconfidence(0/1)\tgeneName\tgeneID\tdogGeneName\tdogGeneID\tseqType\n')
    #secondary_otherOrtho.write('Chr\tstart\tend\tname\tconfidence(0/1)\tgeneName\tgeneID\tdogGeneName\tdogGeneID\tseqType\n')
    #secondary_noOrtho.write('Chr\tstart\tend\tname\tgeneName\n')
    
    for line in reader:
        targets = [x.strip() for x in line['description'].split(',')] 
          
        primary = targets[0]
        chr = line['Chr']
        start = line['start']
        end = line['end']
        try:
            secondary = targets[1:]
        except IndexError:
            secondary = []
        
        fields = primary.split(':')
        gene = fields[0]
        if gene not in reference.keys():
            try:
                gene = hg19Lookup[gene]
            except KeyError:
                pass
        try:
            seqType = fields[3]
        except IndexError:
            pass
        
        
        if gene in reference.keys():

            for entry in reference[gene]:
                if entry.orthologType == 'ortholog_one2one':
                    pone2one.writerow({'Chr':chr,
                                       'start':start, 
                                       'end':end, 
                                       'name':primary, 
                                       'confidence(o/1)':entry.orthologConfidence, 
                                       'geneName':entry.humanGeneName, 
                                       'geneID':entry.geneID, 
                                       'dogGeneName':entry.dogGeneName, 
                                       'dogGeneID':entry.dogGeneID, 
                                       'seqType':seqType})
                    '''
                    primary_one2one.write(chr + '\t' + start + '\t' + end + '\t' + primary + '\t' 
                                          + entry.orthologConfidence + '\t' + entry.humanGeneName 
                                          + '\t' + entry.geneID + '\t' + entry.dogGeneName 
                                          + '\t' + entry.dogGeneID + '\t' + seqType + '\n')
                    '''
                elif entry.orthologType == 'ortholog_one2many':
                    for orthotype in entry.orthologType:
                        pone2many.writerow({'Chr':chr,
                                       'start':start, 
                                       'end':end, 
                                       'name':primary, 
                                       'confidence(o/1)':entry.orthologConfidence, 
                                       'geneName':entry.humanGeneName, 
                                       'geneID':entry.geneID, 
                                       'dogGeneName':entry.dogGeneName, 
                                       'dogGeneID':entry.dogGeneID, 
                                       'seqType':seqType})
                        '''
                        primary_otherOrtho.write(chr + '\t' + start + '\t' + end + '\t' + primary + '\t' 
                                                 + entry.orthologConfidence + '\t' + entry.humanGeneName + '\t' 
                                                 + entry.geneID + '\t' + entry.dogGeneName + '\t' 
                                                 + entry.dogGeneID + '\t' + entry.orthologType 
                                                 + '\t' + seqType + '\n')
                        '''
                else:
                    pnoOrtho.writerow({'Chr':chr, 
                                       'start':start, 
                                       'end':end, 
                                       'name':primary, 
                                       'geneName':primary.split(':')[0].strip()})
                    #primary_noOrtho.write(chr + '\t' + start + '\t' + end + '\t' + primary + '\t' + primary.split(':')[0] + '\n')
        
        else:
            unk.writerow({'Chr':chr, 
                          'start':start, 
                          'end':end, 
                          'name':primary,
                          'twistName':primary.split(':')[0].strip(), 
                          'hg19Name':gene})
            
            #unknown.write(gene + '\t' + primary.split(':')[0] + '\t' + chr + '\t' + start + '\t' + end + '\n')
        
        for target in secondary:
            gene = target.split(':')[0].strip()
            if gene not in reference.keys():
                try:
                    gene = hg19Lookup[gene]
                    
                except KeyError:
                    pass
            
            if gene in reference.keys():
                for entry in reference[gene]:
                    if entry.orthologType == 'ortholog_one2one':
                        
                        sone2one.writerow({'Chr':chr,
                                       'start':start, 
                                       'end':end, 
                                       'name':target, 
                                       'confidence(o/1)':entry.orthologConfidence, 
                                       'geneName':entry.humanGeneName, 
                                       'geneID':entry.geneID, 
                                       'dogGeneName':entry.dogGeneName, 
                                       'dogGeneID':entry.dogGeneID, 
                                       'seqType':seqType})
                        
                        '''
                        secondary_one2one.write(chr + '\t' + start + '\t' + end + '\t' + target + '\t' 
                                              + entry.orthologConfidence + '\t' + entry.humanGeneName 
                                              + '\t' + entry.geneID + '\t' + entry.dogGeneName 
                                              + '\t' + entry.dogGeneID + '\t' + seqType + '\n')
                        '''
                    elif entry.orthologType == 'ortholog_one2many':
                        for orthotype in entry.orthologType:
                            sone2many.writerow({'Chr':chr,
                                       'start':start, 
                                       'end':end, 
                                       'name':target, 
                                       'confidence(o/1)':entry.orthologConfidence, 
                                       'geneName':entry.humanGeneName, 
                                       'geneID':entry.geneID, 
                                       'dogGeneName':entry.dogGeneName, 
                                       'dogGeneID':entry.dogGeneID, 
                                       'seqType':seqType})
                            '''
                            secondary_otherOrtho.write(chr + '\t' + start + '\t' + end + '\t' + target + '\t' 
                                                     + entry.orthologConfidence + '\t' + entry.humanGeneName + '\t' 
                                                     + entry.geneID + '\t' + entry.dogGeneName + '\t' 
                                                     + entry.dogGeneID + '\t' + entry.orthologType 
                                                     + '\t' + seqType + '\n')
                            '''
                    else:
                        snoOrtho.writerow({'Chr':chr, 
                                       'start':start, 
                                       'end':end, 
                                       'name':target, 
                                       'geneName':target.split(':')[0].strip()})
                        
                        #secondary_noOrtho.write(chr + '\t' + start + '\t' + end + '\t' + target + '\t' + target.split(':')[0] + '\n')
            
            else:
                unk.writerow({'Chr':chr, 
                          'start':start, 
                          'end':end, 
                          'name':target,
                          'twistName':target.split(':')[0].strip(), 
                          'hg19Name':gene})
                #unknown.write(chr + '\t' + start + '\t' + end + '\t' + target + '\t' + target.split(':')[0] + '\t' + gene + '\n')

#Using bash to sort and remove duplicate entries for each file    
os.system('for file in *.tmpout; do awk \'NR == 1; NR > 1 {print $0 | "sort -n"}\' $file | uniq > $(basename $file .tmpout).out; rm $file; done;')

os.system('grep \':CDS\' testprimary_one2one.out >primary_one2oneCDS.out')
os.system('grep \':Exon\' testprimary_one2one.out >primary_one2oneExon.out')
os.system('grep \':Intron\' testprimary_one2one.out >primary_one2oneIntron.out')
os.system('grep \':TES\' testprimary_one2one.out >primary_one2oneTES.out')
os.system('grep \':TSS\' testprimary_one2one.out >primary_one2oneTSS.out')
os.system('grep \':CDS\' testsecondary_one2one.out >secondary_one2oneCDS.out')
os.system('grep \':Exon\' testsecondary_one2one.out >secondary_one2oneExon.out')
os.system('grep \':Intron\' testsecondary_one2one.out >secondary_one2oneIntron.out')
os.system('grep \':TES\' testsecondary_one2one.out >secondary_one2oneTES.out')
os.system('grep \':TSS\' testsecondary_one2one.out >secondary_one2oneTSS.out')
os.system('grep \'Intergenic\' testunknownGenes.out > intergenic.out')
os.system('grep -v \'Intergenic\' testunknownGenes.out > unknownGenes.out')

    
