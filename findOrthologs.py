# usage python findOrthologs.py -i TWIST.annotated.bed -db testdb -or
# /seq/vgb/swofford/ref/HumanDogOrthologsHumanInclusive021220.txt -gd
# /seq/vgb/swofford/ref/Homo_sapiens.GRCh37.87.chr.gtf.gz
import re
import csv
import argparse
from csv import DictReader
import gzip
import os
import gffutils

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--inFile')  # TWIST human design file
# ensembl ortholog annotation from biomart
parser.add_argument('-or', '--refFile')
parser.add_argument('-gd', '--geneAnno')  # ensembl gene annotation hg19
parser.add_argument('-db', '--dogDb')

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
    def __init__(
            self,
            dogGeneName,
            humanGeneName,
            dogGeneID,
            geneID,
            dogChr,
            humanChr,
            dogStart,
            dogEnd,
            humanStart,
            humanEnd,
            conservationScore,
            orthologConfidence,
            orthologType):
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


def isNumber(x):
    try:
        float(x)
        return True
    except ValueError:
        return False


def fixChr(inchr):
    if isNumber(inchr):
        return('chr' + str(inchr))
    elif 'chr' in inchr:
        return(inchr)
    elif inchr == 'MT':
        return('chrM')
    elif inchr == 'X':
        return('chrX')
    else:
        return('chrUn_' + str(inchr).split('.')[0])


'''Builds a dictionary of gene names/ids and their associated ortholog objects
to facilitate lookups.
'''


def buildRef(refFile):
    with open(refFile, 'r') as refFile:
        ref = {}
        refReader = csv.DictReader(refFile, delimiter='\t')
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
                ref[geneID].append(
                    ortholog(
                        dogGeneName=dogGeneName,
                        humanGeneName=humanGeneName,
                        dogGeneID=dogGeneID,
                        geneID=geneID,
                        dogChr=dogChr,
                        humanChr=humanChr,
                        dogStart=dogStart,
                        dogEnd=dogEnd,
                        humanStart=humanStart,
                        humanEnd=humanEnd,
                        conservationScore=conservationScore,
                        orthologConfidence=orthologConfidence,
                        orthologType=orthologType))
            except KeyError:
                ref[geneID] = [
                    ortholog(
                        dogGeneName=dogGeneName,
                        humanGeneName=humanGeneName,
                        dogGeneID=dogGeneID,
                        geneID=geneID,
                        dogChr=dogChr,
                        humanChr=humanChr,
                        dogStart=dogStart,
                        dogEnd=dogEnd,
                        humanStart=humanStart,
                        humanEnd=humanEnd,
                        conservationScore=conservationScore,
                        orthologConfidence=orthologConfidence,
                        orthologType=orthologType)]
            try:
                ref[humanGeneName].append(
                    ortholog(
                        dogGeneName=dogGeneName,
                        humanGeneName=humanGeneName,
                        dogGeneID=dogGeneID,
                        geneID=geneID,
                        dogChr=dogChr,
                        humanChr=humanChr,
                        dogStart=dogStart,
                        dogEnd=dogEnd,
                        humanStart=humanStart,
                        humanEnd=humanEnd,
                        conservationScore=conservationScore,
                        orthologConfidence=orthologConfidence,
                        orthologType=orthologType))
            except KeyError:
                ref[humanGeneName] = [
                    ortholog(
                        dogGeneName=dogGeneName,
                        humanGeneName=humanGeneName,
                        dogGeneID=dogGeneID,
                        geneID=geneID,
                        dogChr=dogChr,
                        humanChr=humanChr,
                        dogStart=dogStart,
                        dogEnd=dogEnd,
                        humanStart=humanStart,
                        humanEnd=humanEnd,
                        conservationScore=conservationScore,
                        orthologConfidence=orthologConfidence,
                        orthologType=orthologType)]

        print('Finished building Reference. ' +
              str(int(len(ref.keys()) / 2)) + ' orthologs recorded.')
        return ref


'''Builds a lookup table for searching hg19 for gene ids associated with gene names found in TWIST
but not in hg38
'''


def buildLookup(gtf):
    hg19Lookup = {}
    with gzip.open(gtf, 'rt') as gtfFile:
        # skip header lines
        for line in gtfFile:
            if line[:1] != '#':
                break

        gtfReader = DictReader(
            gtfFile,
            delimiter='\t',
            fieldnames=[
                'seqname',
                'source',
                'feature',
                'start',
                'end',
                'score',
                'strand',
                'frame',
                'attribute'])
        for line in gtfReader:
            biotype = line['attribute'].split('biotype')[1].split(';')[0][2:-1]

            gene_name = line['attribute'].split('gene_name')[1].split(';')[
                0][2:-1]  # .split("\"")[0]
            gene_id = line['attribute'].split('gene_id')[1].split(';')[
                0][2:-1]  # .split("\"")[0]
            hg19Lookup[gene_name] = gene_id

    print('Finished building hg19 lookup. ' +
          str(len(hg19Lookup.keys())) + ' gene/id pairs.')
    return(hg19Lookup)


print('Building ortholog reference from:' + args.refFile)
reference = buildRef(args.refFile)
print('Building hg19 gene/id lookup from:' + args.geneAnno)
hg19Lookup = buildLookup(args.geneAnno)

primeo2o = {}
primeo2m = {}
secondaryo2o = {}
secondaryo2m = {}
noOrtho = {}

with open(args.inFile, 'r') as inFile, \
        open('primary_one2one.tmpout', 'w') as primary_one2one, \
        open('primary_otherOrtho.tmpout', 'w') as primary_otherOrtho, \
        open('primary_noOrtho.tmpout', 'w') as primary_noOrtho, \
        open('secondary_one2one.tmpout', 'w') as secondary_one2one, \
        open('secondary_otherOrtho.tmpout', 'w') as secondary_otherOrtho, \
        open('secondary_noOrtho.tmpout', 'w') as secondary_noOrtho, \
        open('unknownGenes.tmpout', 'w') as unknown:

    reader = csv.DictReader(
        inFile, delimiter='\t', fieldnames=[
            'Chr', 'start', 'end', 'description'])
    fieldnames = [
        'Chr',
        'start',
        'end',
        'name',
        'confidence(o/1)',
        'geneName',
        'geneID',
        'dogGeneName',
        'dogGeneID',
        'seqType']
    fieldnamesBase = ['Chr', 'start', 'end', 'name', 'geneName']
    pone2one = csv.DictWriter(
        primary_one2one,
        fieldnames=fieldnames,
        delimiter='\t')

    pone2many = csv.DictWriter(
        primary_otherOrtho,
        fieldnames=fieldnames,
        delimiter='\t')

    pnoOrtho = csv.DictWriter(
        primary_noOrtho,
        fieldnames=fieldnamesBase,
        delimiter='\t')

    sone2one = csv.DictWriter(
        secondary_one2one,
        fieldnames=fieldnames,
        delimiter='\t')

    sone2many = csv.DictWriter(
        secondary_otherOrtho,
        fieldnames=fieldnames,
        delimiter='\t')

    snoOrtho = csv.DictWriter(
        secondary_noOrtho,
        fieldnames=fieldnamesBase,
        delimiter='\t')

    unk = csv.DictWriter(
        unknown,
        fieldnames=[
            'Chr',
            'start',
            'end',
            'name',
            'twistName',
            'hg19Name'],
        delimiter='\t')

    liftoverList = ['Intron', 'Intergenic', 'TES', 'TSS']

    print('Processing input bed file.')
    for line in reader:
        targets = [x.strip() for x in line['description'].split(',')]

        primary = targets[0]
        chrom = line['Chr']
        start = line['start']
        end = line['end']
        try:
            secondary = targets[1:]
        except IndexError:
            secondary = []

        fields = primary.split(':')
        gene = fields[0]
        featureType = fields[len(fields) - 1]

        if featureType in liftoverList:
            try:
                noOrtho[primary].append([fixChr(chrom), start, end])
            except BaseException:
                noOrtho[primary] = [[fixChr(chrom), start, end]]

        if gene not in reference.keys():
            try:
                gene = hg19Lookup[gene]
            except KeyError:
                pass
        try:
            seqType = fields[3]
        except IndexError:
            pass

        if gene in reference.keys() and featureType not in [
                'Intron', 'Intergenic']:

            for entry in reference[gene]:

                if entry.orthologType == 'ortholog_one2one':

                    try:
                        primeo2o[entry.dogGeneID].add(seqType)
                    except KeyError:
                        primeo2o[entry.dogGeneID] = set(seqType)

                    pone2one.writerow({'Chr': chrom,
                                       'start': start,
                                       'end': end,
                                       'name': primary,
                                       'confidence(o/1)': entry.orthologConfidence,
                                       'geneName': entry.humanGeneName,
                                       'geneID': entry.geneID,
                                       'dogGeneName': entry.dogGeneName,
                                       'dogGeneID': entry.dogGeneID,
                                       'seqType': seqType})

                elif entry.orthologType == 'ortholog_one2many':

                    try:
                        primeo2m[entry.dogGeneID].add(seqType)
                    except KeyError:
                        primeo2m[entry.dogGeneID] = set(seqType)

                    pone2many.writerow({'Chr': chrom,
                                        'start': start,
                                        'end': end,
                                        'name': primary,
                                        'confidence(o/1)': entry.orthologConfidence,
                                        'geneName': entry.humanGeneName,
                                        'geneID': entry.geneID,
                                        'dogGeneName': entry.dogGeneName,
                                        'dogGeneID': entry.dogGeneID,
                                        'seqType': seqType})

                else:
                    if isNumber(chrom):
                        chrw = ('chr' + str(chrom))
                    elif chrom == 'MT':
                        chrw = ('chrM')
                    elif chrom == 'X':
                        chrw = ('chrX')
                    else:
                        chrw = ('chr' + str(chrom))
                    pnoOrtho.writerow({'Chr': chrw,
                                       'start': start,
                                       'end': end,
                                       'name': primary,
                                       'geneName': primary.split(':')[0].strip()})

        else:
            unk.writerow({'Chr': fixChr(chrom),
                          'start': start,
                          'end': end,
                          'name': primary,
                          'twistName': primary.split(':')[0].strip(),
                          'hg19Name': gene})
            try:
                noOrtho[primary].append([fixChr(chrom), start, end])
            except KeyError:
                noOrtho[primary] = [[fixChr(chrom), start, end]]

        for target in secondary:
            fields = target.split(':')
            featureType = fields[len(fields) - 1].strip()
            gene = fields[0].strip()
            if featureType in liftoverList:
                try:
                    noOrtho[target].append([fixChr(chrom), start, end])
                except BaseException:
                    noOrtho[target] = [[fixChr(chrom), start, end]]

            if gene not in reference.keys():
                try:
                    gene = hg19Lookup[gene]

                except KeyError:
                    pass

            if gene in reference.keys():
                for entry in reference[gene]:
                    if entry.orthologType == 'ortholog_one2one':
                        try:
                            secondaryo2o[entry.dogGeneID].add(seqType)
                        except KeyError:
                            secondaryo2o[entry.dogGeneID] = set(seqType)

                        sone2one.writerow({'Chr': chrom,
                                           'start': start,
                                           'end': end,
                                           'name': target,
                                           'confidence(o/1)': entry.orthologConfidence,
                                           'geneName': entry.humanGeneName,
                                           'geneID': entry.geneID,
                                           'dogGeneName': entry.dogGeneName,
                                           'dogGeneID': entry.dogGeneID,
                                           'seqType': seqType})

                    elif entry.orthologType == 'ortholog_one2many':
                        for orthotype in entry.orthologType:
                            try:
                                secondaryo2m[entry.dogGeneID].add(seqType)
                            except KeyError:
                                secondaryo2m[entry.dogGeneID] = set(seqType)

                            sone2many.writerow({'Chr': chrom,
                                                'start': start,
                                                'end': end,
                                                'name': target,
                                                'confidence(o/1)': entry.orthologConfidence,
                                                'geneName': entry.humanGeneName,
                                                'geneID': entry.geneID,
                                                'dogGeneName': entry.dogGeneName,
                                                'dogGeneID': entry.dogGeneID,
                                                'seqType': seqType})

                    else:
                        snoOrtho.writerow({'Chr': fixChr(chrom),
                                           'start': start,
                                           'end': end,
                                           'name': target,
                                           'geneName': target.split(':')[0].strip()})

            else:
                unk.writerow({'Chr': fixChr(chrom),
                              'start': start,
                              'end': end,
                              'name': target,
                              'twistName': target.split(':')[0].strip(),
                              'hg19Name': gene})
                
                try:
                    noOrtho[target].append([fixChr(chrom), start, end])
                except KeyError:
                    noOrtho[target] = [[fixChr(chrom), start, end]]


with open('oneToOnePrimary.tmpbed', 'w') as otopBed, \
        open('oneToManyPrimary.tmpbed', 'w') as otmpBed, \
        open('oneToOneSecondary.tmpbed', 'w') as otosBed, \
        open('oneToManySecondary.tmpbed', 'w') as otmsBed, \
        open('notin95.out', 'w') as weird, \
        open('canFam3CDS.tmpbed', 'w') as CDSbed, \
        open('canFam3rRNA.tmpbed', 'w') as rRNAbed, \
        open('canFam3lncRNA.tmpbed', 'w') as lncbed, \
        open('canFam3snRNA.tmpbed', 'w') as snbed, \
        open('canFam3snoRNA.tmpbed', 'w') as snobed, \
        open('canFam3miRNA.tmpbed', 'w') as mibed:
    db = gffutils.FeatureDB(args.dogDb)
    #CDS = db.features_of_type('CDS')
    EXONS = db.features_of_type('exon')
    for i in EXONS:
        if 'rRNA' in i['gene_biotype']:
            chrom = fixChr(i.chrom)
            start = i.start -1
            end = i.end
            name = i['gene_id'][0]
            biotype = i['gene_biotype'][0]
            strand = i.strand
            try:
                geneName = db[name]['gene_name'][0]
            except KeyError:
                geneName = ''
            
            rRNAbed.write(
                chrom +
                '\t' +
                str(start) +
                '\t' +
                str(end) +
                '\t' +
                geneName + 
                ':' +
                name +
                ':' +
                biotype +
                ':Exon:' +
                strand +
                '\n')
        
        elif 'lncRNA' in i['gene_biotype']:
            chrom = fixChr(i.chrom)
            start = i.start -1
            end = i.end
            name = i['gene_id'][0]
            biotype = i['gene_biotype'][0]
            strand = i.strand
            try:
                geneName = db[name]['gene_name'][0]
            except KeyError:
                geneName = ''
            
            lncbed.write(
                chrom +
                '\t' +
                str(start) +
                '\t' +
                str(end) +
                '\t' +
                geneName + 
                ':' +
                name +
                ':' +
                biotype +
                ':Exon:' +
                strand +
                '\n')
            
        elif 'snRNA' in i['gene_biotype']:
            chrom = fixChr(i.chrom)
            start = i.start -1
            end = i.end
            name = i['gene_id'][0]
            biotype = i['gene_biotype'][0]
            strand = i.strand
            try:
                geneName = db[name]['gene_name'][0]
            except KeyError:
                geneName = ''
            
            snbed.write(
                chrom +
                '\t' +
                str(start) +
                '\t' +
                str(end) +
                '\t' +
                geneName + 
                ':' +
                name +
                ':' +
                biotype +
                ':Exon:' +
                strand +
                '\n')
            
        elif 'snoRNA' in i['gene_biotype']:
            chrom = fixChr(i.chrom)
            start = i.start -1
            end = i.end
            name = i['gene_id'][0]
            biotype = i['gene_biotype'][0]
            strand = i.strand
            try:
                geneName = db[name]['gene_name'][0]
            except KeyError:
                geneName = ''
            
            snobed.write(
                chrom +
                '\t' +
                str(start) +
                '\t' +
                str(end) +
                '\t' +
                geneName + 
                ':' +
                name +
                ':' +
                biotype +
                ':Exon:' +
                strand +
                '\n')
            
        elif 'miRNA' in i['gene_biotype']:
            chrom = fixChr(i.chrom)
            start = i.start -1
            end = i.end
            name = i['gene_id'][0]
            biotype = i['gene_biotype'][0]
            strand = i.strand
            try:
                geneName = db[name]['gene_name'][0]
            except KeyError:
                geneName = ''
            
            mibed.write(
                chrom +
                '\t' +
                str(start) +
                '\t' +
                str(end) +
                '\t' +
                geneName + 
                ':' +
                name +
                ':' +
                biotype +
                ':Exon:' +
                strand +
                '\n')
    
    db = gffutils.FeatureDB(args.dogDb)
    CDS = db.features_of_type('CDS')        
    
    for i in CDS:
        chrom = fixChr(i.chrom)
        start = i.start -1
        end = i.end
        name = i['gene_id'][0]
        biotype = i['gene_biotype'][0]
        strand = i.strand
        try:
            geneName = db[name]['gene_name'][0]
        except KeyError:
            geneName = ''
            
        CDSbed.write(
            chrom +
            '\t' +
            str(start) +
            '\t' +
            str(end) +
            '\t' +
            geneName + 
            ':' +
            name +
            ':' +
            biotype +
            ':CDS:' +
            strand +
            '\n')

    for entry in primeo2o:
        # print(entry)
        try:
            gene = db[entry]
            try:
                name = gene.attributes['gene_name'][0]
            except KeyError:
                name = ''

            if 'Exon' in primeo2o[entry]:
                exons = [x for x in db.children(gene, featuretype='exon')]
                exonIntervals = []
                # BEGINFIXINGHERE
                for chrom, start, end, strand, biotype in zip([x.chrom for x in exons],
                                                   [int(x.start -1) for x in exons],
                                                   [int(x.end) for x in exons],
                                                   [x.strand for x in exons],
                                                   [x['gene_biotype'][0] for x in exons]):
                    exonIntervals.append([chrom, start, end, strand, biotype])
                for interval in exonIntervals:
                    otopBed.write(fixChr(interval[0]) + 
                                  '\t' + 
                                  str(interval[1]) +
                                  '\t' + 
                                  str(interval[2]) + 
                                  '\t' + 
                                  name +
                                  ':' + 
                                  entry + 
                                  ':' + 
                                  str(interval[4]) + 
                                  ':Exon:' + 
                                  interval[3] +
                                  '\n')

            elif 'CDS' in primeo2o[entry]:
                cdss = [x for x in db.children(gene, featuretype='CDS')]
                cdsIntervals = []
                for chrom, start, end, strand, biotype in zip([x.chrom for x in cdss], 
                                                   [int(x.start -1) for x in cdss], 
                                                   [int(x.end) for x in cdss], 
                                                   [x.strand for x in cdss],
                                                   [x['gene_biotype'][0] for x in cdss]):
                    cdsIntervals.append([chrom, start, end, strand, biotype])
                for interval in cdsIntervals:
                    otopBed.write(fixChr(interval[0]) +
                                  '\t' +
                                  str(interval[1]) +
                                  '\t' +
                                  str(interval[2]) +
                                  '\t' +
                                  name +
                                  ':' +
                                  entry +
                                  ':' +
                                  str(interval[4]) +
                                  ':CDS:' +
                                  interval[3] +
                                  '\n')

            elif 'TSS' in primeo2o[entry]:
                exons = [x for x in db.children(gene, featuretype='exon')]
                exonIntervals = []
                for chrom, start, end, strand, biotype in zip([x.chrom for x in exons],
                                                   [int(x.start -1) for x in exons],
                                                   [int(x.end) for x in exons],
                                                   [x.strand for x in exons],
                                                   [x['gene_biotype'][0] for x in exons]):
                    exonIntervals.append([chrom, start, end, strand, biotype])
                for interval in exonIntervals:
                    otopBed.write(fixChr(interval[0]) +
                                  '\t' +
                                  str(interval[1]) +
                                  '\t' +
                                  str(interval[2]) +
                                  '\t' +
                                  name +
                                  ':' +
                                  entry +
                                  ':' +
                                  str(interval[4]) +
                                  ':Exon-TSSinTWIST:' +
                                  interval[3] +
                                  '\n')

            elif 'TES' in primeo2o[entry]:
                exons = [x for x in db.children(gene, featuretype='exon')]
                exonIntervals = []

                for chrom, start, end, strand, biotype  in zip([x.chrom for x in exons],
                                                   [int(x.start -1) for x in exons],
                                                   [int(x.end) for x in exons],
                                                   [x.strand for x in exons],
                                                   [x['gene_biotype'][0] for x in exons]):
                    exonIntervals.append([chrom, start, end, strand, biotype])
                for interval in exonIntervals:
                    otopBed.write(fixChr(interval[0]) +
                                  '\t' +
                                  str(interval[1]) +
                                  '\t' +
                                  str(interval[2]) +
                                  '\t' +
                                  name +
                                  ':' +
                                  entry +
                                  ':' + 
                                  str(interval[4]) +
                                  ':Exon-TESinTWIST:' +
                                  interval[3] +
                                  '\n')

        except KeyError:
            weird.write(entry + '\n')

    for entry in primeo2m:
        # print(entry)
        try:
            gene = db[entry]
            try:
                name = gene.attributes['gene_name'][0]
            except KeyError:
                name = ''

            if 'Exon' in primeo2m[entry]:
                exons = [x for x in db.children(gene, featuretype='exon')]
                exonIntervals = []
                for chrom, start, end, strand, biotype in zip([x.chrom for x in exons],
                                                   [int(x.start -1) for x in exons],
                                                   [int(x.end) for x in exons],
                                                   [x.strand for x in exons],
                                                   [x['gene_biotype'][0] for x in exons]):
                    exonIntervals.append([chrom, start, end, strand, biotype])
                for interval in exonIntervals:
                    otmpBed.write(fixChr(interval[0]) +
                                  '\t' +
                                  str(interval[1]) +
                                  '\t' +
                                  str(interval[2]) +
                                  '\t' +
                                  name +
                                  ':' +
                                  entry +
                                  ':' + 
                                  str(interval[4]) +
                                  ':Exon:' +
                                  interval[3] +
                                  '\n')

            elif 'CDS' in primeo2m[entry]:
                cdss = [x for x in db.children(gene, featuretype='CDS')]
                cdsIntervals = []
                for chrom, start, end, strand, biotype in zip([x.chrom for x in cdss], 
                                                   [int(x.start -1) for x in cdss], 
                                                   [int(x.end) for x in cdss], 
                                                   [x.strand for x in cdss],
                                                   [x['gene_biotype'][0] for x in cdss]):
                    cdsIntervals.append([chrom, start, end, strand, biotype])
                for interval in cdsIntervals:
                    otmpBed.write(fixChr(interval[0]) +
                                  '\t' +
                                  str(interval[1]) +
                                  '\t' +
                                  str(interval[2]) +
                                  '\t' +
                                  name +
                                  ':' +
                                  entry +
                                  ':' +
                                  str(interval[4]) +
                                  ':CDS:' +
                                  interval[3] +
                                  '\n')

            elif 'TSS' in primeo2m[entry]:
                exons = [x for x in db.children(gene, featuretype='exon')]
                exonIntervals = []
                for chrom, start, end, strand, biotype in zip([x.chrom for x in exons], 
                                                   [int(x.start -1) for x in exons], 
                                                   [int(x.end) for x in exons], 
                                                   [x.strand for x in exons],
                                                   [x['gene_biotype'][0] for x in exons]):
                    exonIntervals.append([chrom, start, end, strand, biotype])
                for interval in exonIntervals:
                    otmpBed.write(fixChr(interval[0]) +
                                  '\t' +
                                  str(interval[1]) +
                                  '\t' +
                                  str(interval[2]) +
                                  '\t' +
                                  name +
                                  ':' +
                                  entry +
                                  ':' +
                                  str(interval[4]) +
                                  ':Exon-TSSinTWIST:' +
                                  interval[3] +
                                  '\n')

            elif 'TES' in primeo2m[entry]:
                exons = [x for x in db.children(gene, featuretype='exon')]
                exonIntervals = []
                for chrom, start, end, strand, biotype in zip([x.chrom for x in exons], 
                                                   [int(x.start -1) for x in exons], 
                                                   [int(x.end) for x in exons], 
                                                   [x.strand for x in exons],
                                                   [x['gene_biotype'][0] for x in exons]):
                    exonIntervals.append([chrom, start, end, strand, biotype])
                for interval in exonIntervals:
                    otmpBed.write(fixChr(interval[0]) +
                                  '\t' +
                                  str(interval[1]) +
                                  '\t' +
                                  str(interval[2]) +
                                  '\t' +
                                  name +
                                  ':' +
                                  entry +
                                  ':' +
                                  str(interval[4]) +
                                  ':Exon-TESinTWIST:' +
                                  interval[3] +
                                  '\n')

        except KeyError:
            weird.write(entry + '\n')

    for entry in secondaryo2o:
        # print(entry)
        try:
            gene = db[entry]
            try:
                name = gene.attributes['gene_name'][0]
            except KeyError:
                name = ''

            if 'Exon' in secondaryo2o[entry]:
                exons = [x for x in db.children(gene, featuretype='exon')]
                exonIntervals = []
                # BEGINFIXINGHERE
                for chrom, start, end, strand, biotype in zip([x.chrom for x in exons],
                                                   [int(x.start -1) for x in exons],
                                                   [int(x.end) for x in exons],
                                                   [x.strand for x in exons],
                                                   [x['gene_biotype'][0] for x in exons]):
                    exonIntervals.append([chrom, start, end, strand, biotype])
                for interval in exonIntervals:
                    otosBed.write(fixChr(interval[0]) + 
                                  '\t' + 
                                  str(interval[1]) +
                                  '\t' + 
                                  str(interval[2]) + 
                                  '\t' + 
                                  name +
                                  ':' + 
                                  entry + 
                                  ':' + 
                                  str(interval[4]) + 
                                  ':Exon:' + 
                                  interval[3] +
                                  '\n')

            elif 'CDS' in secondaryo2o[entry]:
                cdss = [x for x in db.children(gene, featuretype='CDS')]
                cdsIntervals = []
                for chrom, start, end, strand, biotype in zip([x.chrom for x in cdss], 
                                                   [int(x.start -1) for x in cdss], 
                                                   [int(x.end) for x in cdss], 
                                                   [x.strand for x in cdss],
                                                   [x['gene_biotype'][0] for x in cdss]):
                    cdsIntervals.append([chrom, start, end, strand, biotype])
                for interval in cdsIntervals:
                    otosBed.write(fixChr(interval[0]) +
                                  '\t' +
                                  str(interval[1]) +
                                  '\t' +
                                  str(interval[2]) +
                                  '\t' +
                                  name +
                                  ':' +
                                  entry +
                                  ':' +
                                  str(interval[4]) +
                                  ':CDS:' +
                                  interval[3] +
                                  '\n')

            elif 'TSS' in secondaryo2o[entry]:
                exons = [x for x in db.children(gene, featuretype='exon')]
                exonIntervals = []
                for chrom, start, end, strand, biotype in zip([x.chrom for x in exons],
                                                   [int(x.start -1) for x in exons],
                                                   [int(x.end) for x in exons],
                                                   [x.strand for x in exons],
                                                   [x['gene_biotype'][0] for x in exons]):
                    exonIntervals.append([chrom, start, end, strand, biotype])
                for interval in exonIntervals:
                    otosBed.write(fixChr(interval[0]) +
                                  '\t' +
                                  str(interval[1]) +
                                  '\t' +
                                  str(interval[2]) +
                                  '\t' +
                                  name +
                                  ':' +
                                  entry +
                                  ':' +
                                  str(interval[4]) +
                                  ':Exon-TSSinTWIST:' +
                                  interval[3] +
                                  '\n')

            elif 'TES' in secondaryo2o[entry]:
                exons = [x for x in db.children(gene, featuretype='exon')]
                exonIntervals = []

                for chrom, start, end, strand, biotype  in zip([x.chrom for x in exons],
                                                   [int(x.start -1) for x in exons],
                                                   [int(x.end) for x in exons],
                                                   [x.strand for x in exons],
                                                   [x['gene_biotype'][0] for x in exons]):
                    exonIntervals.append([chrom, start, end, strand, biotype])
                for interval in exonIntervals:
                    otosBed.write(fixChr(interval[0]) +
                                  '\t' +
                                  str(interval[1]) +
                                  '\t' +
                                  str(interval[2]) +
                                  '\t' +
                                  name +
                                  ':' +
                                  entry +
                                  ':' + 
                                  str(interval[4]) +
                                  ':Exon-TESinTWIST:' +
                                  interval[3] +
                                  '\n')

        except KeyError:
            weird.write(entry + '\n')

    
            
    for entry in secondaryo2m:
        # print(entry)
        try:
            gene = db[entry]
            try:
                name = gene.attributes['gene_name'][0]
            except KeyError:
                name = ''

            if 'Exon' in secondaryo2m[entry]:
                exons = [x for x in db.children(gene, featuretype='exon')]
                exonIntervals = []
                for chrom, start, end, strand, biotype in zip([x.chrom for x in exons],
                                                   [int(x.start -1) for x in exons],
                                                   [int(x.end) for x in exons],
                                                   [x.strand for x in exons],
                                                   [x['gene_biotype'][0] for x in exons]):
                    exonIntervals.append([chrom, start, end, strand, biotype])
                for interval in exonIntervals:
                    otmsBed.write(fixChr(interval[0]) +
                                  '\t' +
                                  str(interval[1]) +
                                  '\t' +
                                  str(interval[2]) +
                                  '\t' +
                                  name +
                                  ':' +
                                  entry +
                                  ':' + 
                                  str(interval[4]) +
                                  ':Exon:' +
                                  interval[3] +
                                  '\n')

            elif 'CDS' in secondaryo2m[entry]:
                cdss = [x for x in db.children(gene, featuretype='CDS')]
                cdsIntervals = []
                for chrom, start, end, strand, biotype in zip([x.chrom for x in cdss], 
                                                   [int(x.start -1) for x in cdss], 
                                                   [int(x.end) for x in cdss], 
                                                   [x.strand for x in cdss],
                                                   [x['gene_biotype'][0] for x in cdss]):
                    cdsIntervals.append([chrom, start, end, strand, biotype])
                for interval in cdsIntervals:
                    otmsBed.write(fixChr(interval[0]) +
                                  '\t' +
                                  str(interval[1]) +
                                  '\t' +
                                  str(interval[2]) +
                                  '\t' +
                                  name +
                                  ':' +
                                  entry +
                                  ':' +
                                  str(interval[4]) +
                                  ':CDS:' +
                                  interval[3] +
                                  '\n')

            elif 'TSS' in secondaryo2m[entry]:
                exons = [x for x in db.children(gene, featuretype='exon')]
                exonIntervals = []
                for chrom, start, end, strand, biotype in zip([x.chrom for x in exons], 
                                                   [int(x.start -1) for x in exons], 
                                                   [int(x.end) for x in exons], 
                                                   [x.strand for x in exons],
                                                   [x['gene_biotype'][0] for x in exons]):
                    exonIntervals.append([chrom, start, end, strand, biotype])
                for interval in exonIntervals:
                    otmsBed.write(fixChr(interval[0]) +
                                  '\t' +
                                  str(interval[1]) +
                                  '\t' +
                                  str(interval[2]) +
                                  '\t' +
                                  name +
                                  ':' +
                                  entry +
                                  ':' +
                                  str(interval[4]) +
                                  ':Exon-TSSinTWIST:' +
                                  interval[3] +
                                  '\n')

            elif 'TES' in secondaryo2m[entry]:
                exons = [x for x in db.children(gene, featuretype='exon')]
                exonIntervals = []
                for chrom, start, end, strand, biotype in zip([x.chrom for x in exons], 
                                                   [int(x.start -1) for x in exons], 
                                                   [int(x.end) for x in exons], 
                                                   [x.strand for x in exons],
                                                   [x['gene_biotype'][0] for x in exons]):
                    exonIntervals.append([chrom, start, end, strand, biotype])
                for interval in exonIntervals:
                    otmsBed.write(fixChr(interval[0]) +
                                  '\t' +
                                  str(interval[1]) +
                                  '\t' +
                                  str(interval[2]) +
                                  '\t' +
                                  name +
                                  ':' +
                                  entry +
                                  ':' +
                                  str(interval[4]) +
                                  ':Exon-TESinTWIST:' +
                                  interval[3] +
                                  '\n')

        except KeyError:
            weird.write(entry + '\n')
    

        # elif 'TSS'
with open('LiftoverFeatureTypes.bed', 'w')as toLiftOver:
    for entry in noOrtho.keys():
        for interval in noOrtho[entry]:
            toLiftOver.write(
                interval[0] +
                '\t' +
                interval[1] +
                '\t' +
                interval[2] +
                '\t' +
                entry +
                '\n')


print('Merging and sorting output bed files.')
#os.system('sort -k1,1 -k2,2n canFam3CDS.tmpbed | bedtools merge -i stdin > canFam3CDS.bed; rm canFam3CDS.tmpbed')
os.system('cat oneToOnePrimary.tmpbed oneToManyPrimary.tmpbed oneToOneSecondary.tmpbed oneToManySecondary.tmpbed > allOrtholog.tmpbed')
#os.system('cat > noOrtholog.tmpbed')
os.system('awk \'{print $1"\t"$2"\t"$3"\t"$4}\' *noOrtho.tmpout > noOrtholog.tmpbed')
os.system('awk \'{print $1"\t"$2"\t"$3"\t"$4}\' unknownGenes.tmpout >> noOrtholog.tmpbed')
os.system('for file in *.tmpbed; do sort -k1,1 -k2,2n $file | uniq > $(basename $file .tmpbed).bed; rm $file; done;') # rm $file; done;')
#os.system('reuse BEDTools')
os.system('bedtools merge -i allOrtholog.bed -c 4 -o distinct > mergedOrthologs.bed')
#os.system('rm allOrtholog.bed')
os.system('cat LiftoverFeatureTypes.bed noOrtholog.bed cancerAddition.bed | sort -k1,1 -k2,2n | uniq | awk \'{print $1"\t"$2"\t"$3"\t"$4":Liftover"}\' > toLiftover.bed')
#os.system('bedtools merge -i noOrtholog.bed -c 4 -o distinct > mergednoOrthologs.bed')
#os.system('rm noOrtholog.bed')
#os.system('rm LiftoverFeatureTypes.bed')
os.system('awk \'{print $0"="$3-$2}\' toLiftover.bed > toLiftoverNumbered.bed')

print('Merging in all CDS regions from /seq/vgb/swofford/ref/Canis_familiaris.CanFam3.1.99.gtf')
os.system('bedtools subtract -a canFam3CDS.bed -b mergedOrthologs.bed >canFamCDSnotInOrthologs.bed')
os.system('cat mergedOrthologs.bed canFamCDSnotInOrthologs.bed | sort -k1,1 -k2,2n | bedtools merge -i stdin -c 4 -o distinct> mergedOrthologsPlusCDS.bed')
#os.system('bedtools merge -i canFamOrthologPlusCDS.bed -c 4 -o distinct >mergedOrthologsPlusCDS.bed')
#os.system('rm canFamOrthologPlusCDS.bed')
print('Lifting over intervals missing orthologs')
os.system('liftOver -minMatch=0.9 toLiftoverNumbered.bed /seq/vgb/swofford/software/liftover/hg19ToCanFam3.over.chain noOrthologLiftedOver.tmpbed noOrthologFailedLiftover.bed')
print('Checking size of lifted over intervals, keeping all that fall between 90-110% of original size')
with open('noOrthologLiftedOver.tmpbed', 'r') as inBed, open('noOrthologLiftedOverSizeChecked.tmpbed', 'w') as outBed, open('wrongSizeLiftover.bed', 'w') as wrongSize:
    for line in inBed:
        origSize = int(line.split('=')[1])
        size = int(line.split('\t')[2]) - int(line.split('\t')[1])
        if 1.1 >= size / origSize >= 0.9:
            outBed.write(line.split('=')[0] + '\n')
        else:
            wrongSize.write(line)


os.system('sort -k1,1 -k2,2n noOrthologLiftedOverSizeChecked.tmpbed > noOrthologLiftedOver.bed')
#os.system('rm noOrthologLiftedOverSizeChecked.tmpbed')
os.system('bedtools merge -i noOrthologLiftedOver.bed -c 4 -o distinct >mergedNoOrthologLiftedOver.bed')
#os.system('rm noOrthologLiftedOver.bed')
print('Merging in passing lifted over intervals. Results file: mergedOrthologCDSLiftover.bed')
os.system('cat mergedNoOrthologLiftedOver.bed mergedOrthologsPlusCDS.bed | sort -k1,1 -k2,2n | bedtools merge -i stdin -c 4 -o distinct > mergedOrthologCDSLiftover.bed')

'''
print('Checking overlap of liftover and stringtie R24 assemblies')
print('Creating stringtie.bed from /seq/vgb-HAL/swofford/RNAseq/assembliesStringtieFlagAndHisatFlag/stringtie/stringMerged/stringtie_Merged.gtf')
os.system(
    'grep exon /seq/vgb-HAL/swofford/RNAseq/assembliesStringtieFlagAndHisatFlag/stringtie/stringMerged/stringtie_Merged.gtf | awk \'{print $1\"\t\"$4\"\t\"$5\"\t\"$9$10$11$12$13$14$15$16$17$18$19$20}\' | sort -k1,1 -k2,2n | bedtools merge -i stdin -c 4 -o distinct > stringtie.bed')
print('Removing overlap between stringtie and dog CDS/ortholog')
os.system('bedtools subtract -a stringtie.bed -b mergedOrthologsPlusCDS.bed >stringtieUnique.bed')
print('Intersecting liftover with stringtie. Output: strintieUniqueLiftoverIntersect.bed')
os.system('bedtools intersect -a stringtieUnique.bed -b mergedNoOrthologLiftedOver.bed | sort -k1,1 -k2,2n | bedtools merge -i stdin -c 4 -o distinct > stringtieUniqueLiftoverIntersect.bed')


print('Pulling full exon list for stringtie hits. output: stringtieTranscriptsFromLiftover.bed')
with open('stringtieUniqueLiftoverIntersect.bed', 'r') as inBed, open('stringtieTranscriptsFromLiftover.bed', 'w') as outBed:
    db = gffutils.FeatureDB('R24db')
    transcripts = set()
    for line in inBed:
        gene = line.split('gene_id')[1].split('\"')[1]
        transcripts.add(gene)
    for gene in transcripts:
        exons = [x for x in db.children(gene, featuretype='exon')]
        exonIntervals = []

        for chrom, start, end, strand in zip([x.chrom for x in exons], [int(
                x.start -1) for x in exons], [int(x.end) for x in exons], [x.strand for x in exons]):
            exonIntervals.append([chrom, start, end, strand])
        for interval in exonIntervals:
            outBed.write(interval[0] +
                         '\t' +
                         str(interval[1]) +
                         '\t' +
                         str(interval[2]) +
                         '\t' +
                         gene +
                         ':ExonFromStringtie:' +
                         interval[3] +
                         '\n')
print('Merging stringtie exons into final. output: mergedOrthologCDSLiftoverStringtie.bed')
os.system('cat mergedNoOrthologLiftedOver.bed mergedOrthologsPlusCDS.bed stringtieTranscriptsFromLiftover.bed | sort -k1,1 -k2,2n | bedtools merge -i stdin -c 4 -o distinct > mergedOrthologCDSLiftoverStringtie.bed')


print('Checking overlap of liftover and Broad improved annotation protein coding gtf')
print('Creating broadAnno.bed from /seq/vgb/swofford/ref/canFam3ImprovedAnnoProteinCoding.gtf')
os.system(
    'grep exon /seq/vgb/swofford/ref/canFam3ImprovedAnnoProteinCoding.gtf | awk \'{print $1"\t"$4"\t"$5"\t"$9$10$11$12$13$14$15$16$17$18$19$20}\' | sort -k1,1 -k2,2n | bedtools merge -i stdin -c 4 -o distinct > broadAnno.bed')
print('Removing overlap between Broad and dog CDS/ortholog')
os.system(
    'bedtools subtract -a broadAnno.bed -b mergedOrthologsPlusCDS.bed >broadUnique.bed')
print('Intersecting liftover with Broad. Output: broadUniqueLiftoverIntersect.bed')
os.system('bedtools intersect -a broadUnique.bed -b mergedNoOrthologLiftedOver.bed | sort -k1,1 -k2,2n | bedtools merge -i stdin -c 4 -o distinct > broadUniqueLiftoverIntersect.bed')


print('Pulling full cds list for broad improved anno hits. output: broadTranscriptsFromLiftover.bed')
#os.system('grep exon /seq/vgb/swofford/ref/canFam3ImprovedAnnoProteinCoding.gtf | awk \'{print $1"\t"$4"\t"$5"\t"$9$10$11$12$13$14$15$16$17$18$19$20}\' | sort -k1,1 -k2,2n | bedtools merge -i stdin -c 4 -o distinct > canFam3ImprovedProteinExon.bed')


with open('broadUniqueLiftoverIntersect.bed', 'r') as inBed, open('broadTranscriptsFromLiftover.bed', 'w') as outBed:
    db = gffutils.FeatureDB('broadImprovedProteinOnlydb')
    transcripts = set()
    for line in inBed:
        gene = line.split('gene_id')[1].split('\"')[1]
        transcripts.add(gene)
    for gene in transcripts:
        exons = [x for x in db.children(gene, featuretype='CDS')]
        exonIntervals = []
        #chr = db[gene].chrom
        #strand = db[gene].strand
        for chrom, start, end, strand in zip([x.chrom for x in exons], [int(
                x.start -1) for x in exons], [int(x.end) for x in exons], [x.strand for x in exons]):
            exonIntervals.append([chrom, start, end, strand])
        for interval in exonIntervals:
            outBed.write(interval[0] +
                         '\t' +
                         str(interval[1]) +
                         '\t' +
                         str(interval[2]) +
                         '\t' +
                         gene +
                         ':CDSFromBroad:' +
                         interval[3] +
                         '\n')
print('Merging Broad anno CDSs into final. output: mergedOrthologCDSLiftoverBroad.bed')
os.system('cat mergedNoOrthologLiftedOver.bed mergedOrthologsPlusCDS.bed broadTranscriptsFromLiftover.bed | sort -k1,1 -k2,2n | bedtools merge -i stdin -c 4 -o distinct > mergedOrthologCDSLiftoverBroad.bed')


#os.system('bedtools intersect -a /seq/vgb-HAL')
# Using bash to sort and remove duplicate entries for each file
os.system(
    'for file in *.tmpout; do awk \'NR == 1; NR > 1 {print $0 | "sort -n"}\' $file | uniq > $(basename $file .tmpout).out; rm $file; done;')


os.system('grep \':CDS\' primary_one2one.out >primary_one2oneCDS.out')
os.system('grep \':Exon\' primary_one2one.out >primary_one2oneExon.out')
os.system('grep \':Intron\' primary_one2one.out >primary_one2oneIntron.out')
os.system('grep \':TES\' primary_one2one.out >primary_one2oneTES.out')
os.system('grep \':TSS\' primary_one2one.out >primary_one2oneTSS.out')
os.system('grep \':CDS\' secondary_one2one.out >secondary_one2oneCDS.out')
os.system('grep \':Exon\' secondary_one2one.out >secondary_one2oneExon.out')
os.system('grep \':Intron\' secondary_one2one.out >secondary_one2oneIntron.out')
os.system('grep \':TES\' secondary_one2one.out >secondary_one2oneTES.out')
os.system('grep \':TSS\' secondary_one2one.out >secondary_one2oneTSS.out')
os.system('grep \'Intergenic\' unknownGenes.out > intergenic.out')
os.system('grep -v \'Intergenic\' unknownGenes.out > unknownGenesNoIntergenic.out')
'''

# Liftover style: SPATA31A3:protein_coding:+:CDS:Liftover