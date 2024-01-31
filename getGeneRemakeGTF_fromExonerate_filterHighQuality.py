import gffutils
#from Bio import AlignIO
#from Bio.AlignIO import MafIO
import os
#import pysam
#from pysam import FastaFile
import argparse
#from textwrap import fill

import pickle
from Bio.Seq import Seq
from fuzzywuzzy import fuzz
import glob
import pysam
from pysam import FastaFile
from builtins import FileNotFoundError










parser = argparse.ArgumentParser()

parser.add_argument('-i', '--geneID', default = '')
parser.add_argument('-f', '--geneFile', default = '')
parser.add_argument('-s', '--speciesFile', default = '/seq/vgb/swofford/zoonomia/speciesList.txt')

args = parser.parse_args()
species = args.speciesFile
species = species.split('/')[1].split('.')[0]

def complement(seq):
    '''Returns the complement of a DNA sequence passed to it. Accepts a
    string as input, and returns a string. Valid bases for input are
    "ATCGW"'''
    compSeq = ''
    codes = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'W': 'W', 'N': 'N'}
    for letter in seq:
        compSeq += codes[letter]
    return(compSeq)


def reverse(seq):
    '''Returns the reverse of a string passed to it. Input is a string,
    output is a string.'''
    return(seq[::-1])


def runGene(geneID):
    #updated version will only run with geneID as input here
    
    print(geneID)       
    positions = []

    gene = db[geneID]
        
    try:
        name = db[gene]['gene_name'][0]
    except KeyError:
        name = geneID
    
    if gene.strand == '-':
        strand = '-1'
    elif gene.strand == '+':
        strand = '1'
    '''check to see if human transcripts have already been written by a previous instance
    all species can share the same human transcript info, so no need to redo this every time
    '''
    if not os.path.exists(name + '/transcripts.pickle'):
        try:
            os.mkdir(name)
            os.mkdir(name + '/humanBeds')
            os.mkdir(name + '/liftoverBeds')
            os.mkdir(name + '/WholeGene')
            os.mkdir(name + '/CDS')
        except:
            pass
        
    
        transcripts = {}
        #get info from each transcript from the human annotation
        for i in db.children(geneID, featuretype='transcript'):
            if i['transcript_biotype'][0] == 'protein_coding':
                cdsStarts = []
                cdsEnds = []
                cdsExonNumber = []
                for j in db.children(i, featuretype='CDS', order_by='start'):
                    cdsStarts.append(j.start -1)
                    cdsEnds.append(j.end)
                    cdsExonNumber.append(j['exon_number'][0])
                    protein_id = j['protein_id'][0]
                exons = []
                for start, end, number in zip(cdsStarts, cdsEnds, cdsExonNumber):
                    exons.append([start,end,number])    
                
                transcripts[i['transcript_id'][0]] = {'exons': exons}
                transcripts[i['transcript_id'][0]]['name'] = i['transcript_name'][0]
                transcripts[i['transcript_id'][0]]['strand'] = i.strand
                transcripts[i['transcript_id'][0]]['chrom'] = 'chr' + str(i.chrom)
                transcripts[i['transcript_id'][0]]['protein'] = protein_id
                
        #I don't remember what's going on with the below and why it's gone, revisit if necessary        
        '''
        proteinList = set()
        for transcript in transcripts.keys():
            proteinList.add(transcripts[transcript]['protein'])
        '''
            
        with open(name + '/humanBeds/' + name + '.bed', 'w') as outBed:
            outBed.write('chr' + str(gene.chrom) + '\t' + str(gene.start - 1) + '\t' 
                      + str(gene.end) + '\n')
        
        pepTideFileNames = {}
        '''
        for protein in proteinList:
            print(protein)
        '''
        transcriptLifts = {}
        seenIntervals = []
        
        '''get the dna and protein sequence for every human transcript''' 
        for transcript in transcripts.keys():
            #print('transcript keys ' + transcript)
            for cds in db.children(transcript, featuretype = 'CDS'):
                #print('transcript has CDS?')
                #if cds['protein_id'][0] == protein:
                transcriptName = transcripts[transcript]['name'] + '_' + transcript
                strand = transcripts[transcript]['strand']
                chrom = transcripts[transcript]['chrom']
                exons = transcripts[transcript]['exons']
                
                intervals = []
                for exon in exons:
                    intervals.append([exon[0],exon[1]])
                    
                #if intervals not in seenIntervals:
                #    print('unseen intervals ' + transcript)
                seenIntervals.append(intervals)
                fileName = name + '/humanBeds/' + transcriptName + '.bed'
                with open(fileName, 'w') as outBed:
                    for exon in exons:
                        outBed.write(chrom + '\t' + str(exon[0]) + '\t'
                                     + str(exon[1]) + '\t'
                                     + transcriptName + '_exon_' + str(exon[2]) + '\t.\t'
                                     + strand + '\n')
                if strand == '-':
                    strandOption = '-1'
                else:
                    strandOption = '1'
                    
                os.system('bedtools getfasta -fi /seq/vgb/swofford/zoonomia/tempFasta/Homo_sapiens.fa -bed ' + fileName + ' > ./' + name + '/CDS/Homo_sapiens.' + transcriptName + '.CDS.fa')
                os.system('union -sequence ./' + name + '/CDS/Homo_sapiens.' + transcriptName + '.CDS.fa -outseq ./' + name + '/CDS/Homo_sapiens.' + transcriptName + '.CDS.union.fa')
                os.system('transeq ./' + name + '/CDS/Homo_sapiens.' + transcriptName + '.CDS.union.fa ' +
                          './' + name + '/CDS/Homo_sapiens.' + transcriptName + '.CDS.union.pep -frame=' + strandOption)
                
                transcriptLifts[transcript] = name + '/CDS/Homo_sapiens.' + transcriptName + '.CDS.union.pep'
        #save dictionaries of transcripts and sequence file locations for use in other species            
        with open(name + '/transcriptLifts.pickle', 'wb') as pickFile:
            pickle.dump(transcriptLifts, pickFile)
        with open(name + '/transcripts.pickle', 'wb') as pickFile:
            pickle.dump(transcripts, pickFile)
        
    

    try:
        with open(name + '/transcriptLifts.pickle', 'rb') as pickFile:
            transcriptLifts = pickle.load(pickFile)
        with open(name + '/transcripts.pickle', 'rb') as pickFile:
            transcripts = pickle.load(pickFile)
    except:
        #this should only happen when the gene id entered does not have any protein coding transcripts listed in the ensembl annotation
        print('No transcripts found, skipping gene')
        return
    
    for i in species:
        for transcript in transcripts.keys():         
            if not os.path.exists(name + '/WholeGene/' + i + '_' + transcript + '.exonerate'):
                if not os.path.exists(name + '/liftoverBeds/' + i + name + '.WholeGene.bed'):
                    os.system('while read -r line; do halLiftover /seq/vgb/200m_project/alignment_files/241-mammalian-2020v2.hal Homo_sapiens ' 
                              + name + '/humanBeds/' + name + '.bed $line ' + name + '/liftoverBeds/"$line"' + name 
                              + '.WholeGene.bed; done < ' + args.speciesFile)
            
        
        #length of pad to add to either end of gene
                pad = 500
            
            #print('smoothing bed results')
            #for i in species:
                with open(name + '/liftoverBeds/' + i + name + '.WholeGene.bed', 'r') as bed, open(name + '/liftoverBeds/' + i + name + '.WholeGenePad.bed', 'w') as outBed:
                    contigs = {}
                    for line in bed:
                        line = line.strip()
                        con, start, end = line.split('\t')
                        try:
                            contigs[con].append([int(start), int(end)])
                        except KeyError:
                            contigs[con] = [[start,end]]
                    for j in contigs.keys():
                        starts = []
                        ends = []
                        for pos in contigs[j]:
                            starts.append(int(pos[0]))
                            ends.append(int(pos[1]))
            
                        start = min(starts)
                        if start - pad > 0:
                            start = start - pad
                        else: start = 0
                        end = max(ends)
                        if end + pad < int(contigLookup[i][j]):
                            end = end + pad
                        else: end = int(contigLookup[i][j])
                        
                        outBed.write(j + '\t' + str(start) + '\t' + str(end) + '\n')
            
                os.system('bedtools getfasta -fi /seq/vgb/swofford/zoonomia/tempFasta/' + i + '.fa -bed ' + name + 
                          '/liftoverBeds/' + i + name + '.WholeGenePad.bed > ./' + name + '/WholeGene/' + i + '.WholeGene.fa')
                
                fileName = './' + name + '/WholeGene/' + i + '.WholeGene.fa' 
        
                for transcript in transcriptLifts.keys():  
                #print('here! ' + transcript) 
                    if not os.path.exists(name + '/WholeGene/' + i + '_' + transcript + '.exonerate'):
                        try:
                            #print('here2!')
                            #os.system('exonerate -m protein2dna -q ' + transcriptLifts[transcript] + ' -t ' + fileName + 
                            os.system('exonerate -m protein2genome --showtargetgff -q ' + transcriptLifts[transcript] + ' -t ' + fileName + 
                                ' --ryo "%tcs" -n 1 > ./' + name + '/WholeGene/' + i + '_' + transcript + '.exonerate')
                        except FileNotFoundError:
                            #print('here3!')
                            continue
            else:
                print(name + '/WholeGene/' + i + '_' + transcript + '.exonerate already exists, skipping liftover')       
    
        '''
        For the GTF entry, we want to make one entry for the gene: 
        chr, source, gene, start, end, . ,strand, . ,
        'gene_id "ENSG..._species_zoonomia"; gene_version(1 for broad, 2 for pfenning?);
        gene_name "name"; gene_source "zoonomia"; gene_biotype "protein_coding";'
        
        one for each transcript:
        chr, source, transcript, start, end, ., strand, ., gene_id "ENSG..._species_zoonomia";
        gene_version "(same as above?)"; transcript_id "ENST...zoonomia"; transcript_version "(same as gene?)";
        gene_name "name"; gene_source "zoonomia"; gene_biotype "protein_coding"; transcript_name; "name-201"; 
        transcript_source "zoonomia"; transcript_biotype "protein_coding";
        
        one for each exon:
        chr, source, exon, start, end, ., strand, ., gene_id; gene_version; transcript_id; transcript_version;
        exon_number "1"; gene_name name; gene_source; gene_biotype; transcript_name; transcript_source;
        transcript_biotype; exon_id "ENSE..." (I think we want to make up our own here, not always going to be 1:1);
        exon_version "1";
        
        Ideally, could do something with the 'basic' tags and/or the support level to see if both methods agree on the same
        transcripts and exons, but that feels like second version stuff
        '''
        
        
        gene_id = geneID
        gene_version = '1'
        gene_name = name
        gene_source = 'zoonomia_broad'
        gene_biotype = 'protein_coding'
        transcript_name = ''
        transcript_source = 'zoonomia_broad'
        transcript_biotype = 'protein_coding'
        tag = 'basic'
        transcript_support_leve = '1'
        exon_number = ''
        exon_version = ''
        contig = ''
        strand = ''
        strandList = []
        
        print('Translating exonerate back to bed')
        bedFiles = []
        transcriptGTF = {}
        highQualTranscriptGTF = {}
        allExons = set()
        #for i in species:
        geneMin = 0
        geneMax = 0
        identity = ''
        similarity = ''
        
        ######Picking only longest transcript to use -remove section if desired to go back to all transcripts
        
        transcriptLengths = {}
        longest = ''
        lengths = []
        for transcript in transcripts.keys():
            exons = transcripts[transcript]['exons']
            intervals = []
            for exon in exons:
                intervals.append(exon[0])
                intervals.append(exon[1])
            length = max(intervals) - min(intervals)
            transcriptLengths[transcript] = length
            lengths.append(length)
        longest = ''
        lengths.sort(reverse = True)
        
        for transcript in transcriptLengths:
            if longest == '' or transcriptLengths[transcript] > transcriptLengths[longest]:
                longest = transcript
                
            
        
        longestHumanProt = ''
        try:
            humanFile = glob.glob(name + '/CDS/Homo_sapiens.*_' + longest + '.CDS.union.pep')[0]
        except IndexError:
            #print('here')
            continue
        with open(humanFile, 'r') as hum:
            humLines = []
            for hline in hum:
                if hline.strip()[0] != '>':
                    humLines.append(hline.strip())
                longestHumanProt = ''.join(humLines)
        
        if not os.path.exists(name + '/WholeGene/' + i + '_' + transcript + '.exonerate'):
            print('longest transcript failed')
            for j in os.listdir(name + '/WholeGene/'):
                if i in j:
                    transcript = j.split('.')[0].split('_')[-1]
                    #print(transcript)      
        
        #try:
        #    transcript = pfenLookup[gene_name][0]
        #except KeyError:
        #    print(gene_name + ' missing')
        #    transcript = longest
        foundTranscript = False
        #i = i for i in transcripts.keys()
        y = 0
        transcript = ''
        toWrite = []
        toWriteHighQual = []
        while not foundTranscript:
            #print(str(y))
            #print(transcriptLengths.keys())
            try:
                for x in transcriptLengths.keys():
                    if transcriptLengths[x] == lengths[y]:
                        transcript = x
                        #print(x)
            except IndexError:
                foundTranscript = True
                continue                
            if transcript not in allTranscripts:
                #print('super here ' + transcript)
                try:
                    humanFile = glob.glob(name + '/CDS/Homo_sapiens.*_' + transcript + '.CDS.union.pep')[0]
                    #print(humanFile)
                except IndexError:
                    #print('here')
                    #print(name + '/CDS/Homo_sapiens.*_' + transcript + '.CDS.union.pep')
                    y += 1
                    continue
                with open(humanFile, 'r') as hum:
                    humLines = []
                    for hline in hum:
                        if hline.strip()[0] != '>':
                            humLines.append(hline.strip())
                        humanProt = ''.join(humLines)
                #for transcript in transcriptLifts.keys():
                    #transcriptGTF[transcript] = {}
                    #strand = transcripts[transcript]['strand']
                transcript_name = transcripts[transcript]['name']
                gene_id = geneID
                gene_max = 0
                gene_min = 0
                pfenMin = 0
                pfenMax = 0
                try:
                    
                    humanFile = glob.glob(name + '/CDS/Homo_sapiens.*_' + transcript + '.CDS.union.pep')[0]
                    with open(name + '/WholeGene/' + i + '_' + transcript + '.exonerate', 'r') as inFile, open(humanFile, 'r') as hum:
                        humLines = []
                        for hline in hum:
                            if hline.strip()[0] != '>':
                                humLines.append(hline.strip())
                        humanProt = ''.join(humLines)
                    
                    with open(name + '/WholeGene/' + i + '_' + transcript + '.exonerate', 'r') as inFile:
                        gtfLines = []
                        gtfStart = False
                        seqStart = False
                        seqLines = []
                        highQual = False
                        keepSet = ['A','C','G','T','a','c','g','t']
                        for eline in inFile:
                            #print(eline)
                            #print(gtfStart)
                            if 'START OF GFF DUMP' in eline:
                                gtfStart = True
                                continue
                            elif gtfStart:
                                if 'END OF GFF DUMP' in eline:
                                    gtfStart = False
                                    seqStart = True
                                    continue
                            
                            if gtfStart and '#' not in eline:
                                gtfLines.append(eline.split('\t'))
                            try:
                                if seqStart and eline.strip()[0] in keepSet:
                                    seqLines.append(eline.strip())
                            except IndexError:
                                break
                                
                        #sequence = Seq(''.join(seqLines)).translate()
                        
                        #if len(sequence) > 0 and sequence[0] == 'M' and '*' not in sequence and len(sequence) >= .9 * len(humanProt) and len(sequence) <= 1.1 * len(humanProt):# and fuzz.ratio(sequence, humanProt) >= 70:
                        #    print('broad result: ' + sequence)
                        #    highQualTranscriptGTF[transcript_name] = gtfLines
                        #elif len(sequence) < .9 * len(humanProt):
                            
                        
                    transcriptGTF[transcript_name] = gtfLines
                                    
                    
                    #toWrite = []
                    #toWriteHighQual = []  
                    intervals = []  
                    tstrand = ''
                    #for tline in transcriptGTF.keys():
                        #######################
                    #    entry = transcriptGTF[tline]
                        
                    for tfields in gtfLines:
                        #for tfields in entry:
                            #contig = tfields[0].split(':')[0]
                        bedStart,bedEnd = tfields[0].split(':')[1].split('-')
                        #print("bed stuff: " + bedStart + ':' + bedEnd)
                        if tfields[2] == 'cds':
                            contig = tfields[0].split(':')[0]
                            start = int(tfields[0].split(':')[1].split('-')[0])
                            end = int(tfields[0].split(':')[1].split('-')[1])
                            intervals.append([contig, int(tfields[3]) + int(bedStart), int(tfields[4]) + int(bedStart)])
                            tstrand = tfields[6]
                    #print('broad intervals = ' + str(intervals))            
                                
                    try:
                        phases = []
                        phase = 0
                        intervals.sort()
                        sequence = ''
                        #gene_min = int(intervals[0][1])
                        #gene_max = int(intervals[-1][2])
                        with pysam.Fastafile('/seq/vgb/swofford/zoonomia/tempFasta/' + i + '.fa') as fasta:
                            for interval in intervals:
                                sequence = ''.join ([sequence, fasta.fetch(interval[0], int(interval[1])-1, interval[2]).upper().strip()])
                        if tstrand == '-':
                            sequence = reverse(complement(sequence))
                        #print('broad sequence: ' + sequence)
                        sequence = Seq(sequence).translate()
                        #print('broad sequence: ' + sequence)
                        if len(sequence) > 0 and sequence[0] == 'M' and '*' not in sequence and len(sequence) >= .9 * len(humanProt) and len(sequence) <= 1.1 * len(humanProt):
                            #print('am I here?')
                            highQualTranscriptGTF[transcript_name] = gtfLines
                            #print(highQualTranscriptGTF)
                        else:
                            print('broad seqBeg = ' + sequence[0])
                            print('broad Too Short = ' + str(len(sequence) >= .9 * len(humanProt)))
                            print('broad Too Long = ' + str(len(sequence) <= 1.1 * len(humanProt)))
                    #except IOError:
                    #except (KeyError, UnboundLocalError, IndexError):
                    except (UnboundLocalError, IndexError):
                        #print('here2')
                        pass
                    #####################        
                            
                            
                        
                         
                    for tline in transcriptGTF.keys():
                        if tline in highQualTranscriptGTF.keys():
                            highQual = True
                        else:
                            highQual = False
                            contig = ''
                            continue
                        entry = transcriptGTF[tline]
                        
                        #print(tfields)
                        #print(tfields[0])
                        for tfields in entry:
                            '''
                            bedStart,bedEnd = tfields[0].split(':')[1].split('-')
                            if geneMin == 0 or geneMin > (int(tfields[3]) + int(bedStart)):
                                geneMin = int(tfields[3]) + int(bedStart)
                            if geneMax < (int(tfields[4]) + int(bedStart)):
                                geneMax = int(tfields[4]) + int(bedStart)
                            '''
                            contig = tfields[0].split(':')[0]
                            bedStart,bedEnd = tfields[0].split(':')[1].split('-')
                            if tfields[2] == 'gene':
                                type = 'transcript'
                                identity = (tfields[8].split(';')[3].strip() + ';' + tfields[8].split(';')[4].strip())
                                identity.replace(' ;',';')
                                ifields = identity.split(';')
                                for j in range(len(ifields)):
                                    entry, score = ifields[j].strip().split()
                                    ifields[j] = entry + ' "' + score + '"' 
                                identity = '; '.join(ifields)
                                start = str(int(tfields[3]) + int(bedStart))
                                end = str(int(tfields[4]) + int(bedStart))
                                positions.append(int(start))
                                positions.append(int(end))
                                #print(positions)
                                tstrand = tfields[6]
                                toWrite.append('\t'.join([contig, 'zoonomia_broad', type, start, end, '.', tstrand, '.',
                                                      'gene_id "' + gene_id + '_' + i + '_zoonomia"; gene_version "'
                                                            + gene_version + '"; gene_name "' + name
                                                            + '"; transcript_id "' + transcript + '_' + i + '_zoonomia"; transcript_version "1' 
                                                            + '"; gene_source "zoonomia"; gene_biotype "protein_coding"; ' + identity + ';\n']))
                                if highQual:
                                    toWriteHighQual.append('\t'.join([contig, 'zoonomia_broad', type, start, end, '.', tstrand, '.',
                                                      'gene_id "' + gene_id + '_' + i + '_zoonomia"; gene_version "'
                                                            + gene_version + '"; gene_name "' + name
                                                            + '"; transcript_id "' + transcript + '_' + i + '_zoonomia"; transcript_version "1' 
                                                            + '"; gene_source "zoonomia"; gene_biotype "protein_coding"; ' + identity + ';\n']))
                                    
                            elif tfields[2] == 'exon':
                                identity = tfields[8].strip()
                                ifields = identity.split(';')
                                for j in range(len(ifields)):
                                    entry, score = ifields[j].strip().split()
                                    ifields[j] = entry + ' "' + score + '"' 
                            
                                identity = '; '.join(ifields)
                                start = str(int(tfields[3]) + int(bedStart))
                                end = str(int(tfields[4]) + int(bedStart))
                                positions.append(int(start))
                                positions.append(int(end))
                                tstrand = tfields[6]
                                strandList.append(tstrand)
                                type = 'exon'
                                #similarity = tfields[8].split(';')[4].strip()
                                toWrite.append('\t'.join([contig, 'zoonomia_broad', type, start, end, '.', tstrand, '.',
                                                      'gene_id "' + gene_id + '_' + i + '_zoonomia"; gene_version "'
                                                            + gene_version + '"; gene_name "' + name
                                                            + '"; transcript_id "' + transcript + '_' + i + '_zoonomia"; transcript_version "1' 
                                                            + '"; gene_source "zoonomia"; gene_biotype "protein_coding"; ' + identity + ';\n']))
                                if highQual:
                                    toWriteHighQual.append('\t'.join([contig, 'zoonomia_broad', type, start, end, '.', tstrand, '.',
                                                      'gene_id "' + gene_id + '_' + i + '_zoonomia"; gene_version "'
                                                            + gene_version + '"; gene_name "' + name
                                                            + '"; transcript_id "' + transcript + '_' + i + '_zoonomia"; transcript_version "1' 
                                                            + '"; gene_source "zoonomia"; gene_biotype "protein_coding"; ' + identity + ';\n']))
                                    
                            
                            elif tfields[2] == 'cds':
                                type = 'CDS'
                                start = str(int(tfields[3]) + int(bedStart))
                                end = str(int(tfields[4]) + int(bedStart))
                                positions.append(int(start))
                                positions.append(int(end))
                                tstrand = tfields[6]
                                strandList.append(tstrand) 
                                toWrite.append('\t'.join([contig, 'zoonomia_broad', type, start, end, '.', tstrand, '.',
                                                          'gene_id "' + gene_id + '_' + i + '_zoonomia"; gene_version "'
                                                                + gene_version + '"; gene_name "' + name
                                                                + '"; transcript_id "' + transcript + '_' + i + '_zoonomia"; transcript_version "1' 
                                                                + '"; gene_source "zoonomia"; gene_biotype "protein_coding";\n']))
                                if highQual:
                                    toWriteHighQual.append('\t'.join([contig, 'zoonomia_broad', type, start, end, '.', tstrand, '.',
                                                          'gene_id "' + gene_id + '_' + i + '_zoonomia"; gene_version "'
                                                                + gene_version + '"; gene_name "' + name
                                                                + '"; transcript_id "' + transcript + '_' + i + '_zoonomia"; transcript_version "1' 
                                                                + '"; gene_source "zoonomia"; gene_biotype "protein_coding";\n']))
                                    
                            
                    posStrand = 0
                    negStrand = 0
                    for j in strandList:
                        if j == '+':
                            posStrand +=1
                        if j == '-':
                            negStrand +=1
                    
                    if posStrand >= negStrand:
                        gstrand = '+'
                    else:
                        gstrand = '-'
                    #print(positions)
                    try:
                        gene_max = max(positions)
                        #print(gene_max)
                        gene_min = min(positions)
                        positions = []
                    except ValueError:
                        gene_max = 0
                        gene_min = 0
                #except IOError:
                except (FileNotFoundError, IndexError):
                    print('no exonerate file found')
                
                if len(toWriteHighQual) > 0:
                    foundTranscript = True
                else:
                    y += 1
            else:
                if len(toWriteHighQual) > 0:
                    foundTranscript = True
                else:
                    y += 1
                        
 
                
        pfenTrans = ''    
        pfenWrite = []
        pfenWriteHQ = []
        try:
            pfenMin = 0
            pfenMax = 0
            #print(pfenLookup[gene_name][1])
            #print('/seq/vgb/swofford/zoonomia/wholeGenomeAnnotation/pfenningGTF/' + species[0] + '/' + pfenLookup[gene_name][1] + '.gtf')
            with open('/seq/vgb/swofford/zoonomia/wholeGenomeAnnotation/pfenningGTF/' + species[0] + '/' + pfenLookup[gene_name][1] + '.gtf', 'r') as pfenFile, pysam.Fastafile('/seq/vgb/swofford/zoonomia/tempFasta/' + i + '.fa') as fasta: 
                intervals = []
                sequence = ''
                strand = ''
                strands = set()
                contigSet = set()
                for pfenline in pfenFile:
                    #intervals = []
                    #sequence = ''
                    if pfenline.split('\t')[2] == 'transcript':
                        pfenTrans = pfenline
                        strand = pfenline.split('\t')[6]
                        
                        #pfenWrite.insert(pfenTrans,0)
                        
                    elif pfenline.split('\t')[2] == 'CDS':
                        if int(pfenline.split('\t')[3]) <= int(pfenline.split('\t')[4]):
                            fields = pfenline.split('\t')
                            strands.add(pfenline.split('\t')[6])
                            contigSet.add(pfenline.split('\t')[0])
                            intervals.append([fields[0], int(fields[3]), int(fields[4])])
                            #meta = fields[8].split(';')
                            #meta[0] = meta[0][:-1] + '_exon_' + meta[4].split('"')[1] + '"'
                            #meta[2] = meta[2][:-1] + '_exon_' + meta[4].split('"')[1] + '"'
                            #meta[0] = meta[0][:-1] + '_' + 'zoonomia"'
                            #meta[2] = meta[2][:-1] + '_' + 'zoonomia"'
                            #fields[8] = ';'.join(meta)
                            #pfenWrite.append(pfenline)
                            pfenWrite.append('\t'.join(fields))
                
                try:
                    phases = []
                    phase = 0
                    intervals.sort()
                    pfenMin = int(intervals[0][1])
                    pfenMax = int(intervals[-1][2])
                    for interval in intervals:
                        #print(fasta.fetch(interval[0], int(interval[1])-1, interval[2]).upper().strip())
                        sequence = ''.join ([sequence, fasta.fetch(interval[0], int(interval[1])-1, interval[2]).upper().strip()])
                    #print(pfenline.split('\t')[6].strip())
                    if strand == '-':
                        #print('minus strand')
                        sequence = reverse(complement(sequence))
                    
                    if strand == '+':
                        for interval in intervals:
                            phases.append((3 - phase) %3 )
                            phase =  ((int(interval[2]) - (int(interval[1]) -1)) + phase) % 3
                        
                    else:
                        intervals.reverse()
                        for interval in intervals:
                            phases.append((3 - phase) %3 )
                            phase =  ((int(interval[2]) - (int(interval[1]) -1)) + phase) % 3
                            
                            #phase =  (3 - (phase + (int(interval[2]) - (int(interval[1]) -1))) % 3) % 3
                        phases.reverse()
                    
                    
                except (KeyError, UnboundLocalError, IndexError):
                    pass
                #print(intervals)
                #print(sequence)
                
                sequence = Seq(sequence).translate()
                #print(sequence)
                try:
                    if contig == '':
                        contig = list(contigSet)[0]
                        gstrand = list(strands)[0]
                    elif list(contigSet)[0] != contig or list(strands)[0] != gstrand:
                        contigSet.add(0)
                except IndexError:
                    contigSet = [0,0]    
                    
                if len(strands) == 1 and len(contigSet) == 1 and len(sequence) > 0 and sequence[0] == 'M' and '*' not in sequence and len(sequence) >= .9 * len(longestHumanProt) and len(sequence) <= 1.1 * len(longestHumanProt): #and fuzz.ratio(sequence, humanProt) >= 70: 
                    #print('fuzzRatio = ' + str(fuzz.ratio(sequence, humanProt)))
                    #meta = fields[8].split(';')
                    #meta[0] = '_'.join(meta[0].split('_')[:-1]) + '"'
                    #fields[8] = ';'.join(meta)
                    #pfenWriteHQ = pfenWrite
                    #print(phases)
                    
                            
                    
                    for entry, phase in zip(pfenWrite, phases):
                        #print('pfenLine = ' + entry)
                        #pfenWriteHQ.append(entry.replace('CDS','exon'))
                        fields = entry.split('\t')
                        #meta = fields[8].split(';')
                        #meta[0] = '_'.join(meta[0].split('_')[:-1]) + '"'
                        #meta[2] = '_'.join(meta[2].split('_')[:-1]) + '"'
                        #fields[8] = ';'.join(meta)
                        #print('fieldsAfter = ' + str(fields))
                        pfenWriteHQ.append('\t'.join(fields).replace('CDS','exon'))
                        fields[7] = str(phase)
                        
                        #pfenWriteHQ.append('\t'.join(fields).replace('CDS','exon'))
                        pfenWriteHQ.append('\t'.join(fields))
                        
                else:
                    try:
                        print('pfen seqBeg = ' + sequence[0])
                    except IndexError:
                        print('pfen seqBeg = NA')
                    print('pfen Too Short = ' + str(len(sequence) >= .9 * len(longestHumanProt)))
                    print('pfen Too Long = ' + str(len(sequence) <= 1.1 * len(longestHumanProt)))
                    print('pfen num strands = ' + str(len(strands)))
                    
                    pfenMin = 0
                    pfenMax = 0
                        
                            
                            
                            
                            
                            
        except (FileNotFoundError, KeyError):
            print(transcript + ' not in pfenning annotation')
                        
                    
            ''''       
            with open(line + '.gtf', 'a') as gtfSummary:
                
                if contig != '' and gene_min != 0 and gene_max != 0:
                    gtfSummary.write(contig + '\t' + gene_source + '\tgene\t' + str(gene_min) 
                                                     + '\t' + str(gene_max) + '\t.\t' + gstrand + '\t.\t'
                                                     + 'gene_id "' + gene_id + '_' + line 
                                                     + '_zoonomia"; gene_version "' + gene_version 
                                                     + '"; gene_name "' + name 
                                                     + '"; gene_source "zoonomia"; gene_biotype "protein_coding";\n')
                    for i in range(len(toWrite)):
                        if toWrite[i].split('\t')[2] == 'transcript':
                            gtfSummary.write(toWrite[i])
                        elif toWrite[i].split('\t')[2] == 'exon':
                            gtfSummary.write(toWrite[i])
                            gtfSummary.write(toWrite[i-1])
            
                #gtfSummary.write(pfenTrans)
                for pline in pfenWrite:
                    gtfSummary.write(pline)
            '''
        with open(i + '_highQual.gtf', 'a') as gtfSummaryHigh:                              
            intervals = []
            strand = ''
            strands = set()
            contigSet = set()
            for broadline in toWriteHighQual:         
                if broadline.split('\t')[2] == 'CDS':
                    strand = broadline.split('\t')[6]
                    strands.add(strand)
                    contigSet.add(broadline.split('\t')[0])
                    if int(broadline.split('\t')[3]) <= int(broadline.split('\t')[4]):
                        fields = broadline.split('\t')
                        intervals.append([fields[0], int(fields[3]), int(fields[4])])
            
            try:
                phases = {}
                phase = 0
                intervals.sort()
                
                if strand == '+':
                    for interval in intervals:
                        phases[str(interval)] = ((3 - phase) %3)
                        phase =  ((int(interval[2]) - (int(interval[1]) -1)) + phase) % 3
                    
                else:
                    intervals.reverse()
                    for interval in intervals:
                        phases[str(interval)] = ((3 - phase) %3)
                        phase =  ((int(interval[2]) - (int(interval[1]) -1)) + phase) % 3
                        
                        #phase =  (3 - (phase + (int(interval[2]) - (int(interval[1]) -1))) % 3) % 3
                    #phases.reverse()
            
                
                
            except KeyError:
                continue
            toWriteHighQualFin = []
            if len(strands) == 1 and len(contigSet) == 1:
                for broadline in toWriteHighQual: 
                    if broadline.split('\t')[2] == 'CDS':        
                        #print(entry)
                        fields = broadline.split('\t')
                        interval = [fields[0], int(fields[3]), int(fields[4])]
                        #print(phases)
                        fields[7] = str(phases[str(interval)])
                        newEntry = '\t'.join(fields)
                        toWriteHighQualFin.append(newEntry)
                    else:
                        toWriteHighQualFin.append(broadline)
                    
           
            #print('length of toWriteHighQual = ' + str(len(toWriteHighQual)))
            broadMin = 0
            broadMax = 0
            
            if len(toWriteHighQualFin) > 0:
                broadMin = gene_min
                broadMax = gene_max
                #print(str(gene_min) + ': ' + str(gene_max))
                
                if pfenMin < broadMin and pfenMin > 0:
                    gene_min = pfenMin
                if pfenMax > broadMax and pfenMax > 0:
                    gene_max = pfenMax
            elif len(pfenWriteHQ) > 0:
                gene_min = pfenMin
                gene_max = pfenMax
            else:
                gene_min = 0
                gene_max = 0
            
            if contig != '' and (len(toWriteHighQualFin) > 0 or len(pfenWriteHQ) > 0) and gene_min != 0 and gene_max != 0:
                #print('yaaaaaaay gene: ' + line + ' pfenWrite: ' + str(pfenWrite) + ' broadWrite: ' + str(toWriteHighQualFin))
                gtfSummaryHigh.write(contig + '\t' + 'zoonomia' + '\tgene\t' + str(gene_min) 
                                                 + '\t' + str(gene_max) + '\t.\t' + gstrand + '\t.\t'
                                                 + 'gene_id "' + gene_id + '_' + i 
                                                 + '_zoonomia"; gene_version "' + gene_version 
                                                 + '"; gene_name "' + name 
                                                 + '"; gene_source "zoonomia"; gene_biotype "protein_coding";\n')
                for j in range(len(toWriteHighQualFin)):
                    if toWriteHighQualFin[j].split('\t')[2] == 'transcript':
                        gtfSummaryHigh.write(toWriteHighQualFin[j])
                    elif toWriteHighQual[j].split('\t')[2] == 'exon':
                        gtfSummaryHigh.write(toWriteHighQualFin[j])
                    else:    
                        #fields = toWriteHighQual[i-1].split('\t')
                        #interval = [fields[0], int(fields[3]), int(fields[4])]
                        #fields[7] = str(phases[str(interval)])
                        #entry = '\t'.join(fields)
                        gtfSummaryHigh.write(toWriteHighQualFin[j])
                        
                        #toWriteHighQual[i-1]
                        #gtfSummaryHigh.write(toWriteHighQual[i-1])
                        
            if len(pfenWriteHQ) >0:
                gtfSummaryHigh.write(pfenTrans)
                for pline in pfenWriteHQ:
                    gtfSummaryHigh.write(pline)            
            
                
        allTranscripts.append(transcript)            
                            
                
    for critter in species:
        try:
            with open(critter + '/genesRun.out', 'a') as genesRun:
                genesRun.write(geneID + '\n')
        except FileNotFoundError:
            os.mkdir(critter)
            with open(critter + '/genesRun.out', 'a') as genesRun:
                genesRun.write(geneID + '\n')            




allTranscripts = []                               
pfenLookup = {}
with open('/seq/vgb/swofford/zoonomia/wholeGenomeAnnotation/forPfenningWholeGenome.out', 'r') as pfenList:
    for line in pfenList:
        trans = line.split('\t')[0]
        prot = line.split('\t')[1]
        genName = line.split('\t')[3].strip()
        #print(genName)
        pfenLookup[genName] = [trans,prot]
                        
species = []
alreadyRun = []
print('getting species list')
with open(args.speciesFile, 'r') as speciesList:
    for line in speciesList:
        species.append(line.strip())

for critter in species:
    if os.path.exists(critter + '/genesRun.out'):
        with open(critter + '/genesRun.out', 'r') as genesRun:
            for gline in genesRun:
                alreadyRun.append(gline.strip())

#if not os.path.exists('contigs.pickle'):         

contigLookup = {}    
print('getting list of contig lengths')
for i in species:    
    with open('/seq/vgb/swofford/zoonomia/tempFasta/' + i + '.fa.fai', 'r') as index:
        contigs = {}
        for line in index:
            contig, contigLength = line.split()[0], line.split()[1]
            contigs[contig] = contigLength
        contigLookup[i] = contigs
    #print('done building lookup, pickling')
#with open('contigs.pickle', 'wb') as pickFile:
    #pickle.dump(contigLookup, pickFile)
#else:
    #with open('contigs.pickle', 'rb') as pickFile:
       # contigLookup = pickle.load(pickFile)

print('getting gene boundaries')

if os.path.exists('/seq/vgb/swofford/zoonomia/hg38.99.db'):
    db = gffutils.FeatureDB('/seq/vgb/swofford/zoonomia/hg38.99.db')
else:
    db = gffutils.create_db('/seq/vgb/swofford/ref/Homo_sapiens.GRCh38.99.chr_patch_hapl_scaff.gtf', 
                            dbfn='hg38.99', 
                            force = True, keep_order = True, 
                            merge_strategy='merge', 
                            sort_attribute_values=True, 
                            disable_infer_genes=True, 
                            disable_infer_transcripts=True)


if args.geneFile != '':
    with open(args.geneFile, 'r') as geneFile:
        for line in geneFile:
            if line.strip() not in alreadyRun:
            #if not os.path.exists(db[line.strip()]['gene_name'][0]):
                runGene(line.strip())
            else:
                print('gene already run: ' + line.strip())
else:
    if args.geneID not in alreadyRun:
    #if not os.path.exists(db[args.geneID]['gene_name'][0]):
        runGene(args.geneID)
    else:
        print('gene already run: ' + args.geneID)
                
