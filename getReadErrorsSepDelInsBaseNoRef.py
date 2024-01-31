import pysam
import argparse
import re
import sys




parser = argparse.ArgumentParser()
parser.add_argument('-i', '--inFile')
parser.add_argument('-o', '--outFile')
parser.add_argument('-v', '--varFile')
parser.add_argument('-r', '--rightOnly', action='store_true')
parser.add_argument('-l', '--leftOnly', action='store_true')
args = parser.parse_args()

def isNumber(num):
    try:
        int(num)
        return True
    except(TypeError):
        return False
    
def buildRef():
    '''Returns a dictionary of ordered oligo variant inserts and names'''
    with open(args.varFile, 'r') as ref:
        variants = {}
        for line in ref:
            fields = re.split('\s+', line.rstrip())
            name = str(fields[0]) + '_' + str(fields[1]) + \
                '_' + str(fields[2]) + '_' + str(fields[9])
            #for full oligo use -> var = str(fields[12], 10 for shortened)
            var = str(fields[10])
            name = name.replace(',', '^')
            variants[name] = var
        ref.close()
        return variants

if args.varFile:
    refVars = buildRef()
else:
    refVars = {}   
# get list of all ref sequence names, these will corespond to regions in the                                                                                                                                                                                                                                 
# reference FASTA used to create the bam files.                                                                                                                                                                                                                                                                  
if args.leftOnly:
    refVars['leftRef'] = 'AACTGGCCGCTTGACGT'
elif args.rightOnly:
    refVars['rightRef'] = 'CACTGCGGCTCCTGCGATCGCGTCGACGAACCTCTAGA'

regions = []    
for name, var in refVars.items():
        regions.append(name)

with pysam.AlignmentFile(args.inFile, 'rb') as bam, open(args.outFile, 'w') as outFile, open(args.outFile + 'dels', 'w') as delFile, open(args.outFile + 'ins', 'w') as insFile, open(args.outFile + 'base', 'w') as baseFile:
#    outFile.write('numDeletions\tdeletedBases\tlongestDeletion\tnumInsertions\tinsertedBases\tlongestInsertion\tnumSubstitutions\n')

    #instantiate array of 0s to track errors per base
    errorSummary = []
    delSummary = []
    insSummary = []
    baseErrorSummary = []
    count = 0
    readCount = 0
    
    for i in range(600):
            baseErrorSummary.append(0)
            errorSummary.append(0)
            delSummary.append(0)
            insSummary.append(0)

    #iterRead = bam.fetch()
    
    #if len(regionRead) < 2:
        #regionRead[0]     
    for region in regions:
            #outFile.write(str(region) + '\n')
        print(str(region))
        actualRef = str(refVars[region])
        refLength = len(actualRef)
        print(str(len(actualRef)))
        iterRead = bam.fetch(str(region))    
            
        for read in iterRead:
            regionRead = []
            readCount += 1
            #md = read.get_tag('MD')
            #cigar = read.cigartuples
            alignment = read.get_aligned_pairs(with_seq=True)
            baseCount = 0
            adjustLeft = 0
            adjustRight = 0
            regionRead = []
            #firstIns = True
            #print(read.get_reference_sequence())
            #print(read.query_sequence)
                
    
    
            #check each position in the alignment for errors
            #count = read.reference_start
            count = 0
            adjustLeft = 0
            lastRef = 0
            firstRef = 'nan'
            firstSeq = 'nan'
            for pair in alignment:
                if isNumber(pair[1]):
                    if firstRef == 'nan':
                        firstRef = pair[1]
                    lastRef = pair[1]
                    #if pair[1] != count:
                    #    print(str(pair[1]) + ' count: ' + str(count))
                        #deletion
                    #    regionRead.append(3)
                    if isNumber(pair[0]):
                        if firstSeq == 'nan':
                            firstSeq = pair[0]
                        if pair[0] + adjustLeft < pair[1]:
                            #deletion
                            for i in range((pair[1] - (pair[0] + adjustLeft))):
                                regionRead.append(3)
                                adjustLeft += 1
                        if str(pair[2]).isupper():
                            regionRead.append(0)
                        else:
                            regionRead.append(1)
                    #else:
                        #deletion
                    #    regionRead.append(3)
                #elif isNumber(pair[0]):
                #    regionRead.append(1)
                count += 1
            for i in range(firstRef):
                regionRead.insert(0,3)        
            #if refLength > lastRef:
                #regionRead.append(3)
             #   for n in range(refLength - lastRef - 1):
             #       regionRead.append(3)
                #print(cigar)
                #print(read.get_tag('MD'))
                #print(str(alignment))
                #print(read.query_sequence)
            '''
            if '^' in (read.get_tag('MD')):
                insertions = re.findall(r"(\d+)\^", read.get_tag('MD'))
                for insertion in insertions:
                    regionRead.insert((int(insertion) - 1), 3)
            '''             
                   
            cigar = read.cigartuples
            cigarCount = 0
            for operation in cigar:
                if operation[0] != 4 and operation[0] != 2:
                    cigarCount += operation[1]
                if operation[0] == 1:
                    #print(regionRead)
                    #print(cigar)
                    #print(read.get_tag('MD'))
                    #print(alignment)
                    #print(cigar)
                    try:
                        regionRead[cigarCount - operation[1]] = 2
                    except:
                        print('cigar index error: ' + str(cigarCount))
                        print('region read length: ' + str(len(regionRead)))
                        print(cigar)
                        for position in regionRead:
                            print(position, end = '')
                        print('\n')
                        print(actualRef)
                        print(read.query_sequence)
                        print(alignment)
                        sys.exit()
                    #print(regionRead)
            #if regionRead[0] == 3:
                #print(read.query_sequence)
                #print(actualRef)
                #for position in regionRead:
                #    print(position, end = '')
                #print('\n')
                #print(cigar)
                #print(read.get_tag('MD'))
                #print(alignment)
                #print('pos: ' + str(read.get_reference_positions()[0]))
                #print('*******')     
            for base in range(len(regionRead)):
                try:
                    if regionRead[base] == 1:             
                        baseErrorSummary[base] += 1
                        errorSummary[base] += 1
                    elif regionRead[base] == 3:
                        delSummary[base] += 1
                        errorSummary[base] += 1
                    elif regionRead[base] ==2:
                        errorSummary[base] += 1
                        insSummary[base] += 1
                except: 
                    print(str(base))
                    print(len(regionRead))
                    sys.exit()
            #if regionRead[0] != 0:          
                #print(read.query_sequence)
                #print(actualRef)
                #for position in regionRead:
                #    print(position, end = '')
                #print('\n')
                #print(cigar)
                #print(read.get_tag('MD'))
                #print(alignment)
                #print('pos: ' + str(read.get_reference_positions()[0]))
                #print('*******')     
            if readCount % 10000 == 0:
                outFile.write('Reads: ' + str(readCount) + '\n')
                delFile.write('Reads: ' + str(readCount) + '\n')
                insFile.write('Reads: ' + str(readCount) + '\n')
                baseFile.write('Reads: ' + str(readCount) + '\n')
                
                for x in range(300):
                    outFile.write(str(x) + '\t' + str(errorSummary[x]) + '\n')
                    delFile.write(str(x) + '\t' + str(delSummary[x]) + '\n') 
                    insFile.write(str(x) + '\t' + str(insSummary[x]) + '\n')
                    baseFile.write(str(x) + '\t' + str(baseErrorSummary[x]) + '\n')
                
        outFile.write('Reads: ' + str(readCount) + '\n')
        delFile.write('Reads: ' + str(readCount) + '\n')
        insFile.write('Reads: ' + str(readCount) + '\n')
        baseFile.write('Reads: ' + str(readCount) + '\n')
                
        for x in range(300):
            outFile.write(str(x) + '\t' + str(errorSummary[x]) + '\n')
            delFile.write(str(x) + '\t' + str(delSummary[x]) + '\n') 
            insFile.write(str(x) + '\t' + str(insSummary[x]) + '\n')
            baseFile.write(str(x) + '\t' + str(baseErrorSummary[x]) + '\n')    
