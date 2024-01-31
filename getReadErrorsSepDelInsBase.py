import pysam
import argparse
import re



parser = argparse.ArgumentParser()
parser.add_argument('-i', '--inFile')
parser.add_argument('-o', '--outFile')
parser.add_argument('-v', '--varFile')
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
            var = str(fields[12])
            name = name.replace(',', '^')
            variants[name] = var
        ref.close()
        return variants


refVars = buildRef()
regions = []
# get list of all ref sequence names, these will corespond to regions in the                                                                                                                                                                                                                                 
# reference FASTA used to create the bam files.                                                                                                                                                                                                                                                              
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
    for i in range(204):
            baseErrorSummary.append(0)
            errorSummary.append(0)
            delSummary.append(0)
            insSummary.append(0)

    for region in regions:
            #outFile.write(str(region) + '\n')
        count += 1
        print(str(region))
        actualRef = str(refVars[region])
        print(str(len(actualRef)))
        iterRead = bam.fetch(str(region))

        regionRead = []
        
        
        
            
        for read in iterRead:
            readCount += 1
            #md = read.get_tag('MD')
            #cigar = read.cigartuples
            alignment = read.get_aligned_pairs(with_seq=True)
            baseCount = 0
            adjustLeft = 0
            adjustRight = 0
            regionRead = []
            firstIns = True
            #print(read.get_reference_sequence())
            #print(read.query_sequence)

            #check each position in the alignment for errors
            for pair in alignment:
                baseCount += 1
                if not isNumber(pair[0]):
                    #deletion 
                    adjustLeft += 1
                    regionRead.append(3)
                    firstIns = True
                elif not isNumber(pair[1]):
                    #insertion
                    adjustRight +=1
                    if firstIns:
                        regionRead.append(2)
                        firstIns = False
                    else:
                        regionRead.append(0)
                    #insertion = True 
                elif (pair[0] + adjustLeft) == (pair[1] + adjustRight):
                    firstIns = True
                    if str(pair[2]).isupper():
                        regionRead.append(0)
                    else:
                        #mismatch
                        regionRead.append(1)
            #if insertion:
            #    print(read.get_reference_sequence())
            #    print(read.query_sequence)
                
            #print(str(regionRead)) 
            for base in range(300):
                #if insertion:
                #    try:
                #        print(regionRead[base], end = '')
                #    except: 
                #        continue
                #print('base =' + str(base))
                    #print(str(regionRead[base]))                
                try:
                    if regionRead[base] == 1:
                        baseErrorSummary[base] += 1
                    elif regionRead[base] == 2:
                        insSummary[base] += 1
                    elif regionRead[base] == 3:
                        delSummary[base] += 1
                        
                    if regionRead[base] != 0:    
                        errorSummary[base] += 1
                except(IndexError):
                    continue
            #if insertion:
            #    print('\n')
#            print(read.get_seq())
#            print(refVars[region])
#            for base in regionRead:
#                print(base, end="")
            
                        #end of oligo, enterred barcode
    #print(str(errorSummary))
    #print(str(readCount))
        outFile.write(str(region) + '\tReads: ' + str(readCount))
        delFile.write(str(region) + '\tReads: ' + str(readCount))
        insFile.write(str(region) + '\tReads: ' + str(readCount))
        baseFile.write(str(region) + '\tReads: ' + str(readCount))
        for base in errorSummary:            
            outFile.write(str(errorSummary.index(base)) + '\t' + str(base) + '\n')    #print(str(regionRead))
        for base in delSummary:            
            delFile.write(str(delSummary.index(base)) + '\t' + str(base) + '\n')
        for base in insSummary:            
            insFile.write(str(insSummary.index(base)) + '\t' + str(base) + '\n')
        for base in baseErrorSummary:            
            baseFile.write(str(baseErrorSummary.index(base)) + '\t' + str(base) + '\n')        
                
  #              and isupper(pair[3]):
  #                  aligned[baseCount] == 0
            
            


'''
    for read in bam:
        try:
            deletions = []
            insertions = []
            deletedBases = 0
            insertedBases = 0
            cigar = read.cigartuples

            md = read.get_tag('MD')
            print(md)
            for item in re.split('[ACTG^]', read.get_tag("MD")):
                print(str(item))
            #print('cigar:\t' + str(cigar) + '\t' + str(md))
        except(TypeError):
            print('errror')

            substitution = 0
            if '^' not in str(md):
                for character in str(md):
                    if character == 'A' or character == 'C' or character == 'G' or character == 'T':
                        substitution += 1
            outFile.write(str(len(deletions)) + '\t' + str(deletedBases) + '\t' + str(largestDel) + '\t' + str(len(insertions)) + '\t' + str(insertedBases) + '\t' + str(largestInsertion) + '\t' + str(substitution) + '\n')
        except(TypeError):
            continue
#        print(str(tags))
#        print('Cigar = ' + str(cigar))
#        print('MD = ' + str(md))
'''
