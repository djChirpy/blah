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

with pysam.AlignmentFile(args.inFile, 'rb') as bam, open(args.outFile, 'w') as outFile:
#    outFile.write('numDeletions\tdeletedBases\tlongestDeletion\tnumInsertions\tinsertedBases\tlongestInsertion\tnumSubstitutions\n')

    #instantiate array of 0s to track errors per base
    errorSummary = []
    for i in range(204):
            errorSummary.append(0)

    for region in regions:
            #outFile.write(str(region) + '\n')
        print(str(region))
        actualRef = str(refVars[region])
        print(str(len(actualRef)))
        iterRead = bam.fetch(str(region))

        regionRead = []
        
        readCount = 0
        
            
        for read in iterRead:
            readCount += 1
            #md = read.get_tag('MD')
            #cigar = read.cigartuples
            alignment = read.get_aligned_pairs(with_seq=True)
            baseCount = 0
            adjustLeft = 0
            adjustRight = 0
            regionRead = []

            #check each position in the alignment for errors
            for pair in alignment:
                baseCount += 1
                if not isNumber(pair[0]):
                    #insertion 
                    adjustLeft += 1
                    regionRead.append(3)
                elif not isNumber(pair[1]):
                    #deletion
                    adjustRight +=1 
                elif (pair[0] + adjustLeft) == (pair[1] + adjustRight):
                    if str(pair[2]).isupper():
                        regionRead.append(0)
                    else:
                        #mismatch
                        regionRead.append(1)
                
            #print(str(regionRead)) 
            for base in range(204):
                #print('base =' + str(base))
                #print(str(regionRead[base]))                
                try:
                    if regionRead[base] != 0:    
                        errorSummary[base] += 1
                except(IndexError):
                    continue
            print(read.get_seq())
            print(refVars[region])
            for base in regionRead:
                print(base, end="")
            
                        #end of oligo, enterred barcode
    #print(str(errorSummary))
    #print(str(readCount))
        outFile.write(str(region))
        for base in errorSummary:            
            outFile.write(str(errorSummary.index(base)) + '\t' + str(base) + '\n')    #print(str(regionRead))
                
                
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
