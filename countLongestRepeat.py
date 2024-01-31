import sys
import re
# takes in oligo order list from XiuLi and spits out fasta


def buildRef(ref):
    '''Returns a dictionary of ordered oligo variant inserts and names'''
    variants = {}
    for line in ref:
        fields = re.split('\t', line.rstrip())
        name = str(fields[0]) + '_' + str(fields[1]) + \
            '_' + str(fields[2]) + '_' + str(fields[9])
        var = str(fields[10])
        name = name.replace(',', '^')
        variants[name] = var
        print(str(len(variants)))
    ref.close()
    return variants


def writeFasta(reads):
    for k, v in reads.items():
        #print('name: ' + k)
        #print('oligo: ' + v)
        outFile.write('>' + k + '\n' + v + '\n')

with open(sys.argv[1], 'rt') as ref, open(sys.argv[2], 'w') as outFile:
    fastaVars = buildRef(ref)
    #writeFasta(fastaVars)
    for k, v in fastaVars.items():
        prev = ''
        highestCount = 0
        count = 0
        GC = 0
        for base in v:
            if base == 'G' or base == 'C':
                GC += 1
            if base == prev:
                count += 1
            else:
                if count > highestCount:
                    highestCount = count
                count = 1
                prev = base
        if count > highestCount:
            highestCount = count
        outFile.write(k + '\t' + str(highestCount) + '\t' + str(GC/len(v)) + '\n')
                
            
