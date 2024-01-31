# 538447106 barcodes
# testing: 1970 saved with full script in 10:40
# testing: 784 saved with no testing singletons 

import sys
import re
from fuzzywuzzy import fuzz, process


def hashingOrder(k):
    '''hashing the reference oligos into bins based on their starting seq.
    assumption is that given this is fairly early in read, quality should
    be good? Returns dictionary of lists, where each list is a bin of ref
    oligos'''
    table = {}
    table['AAAA'] = []
    table['AAAC'] = []
    table['AAAG'] = []
    table['AAAT'] = []
    table['AACA'] = []
    table['AACC'] = []
    table['AACG'] = []
    table['AACT'] = []
    table['AAGA'] = []
    table['AAGC'] = []
    table['AAGG'] = []
    table['AAGT'] = []
    table['AATA'] = []
    table['AATC'] = []
    table['AATG'] = []
    table['AATT'] = []
    table['ACAA'] = []
    table['ACAC'] = []
    table['ACAG'] = []
    table['ACAT'] = []
    table['ACCA'] = []
    table['ACCC'] = []
    table['ACCG'] = []
    table['ACCT'] = []
    table['ACGA'] = []
    table['ACGC'] = []
    table['ACGG'] = []
    table['ACGT'] = []
    table['ACTA'] = []
    table['ACTC'] = []
    table['ACTG'] = []
    table['ACTT'] = []
    table['AGAA'] = []
    table['AGAC'] = []
    table['AGAG'] = []
    table['AGAT'] = []
    table['AGCA'] = []
    table['AGCC'] = []
    table['AGCG'] = []
    table['AGCT'] = []
    table['AGGA'] = []
    table['AGGC'] = []
    table['AGGG'] = []
    table['AGGT'] = []
    table['AGTA'] = []
    table['AGTC'] = []
    table['AGTG'] = []
    table['AGTT'] = []
    table['ATAA'] = []
    table['ATAC'] = []
    table['ATAG'] = []
    table['ATAT'] = []
    table['ATCA'] = []
    table['ATCC'] = []
    table['ATCG'] = []
    table['ATCT'] = []
    table['ATGA'] = []
    table['ATGC'] = []
    table['ATGG'] = []
    table['ATGT'] = []
    table['ATTA'] = []
    table['ATTC'] = []
    table['ATTG'] = []
    table['ATTT'] = []
    table['CAAA'] = []
    table['CAAC'] = []
    table['CAAG'] = []
    table['CAAT'] = []
    table['CACA'] = []
    table['CACC'] = []
    table['CACG'] = []
    table['CACT'] = []
    table['CAGA'] = []
    table['CAGC'] = []
    table['CAGG'] = []
    table['CAGT'] = []
    table['CATA'] = []
    table['CATC'] = []
    table['CATG'] = []
    table['CATT'] = []
    table['CCAA'] = []
    table['CCAC'] = []
    table['CCAG'] = []
    table['CCAT'] = []
    table['CCCA'] = []
    table['CCCC'] = []
    table['CCCG'] = []
    table['CCCT'] = []
    table['CCGA'] = []
    table['CCGC'] = []
    table['CCGG'] = []
    table['CCGT'] = []
    table['CCTA'] = []
    table['CCTC'] = []
    table['CCTG'] = []
    table['CCTT'] = []
    table['CGAA'] = []
    table['CGAC'] = []
    table['CGAG'] = []
    table['CGAT'] = []
    table['CGCA'] = []
    table['CGCC'] = []
    table['CGCG'] = []
    table['CGCT'] = []
    table['CGGA'] = []
    table['CGGC'] = []
    table['CGGG'] = []
    table['CGGT'] = []
    table['CGTA'] = []
    table['CGTC'] = []
    table['CGTG'] = []
    table['CGTT'] = []
    table['CTAA'] = []
    table['CTAC'] = []
    table['CTAG'] = []
    table['CTAT'] = []
    table['CTCA'] = []
    table['CTCC'] = []
    table['CTCG'] = []
    table['CTCT'] = []
    table['CTGA'] = []
    table['CTGC'] = []
    table['CTGG'] = []
    table['CTGT'] = []
    table['CTTA'] = []
    table['CTTC'] = []
    table['CTTG'] = []
    table['CTTT'] = []
    table['GAAA'] = []
    table['GAAC'] = []
    table['GAAG'] = []
    table['GAAT'] = []
    table['GACA'] = []
    table['GACC'] = []
    table['GACG'] = []
    table['GACT'] = []
    table['GAGA'] = []
    table['GAGC'] = []
    table['GAGG'] = []
    table['GAGT'] = []
    table['GATA'] = []
    table['GATC'] = []
    table['GATG'] = []
    table['GATT'] = []
    table['GCAA'] = []
    table['GCAC'] = []
    table['GCAG'] = []
    table['GCAT'] = []
    table['GCCA'] = []
    table['GCCC'] = []
    table['GCCG'] = []
    table['GCCT'] = []
    table['GCGA'] = []
    table['GCGC'] = []
    table['GCGG'] = []
    table['GCGT'] = []
    table['GCTA'] = []
    table['GCTC'] = []
    table['GCTG'] = []
    table['GCTT'] = []
    table['GGAA'] = []
    table['GGAC'] = []
    table['GGAG'] = []
    table['GGAT'] = []
    table['GGCA'] = []
    table['GGCC'] = []
    table['GGCG'] = []
    table['GGCT'] = []
    table['GGGA'] = []
    table['GGGC'] = []
    table['GGGG'] = []
    table['GGGT'] = []
    table['GGTA'] = []
    table['GGTC'] = []
    table['GGTG'] = []
    table['GGTT'] = []
    table['GTAA'] = []
    table['GTAC'] = []
    table['GTAG'] = []
    table['GTAT'] = []
    table['GTCA'] = []
    table['GTCC'] = []
    table['GTCG'] = []
    table['GTCT'] = []
    table['GTGA'] = []
    table['GTGC'] = []
    table['GTGG'] = []
    table['GTGT'] = []
    table['GTTA'] = []
    table['GTTC'] = []
    table['GTTG'] = []
    table['GTTT'] = []
    table['TAAA'] = []
    table['TAAC'] = []
    table['TAAG'] = []
    table['TAAT'] = []
    table['TACA'] = []
    table['TACC'] = []
    table['TACG'] = []
    table['TACT'] = []
    table['TAGA'] = []
    table['TAGC'] = []
    table['TAGG'] = []
    table['TAGT'] = []
    table['TATA'] = []
    table['TATC'] = []
    table['TATG'] = []
    table['TATT'] = []
    table['TCAA'] = []
    table['TCAC'] = []
    table['TCAG'] = []
    table['TCAT'] = []
    table['TCCA'] = []
    table['TCCC'] = []
    table['TCCG'] = []
    table['TCCT'] = []
    table['TCGA'] = []
    table['TCGC'] = []
    table['TCGG'] = []
    table['TCGT'] = []
    table['TCTA'] = []
    table['TCTC'] = []
    table['TCTG'] = []
    table['TCTT'] = []
    table['TGAA'] = []
    table['TGAC'] = []
    table['TGAG'] = []
    table['TGAT'] = []
    table['TGCA'] = []
    table['TGCC'] = []
    table['TGCG'] = []
    table['TGCT'] = []
    table['TGGA'] = []
    table['TGGC'] = []
    table['TGGG'] = []
    table['TGGT'] = []
    table['TGTA'] = []
    table['TGTC'] = []
    table['TGTG'] = []
    table['TGTT'] = []
    table['TTAA'] = []
    table['TTAC'] = []
    table['TTAG'] = []
    table['TTAT'] = []
    table['TTCA'] = []
    table['TTCC'] = []
    table['TTCG'] = []
    table['TTCT'] = []
    table['TTGA'] = []
    table['TTGC'] = []
    table['TTGG'] = []
    table['TTGT'] = []
    table['TTTA'] = []
    table['TTTC'] = []
    table['TTTG'] = []
    table['TTTT'] = []

    for line in k:
        table[line[:4]].append(line)

    return(table)


def hammingDist(seq1, seq2):
    '''Determines the Hamming distance between two strings. Equivalent to the
    minimum number of character substitutions required to change from one
    string to the other. Takes two strings as input, and returns an integer
    value for the Hamming distance'''
    return sum(ch1 != ch2 for ch1, ch2 in zip(seq1, seq2))


def isNumber(num):
    try:
        float(num)
        return True
    except ValueError:
        return False


def buildRef(ref):
    '''Returns a dictionary of ordered oligo variant inserts and names'''
    variants = {}
    for line in ref:
        fields = re.split('\s+', line.rstrip())
        name = str(fields[0]) + '_' + str(fields[1]) + \
            '_' + str(fields[2]) + '_' + str(fields[9])
        var = str(fields[10])
        name = name.replace(',', '^')
        variants[name] = var
        variants[var] = name
    ref.close()
    return variants


def revComplement(seq):
    '''Returns the reverse complement of a DNA sequence passed to it.
    Accepts a string as input, and returns a string. Valid bases for
    input are "ATCGW"'''
    compSeq = ''
    codes = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'W': 'W', 'N': 'N'}
    for letter in seq:
        compSeq += codes[letter]
    return(compSeq[::-1])


def writeBarcode(
        barcode,
        match,
        barcodeCount,
        exactMatch,
        unexactMatch,
        oligoList=[]):
    '''writes a line to the output file consisting of the barcodes seq,
    the matching reference oligo, the total number of reads for the barcode,
    the number of reads that matched the reference exactly, the number of
    reads that didn't match the reference exactly, and the sequences for
    the reads that did not match exactly, in case further review required'''

    outPut.write(match + 
        '\t' +
        barcode +
        '\t' +
        str(barcodeCount) +
        '\t' +
        str(exactMatch) +
        '\t' +
        str(unexactMatch) +
        '\t' +
        refVars[match])
    for oligo in oligoList:
        outPut.write('\t' + oligo)
    outPut.write('\n')


with open(sys.argv[1], 'rt') as reads, open(sys.argv[2], 'r') as reference, open(sys.argv[3], 'w') as outPut:

    barcodeCount = 0
    barcodeOligos = 0
    storedBarcode = 'first'
    namedSeq = set()
    unNamedSeq = set()
    counter = 0
    exact = 0
    bad = False
    refList = []

    # build a list of reference oligos for use with fuzzywuzzy
    refVars = buildRef(reference)
    for k, v in refVars.items():
        if 'chr' not in v:
            refList.append(v)

    hashRef = hashingOrder(refList)

    for line in reads:
        namedSeq = set()
        unNamedSeq = set()
        barcodeCount = 0
        bad = False
        testFuzz = 0
        
        try:
            counter += 1
            fields = re.split('\t', line)
            barcode = fields[0].strip()
            barcodeCount = int(fields[2])
            barcodeOligos = int(fields[1])
        except ValueError:
            print(line)
            continue
        exactMatch = 0
        unexactMatch = 0
        # parse fields into individual oligos
        for field in fields[3:]:
            seq = field.strip()

            if isNumber(seq):
                seqCount = seq
            # sort oligos into those that match design and those that don't
            elif "chr" in seq:
                namedSeq.add(seq)
                exactMatch += 1
            else:
                unNamedSeq.add(seq)
                unexactMatch += 1

        if len(namedSeq) > 1:
            continue

        # If there's exactly one ref matching oligo, great!
        if len(namedSeq) == 1 and len(unNamedSeq) == 0:
            writeBarcode(
                barcode,
                list(namedSeq)[0],
                barcodeCount,
                exactMatch,
                unexactMatch)

        # If there's a ref matching seq, and some not, check to see if non
        # matching oligos are close to matching (fuzzywuzzy) If close,
        # we'll assume they match, and call the barcode "good", but also
        # keep track of the exact sequencing results
        elif len(namedSeq) == 1 and len(unNamedSeq) >= 1:
            ref = refVars[list(namedSeq)[0]]
            for oligo in list(unNamedSeq):
                testFuzz = fuzz.ratio(ref, oligo)
                if testFuzz < 95:
                    bad = True
                    break

            if not bad:
                writeBarcode(
                    barcode,
                    list(namedSeq)[0],
                    barcodeCount,
                    exactMatch,
                    unexactMatch,
                    unNamedSeq)

        # if there's more than one non-ref-matching oligos, first check to
        # see if they closely match each other. If so, pick one and see if
        # it closely matches any of the reference oligos
        elif len(unNamedSeq) > 1:
            #print('multiple unmatched slow if match eachother ')
            seed = list(unNamedSeq)[0]
            for oligo in unNamedSeq:
                testFuzz = fuzz.ratio(seed, oligo)
                if testFuzz < 95:
                    bad = True
                    break
            if not bad:
                try:
                #print('they matched eachother, slow search')
                    testSeq = list(unNamedSeq)[0]
                    bestMatch = process.extract(
                        testSeq,
                        hashRef[
                            testSeq[
                                :4]],
                        limit=2)
                # print(bestMatch)

                    if bestMatch[0][1] > bestMatch[1][1] and bestMatch[0][1] >= 95:
                    #print('writing fuzzy best match')
                        writeBarcode(
                            barcode,
                            refVars[
                                bestMatch[0][0]],
                            barcodeCount, exactMatch, unexactMatch,
                            unNamedSeq)
                

                    elif bestMatch[0][1] == bestMatch[1][1] and bestMatch[0][1] >= 95:
                        hamOne = hammingDist(bestMatch[0][0], testSeq)
                        hamTwo = hammingDist(bestMatch[1][0], testSeq)
                        if hamOne < hamTwo:
                        #print('writing hamOne')
                            writeBarcode(
                                barcode,
                                refVars[
                                    bestMatch[0][0]],
                                barcodeCount, exactMatch, unexactMatch, unNamedSeq)
                        elif hamTwo < hamOne:
                        #print('writing hamTwo')
                            writeBarcode(
                                barcode,
                                refVars[
                                    bestMatch[1][0]],
                                barcodeCount, exactMatch, unexactMatch, unNamedSeq)
                except KeyError:
                    continue
               

        namedSeq = set()
        unNamedSeq = set()
        barcodeCount = 0
        bad = False
        testFuzz = 0
