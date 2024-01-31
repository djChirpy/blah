import pysam
import sys
import re


def isNumber(num):
    try:
        float(num)
        return True
    except ValueError:
        return False


def buildRef():
    '''Returns a dictionary of ordered oligo variant inserts and names.
    Input is original design file from Xue Li'''
    ref = open(sys.argv[2])
    variants = {}
    for line in ref:
        fields = re.split('\s+', line.rstrip())
        name = str(fields[0]) + '_' + str(fields[1]) + \
            '_' + str(fields[2]) + '_' + str(fields[9])
        var = str(fields[10])
        name = name.replace(',', '^')
        variants[name] = [var]
    ref.close()
    return variants


# keep running tabs on all seen barcodes
barcodes = {}
counter = 0

refVars = buildRef()
regions = []
# get list of all ref sequence names, these will corespond to regions in the
# reference FASTA used to create the bam files.
for name, var in refVars.items():
    print(name)
    regions.append(name)


with pysam.AlignmentFile(sys.argv[1], 'rb') as bam, open(sys.argv[3] + '.out', 'w') as outFile, open(sys.argv[3] + '_summary.out', 'w') as summary, open(sys.argv[3] + '_barcodeTemp.out', 'w') as barcodeTemp, open(sys.argv[3] + '_perfectlyAlignedReads.out', 'w') as perfectAlign:
    # Cycle through reference file, one variant at a time.
    for region in regions:
        outFile.write(str(region) + '\n')
        print(str(region))
        actualRef = str(refVars[region])
        iterRead = bam.fetch(str(region))
        oligoBarcodes = set()
        perfectBarcodes = set()
        reads1 = {}
        reads2 = {}
        reads3 = {}
        oligos = {}
        perfectReads = 0
        imperfectReads = 0

        for read in iterRead:

            if read.is_read1 and read.is_proper_pair and read.is_paired:
                name = read.query_name
                reads1[name] = read
            elif read.is_read2 and read.is_proper_pair and read.is_paired:
                name = read.query_name
                reads2[name] = read

        print(str(len(reads1)) + 'read1s in region: ' + str(region))
        print(str(len(reads2)) + 'read2s in region: ' + str(region))

        for name, read in reads1.items():
            try:
                reads3[name] = [read, reads2[name]]
            except(KeyError):
                print("unMatched read found, skipping")
                continue

        print(str(len(reads3)) + 'paired reads in region: ' + str(region))

        for name, reads in reads3.items():
            md = reads[0].get_tag('MD')
            md2 = reads[1].get_tag('MD')
            if reads[0].is_reverse:
                barcode = reads[0].query_sequence[-20:]

            else:
                barcode = reads[1].query_sequence[-20:]
                print("barcode weirdness happening: " + barcode)
            #consensus = align(reads[0].query_sequence, reads[1].query_sequence)

            

            # trying to redo barcode tracking. Dictionary of sets, where each
            # key is a barcode with a set of corresponding regions (variants)
            # also need a barcodes per region thing, but that can be handled
            # in the summary file, I think.

            oligoBarcodes.add(barcode)
            if isNumber(md) and isNumber(md2) and ((int(md) + int(md2)) >= 220):
                perfect = 0
                perfectBarcodes.add(barcode)
                perfectReads += 1
                barcode2 = reads[0].query_sequence[(int(md) + 38):]
                perfectAlign.write(barcode + '\t' + barcode2 + '\t' + region + '\t' + reads[0].query_sequence + '\t' + reads[1].query_sequence + '\n')
            else:
                perfect = 1
                imperfectReads += 1

            barcodeTemp.write(
                barcode +
                '\t' +
                str(perfect) +
                '\t' +
                str(region) +
                '\n')

            # Can probably simplify this into region line, followed by lines of
            # barcode, perfect, md, md2
            oligos[name] = [barcode, perfect, md, md2]
            #oligos[name] = [region, barcode, perfect, md, md2, reads[0].query_sequence, reads[1].query_sequence, consensus, actualRef]

        summary.write(str(region) +
                      '\t' +
                      str(len(oligoBarcodes)) +
                      '\t' +
                      str(len(perfectBarcodes)) +
                      '\t' +
                      str(perfectReads) +
                      '\t' +
                      str(imperfectReads) +
                      '\n')

        for name, data in oligos.items():

            for field in data:
                outFile.write('\t' + str(field))
            outFile.write('\n')

