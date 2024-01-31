import sys
import csv

with open(sys.argv[1]) as clusters, open(sys.argv[2] + ".tree", 'w') as treeMix:
    cluster = csv.DictReader(
        clusters, fieldnames=('CHR',
                              'SNP',
                              'CLST',
                              'A1',
                              'A2',
                              'MAF',
                              'MAC',
                              'NCHROBS'), delimiter="\s+")

    next(cluster)
    groups = {}

    line = next(cluster)
    print(line)
    snp = line['SNP']
    group = line['CLST']
    allele1 = int(line['MAC'])
    allele2 = int(line['NCHROBS']) - allele1

    print(line)
    print(snp)
    print(group)
    print(allele1)
    print(allele2)
