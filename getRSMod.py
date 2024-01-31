import sys
from csv import DictReader, DictWriter

with open(sys.argv[1]) as rsBim, open(sys.argv[2]) as oldBim, open(sys.argv[3], 'w') as newBim:

    rsReader = DictReader(
        rsBim,
        fieldnames=(
            'CHR',
            'name',
            'noIdea',
            'pos',
            'ref',
            'alt'), delimiter='\t')

    bimReader = DictReader(
        oldBim,
        fieldnames=(
            'CHR',
            'name',
            'noIdea',
            'pos',
            'ref',
            'alt'), delimiter='\t')

    bimWriter = DictWriter(
        newBim,
        fieldnames=(
            'CHR',
            'name',
            'noIdea',
            'pos',
            'ref',
            'alt'), delimiter='\t')

    rs = next(rsReader)

    for row in bimReader:
        try:
            
            while int(rs['CHR']) < int(row['CHR']):
                rs = next(rsReader)

            while int(
                rs['CHR']) == int(
                row['CHR']) and int(
                rs['pos']) < int(
                    row['pos']):
                rs = next(rsReader)

            if int(
                row['CHR']) == int(
                rs['CHR']) and int(
                row['pos']) == int(
                    rs['pos']):
                bimWriter.writerow(rs)
                continue
            
            else:
                bimWriter.writerow(row)
                
            #elif int(row['CHR']) < int(rs['CHR']):
            #    bimWriter.writerow(row)
            #    continue

            #elif int(row['CHR']) == int(rs['CHR']) and int(row['pos']) < int(rs['pos']):
            #    bimWriter.writerow(row)
            #    continue
            
        except:
            print(row)
            bimWriter.writerow(row)
