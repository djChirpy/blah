with open('cf3-20k.intervals', 'r') as intervals, open('cf3-noIntervals.intervals', 'w') as outFile:
    for line in intervals:
        try:
            if chr != line.split(':')[0]:
                outFile.write(chr + ':1-'+str(end))

                chr = line.split(':')[0]
                end = line.split(':')[1].split('-')[1]
        except NameError:
                        
            chr = line.split(':')[0]
            end = line.split(':')[1].split('-')[1]

                        
