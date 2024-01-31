import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--inFile', nargs='+')

args = parser.parse_args()

for file in args.inFile:
    type = file.split('_')[0]
    with open(file, 'r') as inFile, open(file[:-4] + '_updated.txt', 'w') as outFile:
        for line in inFile:
            if 'prediction' in line:
                outFile.write('correct\t' + line)
            else:
                pred, label = line.split()[0], line.split()[1]
                if pred == 'other':
                    if label != type:
                        outFile.write('yes\t' + line)
                    else:
                        outFile.write('no\t' + line)
                else:
                    if label == pred:
                        outFile.write('yes\t' + line)
                    else:
                        outFile.write('no\t' + line)
                