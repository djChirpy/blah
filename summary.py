import pandas as pd
import argparse


parser = argparse.ArgumentParser()
parser.add_argument('-i', '--inFile', nargs = '+')
args = parser.parse_args()
args.inFile = os.listdir('/seq/vgb/swofford/temp/umass/geneCounts')

clades = pd.read_csv('/seq/vgb/swofford/zoonomia/discovar_file_locations.txt', sep = '\t', encoding='cp1252')
newNames = []
for i in clades['Species'].values:
    newNames.append('_'.join(i.split()))

lookup = {}
with open('/seq/vgb/swofford/zoonomia/conservedGenes/geneNames.txt', 'r') as look:
    for line in look:
        fields = line.strip().split('\t')
        lookup[fields[0]] = fields[1]
        lookup[fields[1]] = fields[0]
        

busco = []
with open('/seq/vgb/swofford/zoonomia/comparative/BUSCO_humanOrthIDsSeparated.out', 'r') as buscIn:
    for line in buscIn:
        name = line.strip().split()[0]
        try:
            busco.append(lookup[name])
        except:
            print(line)


speciesCounts = {}
for file in args.inFile:
    file = '/seq/vgb/swofford/temp/umass/geneCounts/' + file
    with open(file, 'r') as inFile:
        geneCount = 0
        buscoCount = 0
        try:
            name = file[:-13].split('/')[-1]
            #mini = name[:3] + name.split('_')[1][:1].upper() + name.split('_')[1][1:3]
            #if mini in clades['Assembly Name'].values:
            for line in inFile:
                geneCount += 1
                if line.strip() in busco:
                    buscoCount += 1
            speciesCounts[name] = {}
            speciesCounts[name]['busco'] = buscoCount
            speciesCounts[name]['genes'] = geneCount
        except:
            continue
        #lse:
        #    print(mini)
                      
clades.set_index('Species', inplace = True)
counts = pd.DataFrame.from_dict(speciesCounts)            
species = {}
