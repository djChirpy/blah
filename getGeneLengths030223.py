import pickle
import subprocess

class ortholog:
    def __init__(self, dogGeneName = '', humanGeneName = set(), dogGeneID = set(), geneID = set(), pathways = set()):
        self.dogGeneName = dogGeneName
        self.humanGeneName = humanGeneName
        self.dogGeneID = dogGeneID
        self.geneID = geneID
        self.pathways = pathways


'''builds a reference table lookup between human/dog orthologous genes, parsed from an Ensembl table. File must contain the below columns, extra columns are fine and will be ignored.'''
def buildRef(refFile, pathways):
    with open(refFile, 'r') as refFile, open(pathways, 'r') as pathRef:
        ref = {}
        pathwayList = {}
        genePath = {}
        for line in pathRef:
            genes = line.strip().split('\t')
            pathName = genes[0]
            for gene in genes[2:]:
                try:
                    if 'http' in gene:
                        print(gene)
                    genePath[gene].append(pathName)
                except KeyError:
                    genePath[gene] = [pathName]
        print(str(len(genePath)) + ' genes with pathways found.')                 
        for line in refFile:
            humEns = set()
            dogEns = set()
            humanGeneName = set()
            pathwaysEffected = set()
            if 'Ensembl' in line or 'NCBI' in line:
                fields = line.strip().split()
                meta = fields[3].split(':')
                if meta[1] == '' or meta[3] == '':
                    continue
                dogEnsembl = meta[1].strip()
                if dogEnsembl != '':
                    dogEns.add(dogEnsembl)
                geneID = meta[3].strip()
                if geneID != '':
                    humEns.add(geneID)
                geneName = meta[2].strip()
                if geneName != '':
                    humanGeneName.add(geneName)
                dogGeneName = meta[0].strip()
                ncbiDogID = meta[6].strip()
                if ncbiDogID != '':
                    dogEns.add(ncbiDogID)
                ncbiHumID = meta[5].strip()
                if ncbiHumID != '':
                    humEns.add(ncbiHumID)
                ncbiHumName = meta[9].strip()
                if ncbiHumName != '':
                    humanGeneName.add(ncbiHumName)
                pathwaysEffected = set()
                for i in humanGeneName:
                    try:
                        pathwaysEffected.update(set(genePath[i]))
                    except KeyError:
                        continue                
                orth = ortholog(dogGeneName, humanGeneName, dogEns, humEns, pathways = pathwaysEffected)
                for id in humEns:
                    try:
                        ref[id].pathways.update(pathwaysEffected)
                        ref[id].geneID.update(humEns)
                        ref[id].humanGeneName.update(humanGeneName)
                    except KeyError:
                        ref[id] = orth
                for name in humanGeneName:
                    if name != '' and name != 'NA':
                        try:
                            ref[name].pathways.update(pathwaysEffected)
                            ref[name].geneID.update(humEns)
                            ref[name].humanGeneName.update(humanGeneName)
                        except KeyError:
                            ref[name] = orth 
                for id in dogEns:
                    try:
                        ref[id].pathways.update(pathwaysEffected)
                        ref[id].geneID.update(humEns)
                        ref[id].humanGeneName.update(humanGeneName)
                    except KeyError:
                        ref[id] = orth    
        print('Finished building Reference. ' + str(len(ref)/3) + ' orthologs recorded.')
        return ref


reference = reference = buildRef('/seq/vgb/swofford/dogHumanOrthologs/old/feb2/ncbiGeneNameAdded.bed', '/seq/vgb/swofford/ref/pathways/C2CuratedGeneSetsGeneNames.gmt.txt')

with open('genes.pickle', 'rb') as inFile:
    genes = pickle.load(inFile)

#with open('genes.pickle', 'rb') as inFile:
#    genes = pickle.load(inFile)

canFam3 = {}
hg38 = {}
hg19 = {}
hg19Merged = {}
hg38Merged = {}
canFam3Merged = {}

with open('hg19CDSUniq.bed', 'r') as inFile:
    for line in inFile:
        fields = line.strip().split('\t')
        try:
            hg19[fields[3].split(';')[0]].append(line)
        except KeyError:
            hg19[fields[3].split(';')[0]] = [line]

for i in hg19.keys():
    with open('temp.bed', 'w') as outFile:
        for line in hg19[i]:
            outFile.write(line)
    with open('tempSorted.bed', 'w') as outFile:
        subprocess.run(['sort','-k1,1','-k2,2n','temp.bed'], stdout=outFile)
    with open('tempMerged.bed', 'w') as outFile:
        subprocess.run(['bedtools','merge','-c','4','-o','distinct', '-i', 'tempSorted.bed'], stdout=outFile)
    with open('tempMerged.bed', 'r') as inFile:
        for line in inFile:
            fields = line.strip().split('\t')
            try:
                hg19Merged[fields[3].split(';')[0]].append(int(fields[2]) - int(fields[1]))
            except KeyError:
                hg19Merged[fields[3].split(';')[0]] = [int(fields[2]) - int(fields[1])]
    subprocess.run(['rm','temp.bed'])
    subprocess.run(['rm','tempSorted.bed'])
    subprocess.run(['rm','tempMerged.bed'])

with open('hg38CDSUniq.bed', 'r') as inFile:
    for line in inFile:
        fields = line.strip().split('\t')
        try:
            hg38[fields[3].split(';')[0]].append(line)
        except KeyError:
            hg38[fields[3].split(';')[0]] = [line]

for i in hg38.keys():
    with open('temp.bed', 'w') as outFile:
        for line in hg38[i]:
            outFile.write(line)
    with open('tempSorted.bed', 'w') as outFile:
        subprocess.run(['sort','-k1,1','-k2,2n','temp.bed'], stdout=outFile)
    with open('tempMerged.bed', 'w') as outFile:
        subprocess.run(['bedtools','merge','-c','4','-o','distinct', '-i', 'tempSorted.bed'], stdout=outFile)
    with open('tempMerged.bed', 'r') as inFile:
        for line in inFile:
            fields = line.strip().split('\t')
            try:
                hg38Merged[fields[3].split(';')[0]].append(int(fields[2]) - int(fields[1]))
            except KeyError:
                hg38Merged[fields[3].split(';')[0]] = [int(fields[2]) - int(fields[1])]
    subprocess.run(['rm','temp.bed'])
    subprocess.run(['rm','tempSorted.bed'])
    subprocess.run(['rm','tempMerged.bed'])
    
    
with open('canFam3CDSUniq.bed', 'r') as inFile:
    for line in inFile:
        fields = line.strip().split('\t')
        try:
            canFam3[fields[3].split(';')[0]].append(line)
        except KeyError:
            canFam3[fields[3].split(';')[0]] = [line]

for i in canFam3.keys():
    with open('temp.bed', 'w') as outFile:
        for line in canFam3[i]:
            outFile.write(line)
    with open('tempSorted.bed', 'w') as outFile:
        subprocess.run(['sort','-k1,1','-k2,2n','temp.bed'], stdout=outFile)
    with open('tempMerged.bed', 'w') as outFile:
        subprocess.run(['bedtools','merge','-c','4','-o','distinct', '-i', 'tempSorted.bed'], stdout=outFile)
    with open('tempMerged.bed', 'r') as inFile:
        for line in inFile:
            fields = line.strip().split('\t')
            try:
                canFam3Merged[fields[3].split(';')[0]].append(int(fields[2]) - int(fields[1]))
            except KeyError:
                canFam3Merged[fields[3].split(';')[0]] = [int(fields[2]) - int(fields[1])]
    subprocess.run(['rm','temp.bed'])
    subprocess.run(['rm','tempSorted.bed'])
    subprocess.run(['rm','tempMerged.bed'])            
        
            
            
with open('/seq/vgb/swofford/temp/geneSizesISH.txt', 'w') as outFile:
    outFile.write('ENS\thg19\thg38\tcanFam3\tname\n')
    for i in genes:
        try:
            hum19 = 0
            hum38 = 0
            dog = 0
            for x in reference[i].geneID:
                try:
                    hum19 += sum(hg19Merged[x])
                except KeyError:
                    pass
            for x in reference[i].geneID:
                try:
                    hum38 += sum(hg38Merged[x])
                except KeyError:
                    pass
            for x in reference[i].dogGeneID:
                try:
                    dog += sum(canFam3Merged[x])
                except KeyError:
                    pass      
            outFile.write('\t'.join([i,str(hum19),str(hum38), str(dog), '@'.join(reference[i].humanGeneName)]) + '\n')
        except KeyError:
            continue                                               
        
        
