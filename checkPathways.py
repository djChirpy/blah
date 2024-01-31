import pandas as pd
import numpy as np
import re
import csv
import gzip
import os
import sys

checkFile = 'pathwayGenesToCheckByGroup.txt' #copy of tab from Kate's google doc https://docs.google.com/spreadsheets/d/1jYzQKtV2D5PE--PI3tCQjIkIkiN_eub2WUarUfNyI04/edit#gid=0
refFile = '/seq/vgb/swofford/dogHumanOrthologs/old/feb2/ncbiGeneNameAdded.bed'
pathFile = '/seq/vgb/swofford/ref/pathways/C2CuratedGeneSetsGeneNames.gmt.txt'

''' Class defining a single gene ortholog relationship between dog and human, to facilitate bi-directional lookups later'''
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
        for line in pathRef:
            genes = line.strip().split('\t')
            pathName = genes[0]
            pathwayList[pathName] = []
            for gene in genes[2:]:
                pathwayList[pathName].append(gene)
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

'''get the corresponding dog stableGeneID from human stable geneID'''
def getDogGeneID(geneID):
    try:
        return reference[geneID].dogGeneID
    except (KeyError):
        return 'NA'

'''get the corresponding Human stable gene ID from the stable dog geneID'''
def getHumanGeneID(geneID):
    try:
        return reference[geneID].geneID
    except (KeyError):
        return 'NA'

'''get ortholog object for geneID'''
def getOrtholog(geneID):
    try:
        return reference[geneID]
    except (KeyError):
        return 'NA'

ref = {}
pathwayList = {}
genePath = {}    


df = pd.read_csv('full_updatedLabels_preProcessed_lowPerformingLabelsRemoved.txt', sep = '\t')

if not os.path.exists('pathwayGenes'):
    os.mkdir('pathwayGenes')
    
reference = buildRef(refFile, pathFile)
lookup = {}

with open(checkFile, 'r') as inFile:
    for line in inFile:
        fields = line.strip().split()
        if fields[0] == 'Pathway':
            continue
        pathway = fields[0]
        groups = []
        for i in range(len(fields) - 1):
            groups.append(fields[i + 1])
        lookup[pathway] = groups

for i in lookup.keys():
    keepCol = ['ID', 'label', 'species']
    genes = set()
    for gene in pathwayList[i]:
        for x in list(getHumanGeneID(gene)):
            if 'ENSG' in x and x in df.columns:
                genes.add(x)
                genes.add(x + '_count')
    genes = list(genes)
    groups = lookup[i]
    keepCol = np.append(keepCol, genes)
    keep = df[df['label'].isin(groups)][keepCol]
    meanKeep = keep.groupby('label').mean()
    keep.to_csv('pathwayGenes/' + i + '_' + '_'.join(groups) + '.txt', sep = '\t', index = False)
    meanKeep.to_csv('pathwayGenes/' + i + '_' + '_'.join(groups) + '_mean.txt', sep = '\t', index = False)
    
    


