import numpy as np
import argparse
import matplotlib.pyplot as plt
import shap
import re
import csv
import gzip
import os
import sys
import seaborn as sns
import pandas as pd
import PIL
from PIL import Image

parser = argparse.ArgumentParser()
parser.add_argument('-hc', '--humanCancer')
parser.add_argument('-dc', '--dogCancer')
parser.add_argument('-i', '--inFile', default = '/seq/vgb/swofford/R01/convertingToSKLInput/testVCFs/artifactSigsRemovedNewlineFixed_080522/artifactSigsRemovedNewlineFixed_080522_featureValues_old.txt')
parser.add_argument('-if', '--inFolder', default = '/seq/vgb/swofford/R01/convertingToSKLInput/testVCFs/artifactSigsRemovedNewlineFixed_080522/')
args = parser.parse_args()



refFile = '/seq/vgb/swofford/dogHumanOrthologs/old/feb2/ncbiGeneNameAdded_zeroRemoved_030623.bed'
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

reference = buildRef(refFile, pathFile)

df = pd.read_csv(args.inFile, sep = '\t')

'''
    what features drive the classification of dog BLSA and dog TLSA as PRAD. These are the only significant matches. 
'''

shap = pd.read_csv('/seq/vgb/swofford/R01/convertingToSKLInput/testVCFs/artifactSigsRemovedNewlineFixed_080522/shapley/PRAD_averageShapley.txt', sep = '\t')


dogShap = shap[shap['label'].isin(['dogBLSA','dogTLSA'])]
dogPRAD = dogShap[dogShap['pred'] == 'PRAD']
feature = dogPRAD.columns[4:]

dogBLSA = dogPRAD[dogPRAD['label'] == 'dogBLSA'][features]
dogTLSA = dogPRAD[dogPRAD['label'] == 'dogTLSA'][features]
topDogBLSA = dogBLSA.reindex(dogBLSA.max(numeric_only = True).sort_values(ascending = False).index, axis = 1).columns[:5]
topDogTLSA = dogTLSA.reindex(dogTLSA.max(numeric_only = True).sort_values(ascending = False).index, axis = 1).columns[:5]

'''
    all dog cancers but melanoma and cmt clump more with some human classes than would be expected by chance. 
    Figure out the top 3 human classes by pred probability for all of these dog classes, and then determine the features driving the connections
'''
def image_grid(imgs, rows, cols):
    len(imgs) == rows*cols
    w, h = imgs[0].size
    grid = Image.new('RGB', size=(cols*w, rows*h))
    grid_w, grid_h = grid.size    
    for i, img in enumerate(imgs):
        grid.paste(img, box=(i%cols*w, i//cols*h))
    return grid


dogClasses = ['dogBLSA','dogGlioma','dogHSA','dogOSA','dogTLSA']
labels = pd.unique(df[df['species'] == 'Human']['label'])
probs = pd.read_csv('/seq/vgb/swofford/R01/convertingToSKLInput/testVCFs/artifactSigsRemovedNewlineFixed_080522/artifactSigsRemovedNewlineFixed_080522_probs.txt', sep = '\t')

#making a copy of df wherein any 'count' features are scaled to 0-1 so they can be plotted along with signatures and binary features

normCount = df.copy(deep = True)
countFeatures = ['ID']
cf = []
for i in normCount.columns:
    if 'ount' in i:
        countFeatures.append(i)
        cf.append(i)
countDf = df[countFeatures].copy(deep = True)
normCount.drop(cf, axis = 1, inplace = True)
normCount.set_index('ID', inplace = True)
countDf.set_index('ID', inplace = True)
countDf = (countDf - countDf.min())/(countDf.max()-countDf.min())
normCount = normCount.merge(countDf, left_on='ID', right_on='ID')



dogClump = {}
for i in dogClasses:
    dogProb = probs[probs['label'] == i]
    topProbs = pd.DataFrame(dogProb['pred'].value_counts()[:3])
    topProbs.reset_index(inplace = True)
    topProbs.columns = ['pred','count']
    featdf = df[df['label'] == i]
    featdf.set_index('ID', inplace = True)
    imgs = []
    for p in topProbs['pred'].values:
        humFeat = df[df['label'] == p]
        humFeat.set_index('ID', inplace = True)
        otherFeat = df[-df['label'].isin([i,p])]
        otherFeat.set_index('ID', inplace = True)
        shapdf = pd.read_csv('/seq/vgb/swofford/R01/convertingToSKLInput/testVCFs/artifactSigsRemovedNewlineFixed_080522/shapley/' + p + '_averageShapley.txt', sep = '\t')
        humShap = shapdf[shapdf['label'] == p]
        humShap.set_index('ID', inplace = True)
        otherShap = shapdf[-shapdf['label'].isin([i, p])]
        otherShap.set_index('ID', inplace = True)
        shapdf = shapdf[shapdf['label'] == i]
        shapdf.set_index('ID', inplace = True)
        topShap =  shapdf[shapdf.columns[3:]].reindex( shapdf[shapdf.columns[3:]].max(numeric_only = True).sort_values(ascending = False).index, axis = 1).columns[:5].values
        dogClump[i] = (p,topShap)
        
        plt.clf()
        shap.summary_plot(shapdf[topShap].to_numpy(), features = featdf[topShap].to_numpy(), feature_names = topShap, show=True, sort = False, plot_type = 'violin')
        plt.savefig('/seq/vgb/swofford/temp/clumpiness/' + i + '_' + p + '.png')
        im = PIL.Image.open('/seq/vgb/swofford/temp/clumpiness/' + i + '_' + p + '.png')
        imgs.append(im)
        plt.clf()
   
        #avg = df[topShap].groupby('label').mean()
        #plt.clf()
        shap.summary_plot(humShap[topShap].to_numpy(), features = humFeat[topShap].to_numpy(), feature_names = topShap, show=True, sort = False, plot_type = 'violin')
        plt.savefig('/seq/vgb/swofford/temp/clumpiness/' + p + '_' + i + '.png')
        im = PIL.Image.open('/seq/vgb/swofford/temp/clumpiness/' + p + '_' + i + '.png')
        imgs.append(im)
        plt.clf()
        shap.summary_plot(otherShap[topShap].to_numpy(), features = otherFeat[topShap].to_numpy(), feature_names = topShap, show=True, sort = False, plot_type = 'violin')
        plt.savefig('/seq/vgb/swofford/temp/clumpiness/other' + '_' + p + '.png')
        im = PIL.Image.open('/seq/vgb/swofford/temp/clumpiness/other' + '_' + p + '.png')
        imgs.append(im)
        
        dogFeat = normCount[normCount['label'] == i][topShap]
        dogFeat = dogFeat.unstack().reset_index()
        dogFeat.columns = ['feature','sample','value']
        plt.clf()
        sns.violinplot(data = dogFeat, x = 'value', y = 'feature', cut = 0, scale = 'count')
        plt.xlim(0,1)
        plt.tight_layout()
        plt.savefig('/seq/vgb/swofford/temp/clumpiness/feature_' + i + '.png')
        im = PIL.Image.open('/seq/vgb/swofford/temp/clumpiness/feature_' + i + '.png')
        imgs.append(im)
        
        humanFeat = normCount[normCount['label'] == p][topShap]
        humanFeat = humanFeat.unstack().reset_index()
        humanFeat.columns = ['feature', 'sample', 'value']
        plt.clf()
        sns.violinplot(data = humanFeat, x = 'value', y = 'feature', cut = 0, scale = 'count')
        plt.xlim(0,1)
        plt.tight_layout()
        plt.savefig('/seq/vgb/swofford/temp/clumpiness/feature_' + p + '.png')
        im = PIL.Image.open('/seq/vgb/swofford/temp/clumpiness/feature_' + p + '.png')
        imgs.append(im)
                
        otherFeat = normCount[-normCount['label'].isin([i,p])][topShap]
        otherFeat = otherFeat.unstack().reset_index()
        otherFeat.columns = ['feature', 'sample', 'value']
        plt.clf()
        sns.violinplot(data = otherFeat, x = 'value', y = 'feature', cut = 0, scale = 'count')
        plt.xlim(0,1)
        plt.tight_layout()
        plt.savefig('/seq/vgb/swofford/temp/clumpiness/feature_other_' + i + '_' + p + '.png')
        im = PIL.Image.open('/seq/vgb/swofford/temp/clumpiness/feature_other_' + i + '_' + p + '.png')
        imgs.append(im)
                
        
        
    grid = image_grid(imgs, 6, 3)
    grid.save('/seq/vgb/swofford/temp/clumpiness/' + i + '_all.png', format='png')
    
    
newCol = []
for i in df2.columns:
    if i[:3] == 'ENS':
        fields = i.strip().split('_')
        name = '@'.join(list(reference[fields[0]].humanGeneName))
        outFields = [name]
        for i in fields:
            outFields.append(i)
        outName = '_'.join(outFields)
        newCol.append(outName)
    else:
        newCol.append(i)    
