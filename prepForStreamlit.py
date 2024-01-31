import pandas as pd
import shap
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import PIL



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

'''for shapley stuff website'''

'''feature Values for everyone'''
df = pd.read_csv('/seq/vgb/swofford/R01/convertingToSKLInput/testVCFs/artifactSigsRemovedNewlineFixed_080522/artifactSigsRemovedNewlineFixed_080522_featureValues.txt', sep = '\t')
labels = pd.unique(df['label'])

for cancer in labels:
    shapdf = pd.read_csv('/seq/vgb/swofford/R01/convertingToSKLInput/testVCFs/artifactSigsRemovedNewlineFixed_080522/shapley/' + cancer + '_averageShapley.txt', sep = '\t')
    shapdf = shapdf[shapdf['label'] == cancer].copy(deep = True)
    shapdf.sort_values('ID', inplace = True)
    shapdf.drop(['label', 'ID','pred','prob_' + cancer], inplace = True, axis = 1)
    featdf = df[df['label'] == cancer].copy(deep = True)
    featdf.sort_values('ID', inplace = True)
    featdf.drop(['ID', 'label', 'species', 'factorLabel', 'is_train'], inplace = True, axis = 1)
    
    pos5 = shapdf.reindex(shapdf.max(numeric_only = True).sort_values(ascending = False).index, axis = 1).columns[:5]
    neg5 = shapdf.reindex(shapdf.min(numeric_only = True).sort_values(ascending = True).index, axis = 1).columns[:5]
    
    topTen = set()
    
    for i in pos5:
        topTen.add(i)
    
    for i in neg5:
        topTen.add(i)
    topTen = list(topTen)   
    plt.clf()
    shap.summary_plot(shapdf[topTen].to_numpy(), features = featdf[topTen].to_numpy(), feature_names = topTen, show=True, sort = False, plot_type = 'violin')
    #shap.summary_plot(shapdf.to_numpy(), features = featdf.to_numpy(), feature_names = featdf.columns, show=False, sort = True, plot_type = 'violin', max_display = 10)
    plt.savefig('/seq/vgb/swofford/streamlit/data/shapley/' + cancer + '.jpg', bbox_inches='tight', dpi=300)


'''gaddy grams sigs'''
df = pd.read_csv('/seq/vgb/swofford/R01/convertingToSKLInput/testVCFs/allSamplesProcessedArtifactSigsRemovedLowPerformingClassesRemoved_080522.txt', sep = '\t')
sigs = []
updateCols = ['label', 'species']
for i in df.columns:
    if 'SBS' in i:
        sigs.append(i)
        updateCols.append(i)
df = df[updateCols].copy(deep = True)

for sig in sigs:
    df['logPlotCol'] = np.log(df[sig])
    dfs = [x for _, x in df.groupby('label')]
    meanValsDict = {}
    meanVals = []
    maxN = 0
    for i in dfs:
        if len(i) > maxN:
            maxN = len(i)
    for i in range(len(dfs)):
        val = dfs[i]['logPlotCol'].mean()
        sampleNum = len(dfs[i])
        factor = maxN/sampleNum
        meanVals.append(val)
        meanValsDict[val] = dfs[i]
        dfs[i].sort_values('logPlotCol', ignore_index = True, inplace = True)
        dfs[i]['plotIndex'] = dfs[i].index
        dfs[i]['plotIndex'] = factor * dfs[i]['plotIndex']
    meanVals.sort()
    toMerge = []
    for i in meanVals:
        toMerge.append(meanValsDict[i])
        
    df2 = pd.concat(toMerge, ignore_index = True)
    
    plt.clf()
    g = sns.relplot(data=df2, kind = 'scatter', x = 'plotIndex', y = 'logPlotCol', col = 'label', hue = 'species', height=6, aspect=.1, linewidth=0, alpha = 0.7, facet_kws={'sharey': True, 'sharex': True})
    g.despine(left=True, bottom=True) 
    #for m, ax in zip(df2.groupby('label', sort=False)['logPlotCol'].mean(), g.axes.ravel()):
    #    ax.hlines(m,0,maxN, colors = 'gray', ls='--')
    for label, ax in g.axes_dict.items():
        ax.set_title(label, rotation = 45)
        ax.tick_params(bottom=False, left=False)
        ax.hlines(df2[df2['label'] == label]['logPlotCol'].mean(),0,maxN, colors = 'gray', ls='--')
    g.set(xticklabels=[])
    g.set(xlabel=None)
    g.set(ylabel='log ' + sig)
    plt.savefig('/seq/vgb/swofford/streamlit/data/sigPics/' + sig + '.png', bbox_inches='tight', dpi = 300)
 
 
'''
sig stuff for Diane
'''
shaps = {}
for cancer in labels:
    shapdf = pd.read_csv('/seq/vgb/swofford/R01/convertingToSKLInput/testVCFs/artifactSigsRemovedNewlineFixed_080522/shapley/' + cancer + '_averageShapley.txt', sep = '\t')
    shapdf = shapdf[shapdf['label'] == cancer].copy(deep = True)
    shapdf = shapdf[shapdf['pred'] == cancer].copy(deep = True)
    shapdf.sort_values('ID', inplace = True)
    shapdf.drop(['label', 'ID','pred','prob_' + cancer], inplace = True, axis = 1)
    featdf = df[df['label'] == cancer].copy(deep = True)
    featdf.sort_values('ID', inplace = True)
    featdf.drop(['ID', 'label', 'species', 'factorLabel', 'is_train'], inplace = True, axis = 1)
    
    pos10 = shapdf.reindex(shapdf.max(numeric_only = True).sort_values(ascending = False).index, axis = 1).columns[:10]
    neg10 = shapdf.reindex(shapdf.min(numeric_only = True).sort_values(ascending = True).index, axis = 1).columns[:10]
    
    topTwenty = set()
    
    for i in pos10:
        topTwenty.add(i)
    
    for i in neg10:
        topTwenty.add(i)
    topTwenty = list(topTwenty)
    shapdf = pd.read_csv('/seq/vgb/swofford/R01/convertingToSKLInput/testVCFs/artifactSigsRemovedNewlineFixed_080522/shapley/' + cancer + '_averageShapley.txt', sep = '\t')
    shapdf = shapdf[shapdf['label'] == cancer].copy(deep = True)
    shapdf = shapdf[shapdf['pred'] == cancer].copy(deep = True)
    keep = ['ID','label','pred','prob_' + cancer]
    for i in topTwenty:
        keep.append(i)
    shapdf = shapdf[keep].copy(deep = True)
    shapdf.columns = shapdf.columns.str.replace('prob_' + cancer, 'prob')   
    shaps[cancer] = shapdf

dfs = []
for i in shaps.keys():
    dfs.append(shaps[i])

test = pd.concat(dfs)
test = test.groupby(['label','pred']).mean()
test.columns = test.columns.str.replace('prob','avgProbOfAccuratePred')
test.to_csv('topTwenyFeaturesPredsEqualLabsMeanShapValues.txt', sep = '\t')

keep = set()
for i in shaps.keys():
    for x in shaps[i].columns:
        if 'prob' not in x:
            keep.add(x)
            
            
'''
    Trying to get plot with:
    class(topFeature) classTopFeatureShapViolin allOtherClassTopFeatureShapViolin
'''
def getGeneName(ens):
    return '@'.join(list(reference[ens].humanGeneName))


df = pd.read_csv('/seq/vgb/swofford/R01/convertingToSKLInput/testVCFs/artifactSigsRemovedNewlineFixed_080522/artifactSigsRemovedNewlineFixed_080522_featureValues.txt', sep = '\t')
labels = pd.unique(df['label'])
reference = buildRef('/seq/vgb/swofford/dogHumanOrthologs/old/feb2/ncbiGeneNameAdded_zeroRemoved_030623.bed', '/seq/vgb/swofford/ref/pathways/C2CuratedGeneSetsGeneNames.gmt.txt')
df = df[df['species'] != 'dog']

plt.figure().clear()
plt.close()
plt.cla()
plt.clf()

#cancer = 'AS'
#fig, (ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8, ax9, ax10, ax11, ax12, ax13, ax14, ax16, ax17) = plt.subplots(nrows = 16, ncols = 1)
#for i in range(8):
#cancer = labels[i]
for cancer in labels:
    if 'dog' not in cancer:
        fig, (ax1, ax2) = plt.subplots(nrows = 2, ncols = 1)
        shapdf = pd.read_csv('/seq/vgb/swofford/R01/convertingToSKLInput/testVCFs/artifactSigsRemovedNewlineFixed_080522/shapley/' + cancer + '_averageShapley.txt', sep = '\t')
        shapdf = shapdf[~shapdf['label'].str.contains('dog')].copy(deep = True)
        shapdf.loc[shapdf.label != cancer, 'label'] = 'other'
        #get top shapley value
        tempShap = shapdf[shapdf['label'] == cancer].copy(deep = True)
        tempShap = tempShap[tempShap['pred'] == cancer]
        tempShap.sort_values('ID', inplace = True)
        keepIDs = tempShap['ID'].values
        tempShap.drop(['label', 'ID','pred','prob_' + cancer], inplace = True, axis = 1)
        
        top = tempShap.reindex(tempShap.abs().max(numeric_only = True).sort_values(ascending = False).index, axis = 1).columns[0]
        
        other = shapdf[shapdf.label == 'other'].copy(deep = True)
        other.sort_values('ID', inplace = True)
        other.drop(['label','ID','pred','prob_' + cancer], axis = 1, inplace = True)
        
        featdf = df.copy(deep = True)
        featdf.loc[featdf.label != cancer, 'label'] = 'other'
        featdf[['label', top]]
        otherFeat = featdf[featdf['label'] == 'other'].copy(deep = True)
        otherFeat.sort_values('ID', inplace = True)
        otherFeat.drop(['ID', 'label', 'species', 'factorLabel', 'is_train'], inplace = True, axis = 1)
        featdf = featdf.loc[featdf['ID'].isin(keepIDs)]
        featdf.sort_values('ID', inplace = True)
        featdf.drop(['ID', 'label', 'species', 'factorLabel', 'is_train'], inplace = True, axis = 1)
        plt.subplot(2,1,1)
        shap.summary_plot(tempShap[[top]].to_numpy(), features = featdf[[top]].to_numpy(), feature_names = [''], show=False, plot_size = None)
        plt.subplot(2,1,2)
        shap.summary_plot(other[[top]].to_numpy(), features = otherFeat[[top]].to_numpy(), feature_names = [''], show=False, plot_size = None, color_bar_label = None)
        ax1.set_title(cancer, fontdict = {'fontsize':'xx-large', 'fontweight':'extra bold'})
        ax2.set_title('other', fontdict = {'fontsize':'xx-large', 'fontweight':'extra bold'})
        ax1.set_xlabel('')
        if 'ENS' in top:
            if 'count' in top:
                name = getGeneName(top.split('_')[0]) + '_count'
            else:
                name = getGeneName(top)
        else:
            name = top
        ax2.set_xlabel(name, fontdict = {'fontsize':'xx-large', 'fontweight':'extra bold'})    #ax3.set_xlabel('')
        plt.subplots_adjust(wspace = 6.0)
        plt.rcParams["axes.edgecolor"] = "black"
        plt.rcParams["axes.linewidth"] = 1
        plt.tight_layout()
        
        
        #plt.subplots_adjust(wspace=6.0)
        #plt.tight_layout()
        
        
        plt.savefig('/seq/vgb/swofford/temp/topShapImages/' + cancer + '.png')   
    
    
def image_grid(imgs, rows, cols):
    assert len(imgs) == rows*cols
    w, h = imgs[0].size
    grid = Image.new('RGB', size=(cols*w, rows*h))
    grid_w, grid_h = grid.size    
    for i, img in enumerate(imgs):
        grid.paste(img, box=(i%cols*w, i//cols*h))
    return grid

cancerImages = []
for cancer in labels:
    if 'dog' not in cancer:
        im = PIL.Image.open('/seq/vgb/swofford/temp/topShapImages/' + cancer + '.png')
        im = PIL.ImageOps.expand(im, border = 5, fill = 'black')
        cancerImages.append(im)
grid = image_grid(cancerImages, 8, 4)
grid.save('/seq/vgb/swofford/temp/topShapImages/all.png', format='png')
    
    pos5 = shapdf.reindex(shapdf.max(numeric_only = True).sort_values(ascending = False).index, axis = 1).columns[:5]
    neg5 = shapdf.reindex(shapdf.min(numeric_only = True).sort_values(ascending = True).index, axis = 1).columns[:5]
    
    topTen = set()
    
    for i in pos5:
        topTen.add(i)
    
    for i in neg5:
        topTen.add(i)
    topTen = list(topTen)   
    plt.clf()
    shap.summary_plot(shapdf[topTen].to_numpy(), features = featdf[topTen].to_numpy(), feature_names = topTen, show=True, sort = False, plot_type = 'violin')
    #shap.summary_plot(shapdf.to_numpy(), features = featdf.to_numpy(), feature_names = featdf.columns, show=False, sort = True, plot_type = 'violin', max_display = 10)
    plt.savefig('/seq/vgb/swofford/streamlit/data/shapley/' + cancer + '.jpg', bbox_inches='tight', dpi=300)
    
    

hg38CDSLabeled.bed
 