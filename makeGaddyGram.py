import seaborn as sns
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-c', '--categoryColumn', default = 'label')
parser.add_argument('-i', '--inFile', default = '/seq/vgb/swofford/R01/convertingToSKLInput/testVCFs/inProcess071923_addingAllFeaturesToSingleFile_nonDamageAdded_hotSpotAdded.txt.gz')
parser.add_argument('-p', '--plotColumn', default = 'mutCount')
parser.add_argument('-y', '--ylabel', default = 'log n Mutations')
parser.add_argument('-o', '--outFile', default = 'gaddyTest.png')
parser.add_argument('-hu', '--hueCol', default = 'species')
args = parser.parse_args()

df = pd.read_csv(args.inFile, sep = '\t', usecols = [args.categoryColumn, args.plotColumn, args.hueCol])
df['logPlotCol'] = np.log(df[args.plotColumn])
#taking the log of any zero values results in - infinity, so replace those with zeros
df.replace(-np.inf, 0,inplace = True)
dfs = [x for _, x in df.groupby(args.categoryColumn)]
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
if args.hueCol == 'species':
    g = sns.relplot(data=df2, kind = 'scatter', x = 'plotIndex', 
                    y = 'logPlotCol', col = args.categoryColumn, 
                    hue = args.hueCol, height=6, aspect=.1, 
                    linewidth=0, alpha = 0.7, palette = {'dog':'#7570B3', 'Human':'#1B9E77'}, 
                    facet_kws={'sharey': True, 'sharex': True})

else:   
    g = sns.relplot(data=df2, kind = 'scatter', x = 'plotIndex', 
                    y = 'logPlotCol', col = args.categoryColumn, 
                    hue = args.hueCol, height=6, aspect=.1, linewidth=0, 
                    alpha = 0.7, facet_kws={'sharey': True, 'sharex': True})
g.despine(left=True, bottom=True) 

for label, ax in g.axes_dict.items():
    ax.set_title(label, rotation = 45)
    ax.tick_params(bottom=False, left=False)
    ax.hlines(df2[df2[args.categoryColumn] == label]['logPlotCol'].mean(),0,maxN, colors = 'gray', ls='--')
g.set(xticklabels=[])
g.set(xlabel=None)
g.set(ylabel=args.ylabel)
plt.savefig(args.outFile, bbox_inches='tight')


''' 
df['denPlotLab'] = df['species']
df.loc[df['label'] == 'MELA', 'denPlotLab'] = 'MELA'
UVSigs = ['SBS7a', 'SBS7b', 'SBS7c', 'SBS7d']
for i in UVSigs:
    plt.clf()
    sns.kdeplot(data=df[df['species'] == 'dog'], x=i, hue='label', fill = True, common_norm = True, alpha=0.4)
    plt.savefig(i + '_Dogsdensity.pdf', bbox_inches='tight')
    
for i in UVSigs:
    plt.clf()
    sns.kdeplot(data=df, x=i, hue='denPlotLab', fill = True, common_norm = False, alpha=0.4)
    plt.savefig(i + '_density.pdf', bbox_inches='tight')
'''
