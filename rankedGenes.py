import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt



df = pd.read_csv('driverFrequenciesAllClasses.txt', sep = '\t')
rankedGenes = {}
labels = pd.unique(df['label'])
features = df.columns[1:]


for i in labels:
    rankedGenes[i] = df[df['label'] == i][features].rank(axis=1, ascending = False)


for i in rankedGenes:
    rankedGenes[i]['label'] = i
    rankedGenes[i] = rankedGenes[i].set_index('label')

keep = pd.DataFrame()
topTenSomewhere = set()

for i in labels:
    test = df[df['label'] == i][features].rank(axis=1, ascending = False)
    for x in test.columns:
        if test[x].values <= 10:
            topTenSomewhere.add(x)

for i in rankedGenes:
    keep = pd.concat([keep, rankedGenes[i][list(topTenSomewhere)]])

keep = keep.reset_index()

for i in topTenSomewhere:
    keep.loc[keep[i] > 10, i] = 0


keep = keep.set_index('label')
keep = keep.reindex(keep.mean().sort_values(ascending = False).index, axis = 1)
    
#for i in topTenSomewhere:
#    keep.loc[keep[i] > 0, i] = abs(11-keep[i])
    
uns = keep.unstack().reset_index()
uns.columns = ['gene', 'class', 'rank']
uns['size'] = abs(11 - uns['rank'])
uns.loc[uns['size'] == 11, 'size'] = 0
uns.loc[uns['size'] != 0, 'size'] = uns['size'] * 10


plt.clf()
sns.scatterplot(data = uns, x = 'gene', y = 'class', hue = 'rank', size = 'size')
plt.xticks(rotation=90)
plt.gcf().set_size_inches(30,15)
plt.grid(visible = True, which = 'both', axis = 'x', color='gray', linestyle = '-', linewidth=0.075)
plt.savefig('topTenRanked.pdf', bbox_inches = 'tight')
plt.close()



    
    