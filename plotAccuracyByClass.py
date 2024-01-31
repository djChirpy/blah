import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scipy as sp

plt.figure()
plt.figure().clear()
plt.close()
plt.cla()
plt.clf()


sns.set(rc={'figure.figsize':(20,10)})

probs = pd.read_csv('/seq/vgb/cancer_r01/allxallCode/test092823Orig/test092823Orig_probs.txt', sep = '\t')

correct = {}

labels = pd.unique(probs.label)

for i in labels:
    if 'dog' not in i:
        temp = probs[probs['label'] == i]
        correct[i] = (len(temp),len(temp[temp['pred'] == i])/len(temp['label']))
        
 

temp = pd.DataFrame.from_dict(correct, orient = 'index')

temp.columns = ['count','fractionCorrect']

ax = sns.regplot(data = temp, x = 'count', y = 'fractionCorrect')

r, p = sp.stats.pearsonr(temp['count'], temp['fractionCorrect'])

for i in temp.index:
    if i == 'HGG' or i == 'PACA':
        plt.text(x = temp.loc[i]['count'] + 10, y = temp.loc[i]['fractionCorrect'] +.01, s = i, size='small')
    else:
        plt.text(x = temp.loc[i]['count'] + 10, y = temp.loc[i]['fractionCorrect'] -.01, s = i, size='small')

plt.text(1300, .3, 'r={:.2f}, p={:.2g}'.format(r, p))

plt.ylim(0,1)    
plt.ylabel('Prediction Accuracy')
plt.xlabel('Number of Samples in Class')

plt.tight_layout()

plt.savefig('accuracyBySampleNumber.pdf')


