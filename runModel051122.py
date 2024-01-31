# source activate /seq/vgb/swofford/conda/my3env
# python runModel.py -i yourFakeDataFile -f features file (ex. importantFeatures_fold_0.pickle) -m modelFile (ex. xgb_model_fold_0.pickle) -o outFile
# will outPut a file with columns label, pred, ID, species
# the example_fold_*.txt files are files with all of the real samples and their feature values for each of the 5 fold's selected features.



import pandas as pd
import numpy as np
import xgboost as xgb
import argparse
import pickle
from sklearn.metrics import confusion_matrix, ConfusionMatrixDisplay
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--inFile', default = '')
#parser.add_argument('-m', '--modelFile')
#parser.add_argument('-f', '--featureFile')
parser.add_argument('-o', '--outPrefix')
args = parser.parse_args()

models = {}
features = {}
allFeatures = set()
preds = {}
probs = {}
dogPreds = {}
dogProbs = {}

outFeatures = ['ID','label','species','factorLabel','is_train']
              
for i in range(5):
    model = xgb.XGBClassifier()
    model.load_model('/seq/vgb/swofford/R01/convertingToSKLInput/testVCFs/finalOutput051122/xgb_model_fold_' + str(i) + '.json')
    models[i]= model
    with open('/seq/vgb/swofford/R01/convertingToSKLInput/testVCFs/finalOutput051122/importantFeatures_fold_' + str(i) + '.pickle', 'rb') as inFile:
        features[i] = pickle.load(inFile)
        for x in features[i]:
            allFeatures.add(x)

if args.inFile == '':
    df = pd.read_csv('/seq/vgb/swofford/R01/convertingToSKLInput/testVCFs/full_updatedLabels_preProcessed.txt', sep = '\t')
    outFeatures = np.append(outFeatures,list(allFeatures))
    df[outFeatures].to_csv(args.outPrefix + '_featureValues.txt', sep = '\t', index=False)
else:
    df = pd.read_csv(args.inFile, sep = '\t')
labels = np.unique(df['label'])
dogs = df[df['species'] == 'dog'].copy(deep = True)
dogs.set_index('ID', inplace = True)
df = df[df['species'] != 'dog'].copy(deep = True)
df.set_index('ID', inplace = True)
for i in range(5):
    test = df[df['is_train'] == i].copy(deep = True)
    test['pred'] = labels[models[i].predict(test[features[i]])]
    preds[i] = test[['label','pred']].copy(deep = True)
    dogs['pred'] = labels[models[i].predict(dogs[features[i]])]
    dogPreds[i] = dogs[['label','pred']].copy(deep = True)
    prob = pd.DataFrame(models[i].predict_proba(test[features[i]]))
    prob.columns = labels[prob.columns]
    columns = prob.columns
    prob['label'] = test['label'].values
    prob['pred'] = preds[i]['pred'].values
    prob['ID'] = test.index.values
    prob = prob[np.append(['ID','label','pred'],columns)]
    prob = prob.set_index('ID')
    probs[i] = prob.copy(deep = True)
    dogProb = pd.DataFrame(models[i].predict_proba(dogs[features[i]]))
    dogProb.columns = labels[dogProb.columns]
    columns = dogProb.columns
    dogProb['label'] = dogs['label'].values
    dogProb['pred'] = dogPreds[i]['pred'].values
    dogProb['ID'] = dogs.index.values
    dogProb = dogProb[np.append(['ID','label','pred'],columns)]
    dogProb = dogProb.set_index('ID')
    dogProbs[i] = dogProb.copy(deep = True)

allDogProbs = pd.concat((dogProbs[0],dogProbs[1],dogProbs[2],dogProbs[3],dogProbs[4]), axis=0)
#tempDogLabels = dogs['label'].copy(deep=True)
allDogProbs = allDogProbs.groupby(allDogProbs.index).mean()
allDogProbs['pred'] = allDogProbs.idxmax(axis = 1)
allDogProbs['label'] = dogs['label']
allProbs = pd.concat((probs[0],probs[1],probs[2],probs[3],probs[4], allDogProbs), axis=0)
allProbs.groupby(allProbs.index).mean()
allProbs.to_csv(args.outPrefix + '_probs.txt', sep = '\t')

cm = confusion_matrix(allProbs['label'], allProbs['pred'], labels=labels)
ConfusionMatrixDisplay.from_predictions(allProbs['label'],allProbs['pred'], 
                                        display_labels=labels, normalize = 'true', 
                                        include_values = False, xticks_rotation = 'vertical')
figure = plt.gcf()
figure.set_size_inches(14,14)
plt.savefig('confusionMatrixAll.pdf')
#meanDogProbs = pd.concat((dogProbs[0],dogProbs[1],dogProbs[2],dogProbs[3],dogProbs[4]), axis=1).mean(axis=1)
#humProbs = pd.concat((probs[0],probs[1],probs[2],probs[3],probs[4]))





'''
ConfusionMatrixDisplay.from_predictions(test['label'],labels[preds], display_labels=labels, normalize = 'true', include_values = False, xticks_rotation = 'vertical')


with open(args.modelFile, 'rb') as modelIn, open(args.featureFile, 'rb') as featureIn, open('/seq/vgb/swofford/R01/convertingToSKLInput/testVCFs/labels.pickle', 'rb') as labelsIn:
    xgb_model = pickle.load(modelIn)
    features = pickle.load(featureIn)
    labels = pickle.load(labelsIn)   

if args.inFile == '':
    df = pd.read_csv('/seq/vgb/swofford/R01/convertingToSKLInput/testVCFs/full_updatedLabels_preProcessed.txt', sep = '\t')
    outFeatures = np.append(outFeatures,features)
    df[outFeatures].to_csv(args.outFile, sep = '\t', index=False)
    

else:
    df = pd.read_csv(args.inFile, sep = '\t')
    preds = xgb_model.predict(df[features])
    probs = pd.DataFrame(xgb_model.predict_proba(df[features]))
    probs.columns = labels[probs.columns]
    outCols = np.append(['label','pred','ID','species'], probs.columns)
    probs['label']  = df['label'].values
    probs['pred'] = labels[preds]
    probs['ID'] = df['ID'].values
    probs['species'] = df['species'].values
    #outPut = pd.DataFrame({'label': df['label'], 'pred': labels[preds], 'ID': df['ID'], 'species':df['species']})
    #outPut.to_csv(args.outFile, sep = '\t', index=False)
    #probs[outCols].to_csv(args.outFile, sep = '\t', index = False)

 below is just keeping track of manually run stuff in case we want to do it again
triFeatures = {}
noTriFeatures = {}
triModels = {}
noTriModels = {}
for i in range(5):
    inFile = open('importantFeatures_fold_' + str(i) + '.pickle','rb')
    triFeatures[i] = pickle.load(inFile)
    inFile.close()
    inFile = open('./old/importantFeatures_fold_' + str(i) + '.pickle','rb')
    noTriFeatures[i] = pickle.load(inFile)
    inFile.close()
    inFile = open('xgb_model_fold_' + str(i) + '.pickle', 'rb')
    triModels[i] = pickle.load(inFile)
    inFile.close()
    inFile = open('./old/xgb_model_fold_' + str(i) + '.pickle', 'rb')
    noTriModels[i] = pickle.load(inFile)
    inFile.close()

dogs = df[df['species'] == 'dog']
for i in range(5):
    test = df[df['is_train'] == i]
    test = test[test['species'] != 'dog']
    test = test.append(dogs)
    model = noTriModels[i]
    features = noTriFeatures[i]
    preds = labels[model.predict(df[features])]
    probs = pd.DataFrame(model.predict_proba(df[features]))
    probs.columns = labels[probs.columns]
    outCols = np.append(['label','pred','ID','species','is_train'], probs.columns)
    probs['label']  = df['label'].values
    probs['pred'] = preds
    probs['ID'] = df['ID'].values
    probs['species'] = df['species'].values
    probs['is_train'] = df['is_train'] != i
    probs.loc[probs.species == 'dog', 'is_train'] = False
    probs[outCols].to_csv('noTriPredictions_trainIncluded_' + str(i) + '.txt', sep = '\t', index = False)
    model = triModels[i]
    features = triFeatures[i]
    preds = labels[model.predict(df[features])]
    probs = pd.DataFrame(model.predict_proba(df[features]))
    probs.columns = labels[probs.columns]
    outCols = np.append(['label','pred','ID','species', 'is_train'], probs.columns)
    probs['label']  = df['label'].values
    probs['pred'] = preds
    probs['ID'] = df['ID'].values
    probs['species'] = df['species'].values
    probs['is_train'] = df['is_train'] != i
    probs.loc[probs.species == 'dog', 'is_train'] = False
    probs[outCols].to_csv('triPredictions_trainIncluded_' + str(i) + '.txt', sep = '\t', index = False)
'''
    
    
    
     

    
            