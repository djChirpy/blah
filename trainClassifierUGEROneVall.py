import pandas as pd
import numpy as np
import xgboost as xgb
import matplotlib.pyplot as plt
from sklearn.metrics import plot_confusion_matrix
from sklearn.utils import class_weight
import matplotlib.pyplot as plt
import pickle
from sklearn.metrics import ConfusionMatrixDisplay, confusion_matrix, accuracy_score
import argparse
import os




parser = argparse.ArgumentParser()
parser.add_argument('-tr', '--trainFile', nargs=1)
parser.add_argument('-te', '--testFile', nargs=1)
parser.add_argument('-df', '--dogFile', nargs=1)
parser.add_argument('-o', '--outPath', nargs=1)
parser.add_argument('-n', '--nThreads', nargs=1, default = 1)
parser.add_argument('-e', '--nEstimators', nargs=1, default = 500)
parser.add_argument('-d', '--maxDepth', nargs=1, default = 10)
parser.add_argument('-s', '--stoppingRounds', nargs=1, default = 10)
parser.add_argument('-sub', '--subSample', nargs=1, default = 0)

args = parser.parse_args()

if not os.path.exists(args.outPath):
    os.mkdir(args.outPath)
    

#xgboost in native form stuff
train = pd.read_csv(args.trainFile, sep = '\t')
test = pd.read_csv(args.testFile, sep = '\t')
dog = pd.read_csv(args.dogFile, sep = '\t')

#class_weights = list(class_weight.compute_class_weight('balanced', np.unique(train['label']),train['label']))
y_train = pd.factorize(train['label'])[0]
y_test = pd.factorize(test['label'])[0]

'''
w_array = np.ones(y_train.shape[0], dtype = 'float')
for i, val in enumerate(y_train):
    w_array[i] = class_weights[val]
w_array = 1/w_array
    
train.columns = train.columns.str.replace("[\[,\]<]","")
features = train.columns[3:]
test.columns = test.columns.str.replace("[\[,\]<]","")


xgb_model = xgb.XGBClassifier(objective='multi:softprob', random_state=42, 
                              subsample = args.subSample, nthread = args.nThreads, max_depth = args.maxDepth, 
                              n_estimators = args.nEstimators)

eval_set = [(train[features], y_train), (test[features], y_test)]

xgb_model.fit(train[features], y_train, early_stopping_rounds=10,
    eval_metric=["merror", "mlogloss"], eval_set=eval_set, verbose=True, sample_weight=w_array)

#xgb_model = xgb.XGBClassifier(objective='multi:softprob', random_state=42, verbosity = 2, subsample = .8)
#xgb_model.fit(df[features], y_train, sample_weight=w_array)
#test = pd.read_csv('updatedOrthologsNoCountsCumulativeForLearning.out_humanTest.out', sep = '\t')
#test.columns = test.columns.str.replace("[\[,\]<]","")
#y_test = pd.factorize(test['label'])[0]
#y_pred = xgb_model.predict(test[features])

plot_confusion_matrix(xgb_model, test[features], y_test, normalize = 'true', 
                      xticks_rotation = 'vertical', display_labels = labels, include_values = False)
labels = np.unique(train['label'])
figure = plt.gcf()
figure.set_size_inches(14,14)
plt.savefig('xgbConfusion012821_1_true.pdf')
#48% accuracy with default values (above)
xgb_model = xgb.XGBClassifier(objective='multi:softprob', random_state=42, verbosity = 2, subsample = .8, nthread = 4, max_depth = 10, n_estimators = 500)

eval_set = [(df[features], y_train), (test[features], y_test)]
xgb_model.fit(df[features], y_train, early_stopping_rounds=10,
    eval_metric=["error", "mlogloss"], eval_set=eval_set, verbose=True, sample_weight=w_array)
'''
'''OneVall don all at once for speeds'''
#featureFile = open('1611FeaturesSharedWithDog.pickle', 'rb')
#features = pickle.load(featureFile)
#featureFile.close()
labelTest = np.unique(train['label'])
label = train['label']
train.columns = train.columns.str.replace("[\[,\]<]","")
features = train.columns[3:]
test.columns = test.columns.str.replace("[\[,\]<]","")
dog.columns = dog.columns.str.replace("[\[,\]<]","")
train['category'] = 'train'
test['category'] = 'test'
dog['category'] = 'dog'
from sklearn.metrics import accuracy_score
for j in labelTest:
    oneVall = []
    for lab in train['label']:
        if lab != j:
            oneVall.append(lab)
    
            
    useLabel = train['label'].replace(to_replace = oneVall, value = 'other')
    testLabel = test['label'].replace(to_replace = oneVall, value = 'other')
    y, lookup = pd.factorize(useLabel)
    useTestLabel = test['label'].replace(to_replace = oneVall, value = 'other')
    #class_weights = list(class_weight.compute_class_weight('balanced', np.unique(useLabel),useLabel))
    #print(class_weights)
    y_train = pd.factorize(useLabel)[0]
    #print("y_train = " + str(y_train))
    y_test = pd.factorize(useTestLabel)[0]
    #w_array = np.ones(y_train.shape[0], dtype = 'float')
    #for i, val in enumerate(y_train):
    #    w_array[i] = class_weights[val]
    #w_array = 1/w_array
    eval_set = [(train[features], y_train), (test[features], y_test)]
    #xgb_model = xgb.XGBClassifier(objective='multi:softprob', random_state=42, subsample = .8, nthread = 5, max_depth = 10, n_estimators = 1000, learning_rate=0.05)
    xgb_model = xgb.XGBClassifier(objective='binary:logistic', random_state=42, subsample = .8, nthread = 10, max_depth = 10, n_estimators = 1000, learning_rate=0.05)
    #xgb_model.fit(train[features], y_train, early_stopping_rounds=10, eval_metric=["error", "logloss"], eval_set=eval_set, verbose=True, sample_weight=w_array)
    xgb_model.fit(train[features], y_train, early_stopping_rounds=10, eval_metric=["error", "logloss"], eval_set=eval_set, verbose=True, sample_weight=class_weight.compute_sample_weight("balanced", y_train))
    featureImportances = xgb_model.feature_importances_
    featureImportances.sort()
    cutoff = 0
    i = -200
    while cutoff == 0:
        cutoff = min(featureImportances[i:])
        i += 1
    keep = np.array(['label', 'category','ID','species'])
    topFeature = []
    featureVal = []
    for i, val in enumerate(xgb_model.feature_importances_):
        if val >= cutoff:
            topFeature.append(features[i])
            featureVal.append(val)
    featureImportance = pd.DataFrame({'feature':topFeature, 'importance':featureVal})
    featureImportance.sort_values(by = 'importance', ascending = False, inplace = True)
    featureImportance.to_csv('oneVall/featureImportances/' + str(j) + '_top200Features.txt', sep = '\t', index=False)
    keep = np.append(keep, featureImportance['feature'].to_numpy())
    fullPredict = lookup[xgb_model.predict(train[features].append(test[features], ignore_index = True).append(dog[features], ignore_index=True))]
    fullSet = pd.DataFrame({'prediction':fullPredict})
    fullSet = fullSet.join(train[keep].append(test[keep], ignore_index=True).append(dog[keep], ignore_index=True))
    testLabel = fullSet['label'].replace(to_replace = oneVall, value = 'other')
    fullSet['correct'] = testLabel == fullSet['prediction']
    keep = np.append(['label', 'category','ID','species', 'correct', 'prediction'], featureImportance['feature'].to_numpy())
    fullSet[keep].to_csv('oneVall/sampleFeatureValues/' + str(j) + '_allSamplesTop200Features.txt', sep = '\t', index = False)
    preds = xgb_model.predict(test[features])
    accuracy = accuracy_score(y_test, preds)
    print(str(accuracy))
    comboTest = test.append(dog, ignore_index=True)    
    fullPreds = xgb_model.predict(comboTest[features])
    outFile = open('oneVall/models/' + str(j) + '_' + str(accuracy) + '.pickle', 'wb')
    pickle.dump(xgb_model, outFile)
    outFile.close()
    fullPred_label = lookup[fullPreds]
    preds = fullPred_label.to_numpy()
    actual = comboTest['label'].to_numpy()
    y, matrixLabels = pd.factorize(comboTest['label'])
    matrixLabels = np.append(np.array(matrixLabels),'other')
    labeledMatrix = confusion_matrix(actual,preds,labels=matrixLabels)
    disp = ConfusionMatrixDisplay(confusion_matrix=labeledMatrix, display_labels=matrixLabels)
    disp.plot(xticks_rotation='vertical')
    figure = plt.gcf()
    figure.set_size_inches(21,21)
    plt.savefig('oneVall/confusionMatrix/' + str(j) + '.pdf')

'''
from sklearn.ensemble import RandomForestClassifier
#from sklearn.model_selection import cross_val_score
#from sklearn import preprocessing
from sklearn.feature_selection import SelectFromModel
import pandas as pd
import numpy as np
import argparse
import h2o
from h2o.automl import H2OAutoML
import os

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--inFile', nargs=1)
parser.add_argument('-o', '--outPrefix', nargs=1)
parser.add_argument('-v', '--oneVall', nargs=1, default='')
parser.add_argument('-l', '--limit', nargs=1, default=0)
args = parser.parse_args()

#read in file from parsing script
limit = args.limit

if not os.path.exists(args.inFile[0] + '_human.out'):
    df = pd.read_csv(args.inFile[0], sep = '\t')
    features = df.columns
    df['isDog'] = df.isin({'label':['dogHSA','dogOSA','dogMELANOMA','dogLSA']})['label']
    dogs = df[df['isDog'] == True]
    df = df[df['isDog'] == False]
    df.sort_values(by = 'label', inplace=True)
    dogs.sort_values(by = 'label', inplace=True)
    df['is_train'] = np.random.uniform(0, 1, len(df)) <= .80
    train, test = df[df['is_train']==True], df[df['is_train']==False]
    dogs[features].to_csv(args.inFile[0] + '_dog.out', sep = '\t', index = False)
    train[features].to_csv(args.inFile[0] + '_humanTrain.out', sep = '\t', index = False)
    test[features].to_csv(args.inFile[0] + '_humanTest.out', sep = '\t', index = False)
    df[features].to_csv(args.inFile[0] + '_human.out', sep = '\t', index = False)
    del dogs
    del train
    del test
    
else:
    df = pd.read_csv(args.inFile[0] + '_human.out', sep = '\t')
    #features = df.columns
 --working on subsetting/duplicating samples for ml input
if limit != 0:
    fractionToKeep = {}    
    labels = df.label.unique()
    for label in labels:
            fractionToKeep[label] = int(limit)/df.label.value_counts()[label]
    
    subsetting = {}
    for label in labels:
        df['temp'] = df['label'] == label
        subsetting[label] = df[df['temp'] == True]
        
        if len(subsetting[label]) >0.5*limit && len(subsetting[label]) < limit:
            #randomly select samples to duplicate
            subsetting[label]['temp'] = np.random.uniform(0, 1, len(subsetting[label]) )
        subsetting[label]['temp'] = np.random.uniform(0, 1, len(df)) <= .80
    
            
    
outPrefix = args.outPrefix[0]
modelDir = '/seq/vgb/swofford/R01/convertingToSKLInput/testVCFs/H2OModelsUGER/' + outPrefix

if not os.path.exists(modelDir):
    os.mkdir(modelDir)
    os.mkdir(modelDir + '/H2OModels1stLevel')
    os.mkdir(modelDir + '/H2OModels2ndLevel')

#create list of features
#features = df.columns[1:]

#set random seed for reproducibility testing
#np.random.seed(0)

#df['isDog'] = df.isin({'label':['dogHSA','dogOSA','dogMELANOMA','dogLSA']})['label']
#df['keep'] = df.isin({'label':['LSA','AS','MELA','OSA']})['label']
#separate dog (not used for training) data from human
#dogs = df[df['isDog'] == True]
#df = df[df['isDog'] == False]


#sort for labeling 

features = df.columns[1:]
useLabel = []

if args.oneVall == '':
    df['useLabel'] = df['label']
else:
    for entry in df['label']:
        if entry == args.oneVall[0]:
            useLabel.append(entry)
        else:
            useLabel.append('other')

df['useLabel'] = useLabel

#add column of factorized labels for training
df['factorLabel'] = pd.factorize(df['useLabel'])[0]

#create list of names to recover labels from factor labels
#names = pd.unique(df['label'])
        
# Train random forest on full data set to use as a feature selector
clf = RandomForestClassifier(n_jobs=-1, n_estimators = 2000, bootstrap=False, verbose = 2, class_weight = 'balanced')
clf.fit(df[features], df['factorLabel'])

#find features with greater than threshold importance
#0.00004 leaves about 6-7k features, I'm not sure yet
#about how few features best categorize classes
sfm = SelectFromModel(clf, threshold=0.00004)
#sfm = SelectFromModel(clf, threshold=0.00004)
sfm.fit(df[features], df['factorLabel'])

important_features = []
for feature_list_index in sfm.get_support(indices=True):
    important_features.append(features[feature_list_index])
    
with open(modelDir + '/importantFeatures.out', 'w') as featureOut:
    for feature in important_features:
        featureOut.write(feature + '\n')
    
#sfm2 = SelectFromModel(clf, threshold=0.0001)
#sfm2.fit(df[features], df['factorLabel'])
#important_feature3 = []

#for feature_list_index in sfm3.get_support(indices=True):
#    important_features3.append(features[feature_list_index])
    
#print('Features: ' + str(len(features)))
#print('Important features: ' + str(len(important_features)))
#print('Extra important: ' + str(len(important_features2)))

#trying to free up memory by clearing
#no-longer necessary     
del clf
del sfm

important_features.append('label')

#save separated dog df with important features
#dogs[important_features].to_csv(modelDir + '/tempDogs.temp', sep = '\t', index = False)
#dogLabels = dogs['label'].to_list()

#del dogs

#add factor Label to list of features to keep 

#df = df[important_features]

important_features.append('useLabel')

#df[important_features].to_csv(modelDir + '/tempDf.temp', sep = '\t', index = False)

del df

h2o.init()

train = h2o.import_file(args.inFile[0] + '_humanTrain.out')
test = h2o.import_file(args.inFile[0] + '_humanTest.out')
dogs = h2o.import_file(args.inFile[0] + '_dog.out')

#useLabel = []

if args.oneVall == '':
    train['useLabel'] = train['label']
    test['useLabel'] = test['label']
    
else:
    useLabel = []
    for entry in train['label'].as_data_frame()['label']:
        if entry == args.oneVall[0]:
            useLabel.append(entry)
        else:
            useLabel.append('other')
    train['useLabel'] = h2o.H2OFrame(useLabel)
    
    useLabel = []
    for entry in test['label'].as_data_frame()['label']:
        if entry == args.oneVall[0]:
            useLabel.append(entry)
        else:
            useLabel.append('other')
    test['useLabel'] = h2o.H2OFrame(useLabel)


    



#os.remove('tempDf.temp')
#os.remove('tempDogs.temp')

#train,test = df.split_frame(ratios=[.80,.20]) #70,20 best so far 60.7%
#train.as_data_frame().to_csv(modelDir + '/trainDf.out', sep = '\t', index = False)
#train2.as_data_frame().to_csv(modelDir + '/train2Df.out', sep = '\t', index = False)
#test.as_data_frame().to_csv(modelDir + '/testDf.out', sep = '\t', index = False)
x = train.columns
y = 'useLabel'
x.remove(y)
x.remove('label') 
train[y] = train[y].asfactor()
#train2[y] = train2[y].asfactor()
test[y] = test[y].asfactor()
aml = H2OAutoML(seed=1, max_runtime_secs=10000,
                balance_classes = True, 
                #exclude_algos = ['GLM', 'DRF', 'GBM'], 
                export_checkpoints_dir = modelDir + '/H2OModels1stLevel/')

#aml = H2OAutoML(seed=1, max_runtime_secs=20000,
#                balance_classes = True,  
#                export_checkpoints_dir = '/seq/vgb/swofford/R01/convertingToSKLInput/testVCFs/H2OModels2nd/')

aml.train(x=x, y=y, training_frame=train)

#Save classifier metrics to file
aml.leaderboard.as_data_frame().to_csv(modelDir + '/leaderBoard.out', sep = '\t', index = False)


#preds = aml.leader.predict(train[x])
preds = aml.leader.predict(test[x])


aml2 = H2OAutoML(seed=1, max_runtime_secs=13000,
                balance_classes = True,  
                export_checkpoints_dir = modelDir + '/H2OModels2ndLevel/')



labelPreds = test['useLabel'].cbind(preds)
#x2 = labelPreds.columns

#x2.remove('predict')
#x2.remove('label')
#x2.remove('useLabel')

labelPreds[y] = labelPreds[y].asfactor()
labelPreds['predict'] = labelPreds['predict'].asfactor()

aml2.train(x=x2, y=y, training_frame=labelPreds)

aml2.leaderboard.as_data_frame().to_csv(modelDir + '/leaderBoard2ndLevel.out', sep = '\t', index = False)


testPreds1st = aml.leader.predict(test[x])
testPreds1st = test['useLabel'].cbind(testPreds1st)
testPreds1st[y] = testPreds1st[y].asfactor()

testPreds2nd = aml2.leader.predict(testPreds1st[x2])
testPreds2nd = testPreds1st['useLabel'].cbind(testPreds2nd)


labelPreds.as_data_frame().to_csv(modelDir + '/testPredictions.out', sep = '\t', index = False)

dogPreds = aml.leader.predict(dogs[x])
dogPreds = dogs['label'].cbind(dogPreds)
#dogPreds[y] = dogPreds[y].asfactor()

#dogPreds2nd = aml2.leader.predict(dogPreds[x2])
#dogPreds2nd = dogPreds['label'].cbind(dogPreds2nd)

dogPreds.as_data_frame().to_csv(modelDir + '/dogPredictions.out', sep = '\t', index = False)


    
h2o.save_model(aml.leader, modelDir, force=True)
#h2o.save_model(aml2.leader, modelDir, force=True)
modelfile = aml.leader.download_mojo(path=modelDir, get_genmodel_jar=True)
#print("Model saved to " + modelfile)
#modelfile = aml2.leader.download_mojo(path=modelDir, get_genmodel_jar=True)
#print("2nd Layer Model saved to " + modelfile)


#best so far '/seq/vgb/swofford/R01/convertingToSKLInput/testVCFs/XGBoost_1_AutoML_20191112_234137'    



taking notes on a direct xgboost application

After loading the human data, the following splits into x (data) y (labels)
x,y = df.iloc[:,1:],df.iloc[:,0]
convert labesl to codes
ytest = y.astype('category').cat.codes

remove brackets and other characters from feature names, because it hates it
x.columns = x.columns.str.replace("[\[,\]<]","")
Convert to xgboost special format
data_dmatrix = xgb.DMatrix(data=x, label=ytest)


from sklearn.model_selection import train_test_split
x_train, x_test, y_train, y_test = train_test_split(x,y, test_size=0.2, random_state=123)
xg_class = xgb.XGBClassifier()
xg_class.fit(x_train,y_train)



xgboost in native form stuff
df = pd.read_csv('updatedOrthologsNoCountsCumulativeForLearning.out_humanTrain.out', sep = '\t')
from sklearn.utils import class_weight
class_weights = list(class_weight.compute_class_weight('balanced', np.unique(df['label']),df['label']))
y_train = pd.factorize(df['label'])[0]
for i, val in enumerate(y_train):
     w_array[i] = class_weights[val-1]
df.columns = df.columns.str.replace("[\[,\]<]","")
features = df.columns[1:]
xgb_model = xgb.XGBClassifier(objective='multi:softprob', random_state=42, verbosity = 2, subsample = .8)
xgb_model.fit(df[features], y_train, sample_weight=w_array)
test = pd.read_csv('updatedOrthologsNoCountsCumulativeForLearning.out_humanTest.out', sep = '\t')
test.columns = test.columns.str.replace("[\[,\]<]","")
y_test = pd.factorize(test['label'])[0]
y_pred = xgb_model.predict(test[features])
import matplotlib.pyplot as plt
from sklearn.metrics import plot_confusion_matrix
plot_confusion_matrix(xgb_model, test[features], y_test, normalize = 'all', xticks_rotation = 'vertical', display_labels = labels, include_values = False)
figure = plt.gcf()
figure.set_size_inches(14,14)
plt.savefig('xgbConfusionTest.pdf')
48% accuracy with default values (above)
xgb_model = xgb.XGBClassifier(objective='multi:softprob', random_state=42, verbosity = 2, subsample = .8, nthread = 4, max_depth = 10, n_estimators = 500)
xgb_model = xgb.XGBClassifier(objective='multi:softprob', random_state=42, 
    subsample = .8, nthread = 4, max_depth = 10, 
    n_estimators = 500)
eval_set = [(df[features], y_train), (test[features], y_test)]
xgb_model.fit(df[features], y_train, early_stopping_rounds=10,
    eval_metric=["error", "mlogloss"], eval_set=eval_set, verbose=True)


'''

