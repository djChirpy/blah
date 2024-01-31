import pandas as pd
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
from matplotlib.patches import Patch
from sklearn.manifold import TSNE 
import plotly.express as px
from nancorrmp.nancorrmp import NaNCorrMp

parser = argparse.ArgumentParser()
parser.add_argument('-hc', '--humanCancer', default = 'all')
parser.add_argument('-ag', '--ageAtDiag', default = '/seq/vgb/jason/cancer_mutations/cancer_machine_learning/age_map/sampleAges.tsv')
parser.add_argument('-i', '--inFile', default = '/seq/vgb/swofford/R01/convertingToSKLInput/testVCFs/allSamplesProcessedArtifactSigsRemovedLowPerformingClassesRemoved_080522.txt')
parser.add_argument('-if', '--inFolder', default = '/seq/vgb/swofford/R01/convertingToSKLInput/testVCFs/artifactSigsRemovedNewlineFixed_080522/')
parser.add_argument('-ip', '--inProb', default = '/seq/vgb/swofford/R01/convertingToSKLInput/testVCFs/artifactSigsRemovedNewlineFixed_080522/artifactSigsRemovedNewlineFixed_080522_probs.txt')
parser.add_argument('-fn', '--featureNum', default = 10) #number of features for heatmap
parser.add_argument('-sub', '--subTypes', default = '/seq/vgb/cancer_r01/ML_annotations/all_human_ids_categories_03032022.txt')
parser.add_argument('-t', '--plotType', default = 'both') #other options: heatmap, tsne
parser.add_argument('-ft', '--featureType', default = 'shapley') #other options: feature
parser.add_argument('-tc', '--toCluster', default = 'classified') #other options: truePos, allSamples
parser.add_argument('-o', '--outPrefix', default = 'featureCorrelations') #folder where graphs will be produced
parser.add_argument('-it', '--iterations', default = 10000) #number of tsne iterations

args = parser.parse_args()
args.inFile = 'driverGenesBinary092322/driverGenesBinary092322_featureValues.txt'
df = pd.read_csv(args.inFile, sep = '\t')
labels = pd.unique(df['label'])


features = df.columns[3:-2]


for i in labels:
    test, p_test = NaNCorrMp.calculate_with_p_value(df[df['label'] == i][features], n_jobs = -1)
    test.fillna(0).to_csv('./featureCorrelations/driverBinary/' + i + '_corr.txt', sep = '\t')
    p_test.fillna(1).to_csv('./featureCorrelations/driverBinary/' + i + '_pval.txt', sep = '\t')

for i in labels:
    test, p_test = NaNCorrMp.calculate_with_p_value(df[df['pred'] == i][features], n_jobs = -1)
    test.fillna(0).to_csv('./featureCorrelations/driverBinary/' + i + '_pred_corr.txt', sep = '\t')
    p_test.fillna(1).to_csv('./featureCorrelations/driverBinary/' + i + '_pred_pval.txt', sep = '\t')

for i in labels:
    if not os.path.exists('./featureCorrelations/full/' + i + '_corr.txt'):
        test, p_test = NaNCorrMp.calculate_with_p_value(df[df['label'] == i][features], n_jobs = -1)
        test.fillna(0).to_csv('./featureCorrelations/full/' + i + '_corr.txt', sep = '\t')
        p_test.fillna(1).to_csv('./featureCorrelations/full/' + i + '_pval.txt', sep = '\t')
