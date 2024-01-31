import gffutils
import argparse
import os
import re
import pickle
import statistics as stats

parser = argparse.ArgumentParser()

parser.add_argument('-i', '--inFile', nargs='+', default = '/seq/vgb/swofford/zoonomia/wholeGenomeAnnotation/getGene/fileTest/*.gtf')
args = parser.parse_args()

buscoFile = 'buscoEnsemblIDs.txt'
humdb = db = gffutils.FeatureDB('/seq/vgb/swofford/zoonomia/hg38.99.db')

speciesDB = {}
busco = {}

genes = {}
with open(buscoFile, 'r') as inFile:
    for line in inFile:
        fields = re.sub(r"\s+", '_', line.strip()).split(',')
        if len(fields) == 2:
            genes[fields[0]] = [fields[1]]
        elif len(fields) > 2:
            for i in range(len(fields) - 2):
                genes[fields[i]] = [fields[i]]
            genes[fields[-2]] = [fields[-1]]
        elif len(fields) == 1:
            #print(fields)
            name = humdb[fields[0]]['gene_name'][0]
            try:
                genes[name].append(fields[0])
            except KeyError:
                genes[name] = [fields[0]]

for file in args.inFile:
    name = file.split('/')[-1:][0]
    species = name[:-13]
    if os.path.exists(name + '.db'):
        db = gffutils.FeatureDB(name + '.db')
        speciesDB[species] = db
    else:
        db = gffutils.create_db(file, 
                            dbfn=name + '.db', 
                            force = True, keep_order = True, 
                            merge_strategy='merge', 
                            sort_attribute_values=True, 
                            disable_infer_genes=True, 
                            disable_infer_transcripts=True)
        speciesDB[species] = db
    
    for i in genes.keys():
        ens = genes[i][0]
        if 'ENSG' not in ens:
            continue
        for critter in speciesDB.keys():
            zoo = ens + '_' + critter + '_zoonomia'
            try:
                gene = speciesDB[critter][zoo]
                exons = set()
                for exon in speciesDB[critter].children(zoo, featuretype = 'exon'):
                    exons.add((exon.end, (exon.start - 1)))
                    #exons.append(exon.end - (exon.start - 1))
                exons = list(exons)
                exonLens = []
                for x in exons:
                    exonLens.append(x[0] - x[1])
                exons = exonLens
                try:
                    busco[i]['numExons'].append(len(exons))
                    busco[i]['cdsLen'].append(sum(exons))
                except KeyError:
                    busco[i] = {}
                    busco[i]['numExons'] = len(exons)
                    busco[i]['cdsLen'] = sum(exons)
                #print(gene['gene_name'][0])
                #print(exons)
            except gffutils.exceptions.FeatureNotFoundError:
                pass
                #print('no ' + ens + ' in ' + critter)
    totalExons = 0
    totalCDS = 0
    numGenes = len(busco.keys())
    for i in busco:
        totalExons += busco[i]['numExons'] #errorLine
        totalCDS += busco[i]['cdsLen']
    meanExons = totalExons/len(busco)
    meanCDS = totalCDS/len(busco)
    specName = name[:-13]
    #need number of genes, exons, total coding sequence per genome
    with open(specName + '_summary.txt', 'w') as outFile:
        outFile.write(specName + '\t' + str(numGenes) + '\t' + str(totalExons) 
                      + '\t' + str(meanExons) + '\t' + str(totalCDS) 
                      + '\t' + str(meanCDS) + '\n')
        


#/seq/vgb/swofford/zoonomia/wholeGenomeAnnotation/getGene/fileTest/Ziphius_cavirostris_highQual.gtf


#gtfDir = '/seq/vgb/swofford/zoonomia/wholeGenomeAnnotation/getGene/fileTest/*.gtf'
