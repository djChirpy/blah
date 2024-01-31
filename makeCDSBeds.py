import gffutils
import os
import gzip
import subprocess

#orthFile = '/seq/vgb/swofford/ref/EnsemblHumanDogOrthologs.txt'
#orthFile = '/seq/vgb/swofford/dogHumanOrthologs/combined.bed'
#orthFile = '/seq/vgb/swofford/dogHumanOrthologs/combinedOrthologs.bed'
orthFile = '/seq/vgb/swofford/dogHumanOrthologs/old/feb2/ncbiGeneNameAdded_zeroRemoved_030623.bed'
humanIDs = set()
dogIDs = set()
#zooIDs = set()

with open(orthFile, 'r') as inFile:
    for line in inFile:
        #if 'Ensembl' in line or 'NCBI' in line:
        fields = line.split()
        try:
            meta = fields[3].split(':')
        except:
            print(line)
            continue
        if meta[1] != '':
            dogIDs.add(meta[1].strip())
        if meta[3] != '':
            humanIDs.add(meta[3].strip())
        if meta[5] != '':
            humanIDs.add(meta[5].strip())
        if meta[6] != '':
            dogIDs.add(meta[6].strip())
        '''
        if meta[3] != '' and meta[1] != '':
            humanIDs.add(meta[3])
            dogIDs.add(meta[1])
                
            #elif meta[2] != '':
                #zooIDs.add(meta[2])
            #    humanIDs.add(meta[2])
        
        if 'ortholog_one2one' in line:
            fields = line.split()
            humanIDs.add(fields[0])
            dogIDs.add(fields[2])
        '''

print('human = ' + str(len(humanIDs)))
print('dog = ' + str(len(dogIDs)))
#print('zoo = ' + str(len(zooIDs)))
            

if os.path.exists('/seq/vgb/swofford/ref/hg38.104.db'):
    hg38db = gffutils.FeatureDB('/seq/vgb/swofford/ref/hg38.104.db')
else:
    hg38db = gffutils.create_db('/seq/vgb/swofford/ref/Homo_sapiens.GRCh38.99.chr_patch_hapl_scaff.gtf', 
                            dbfn='hg38.99', 
                            force = True, keep_order = True, 
                            merge_strategy='merge', 
                            sort_attribute_values=True, 
                            disable_infer_genes=True, 
                            disable_infer_transcripts=True)

if os.path.exists('/seq/vgb/swofford/ref/hg19.87.db'):
    hg19db = gffutils.FeatureDB('/seq/vgb/swofford/ref/hg19.87.db')
else:
    hg19db = gffutils.create_db('/seq/vgb/swofford/ref/Homo_sapiens.GRCh37.87.chr_patch_hapl_scaff.gtf', 
                            dbfn='/seq/vgb/swofford/ref/hg19.87.db', 
                            force = True, keep_order = True, 
                            merge_strategy='merge', 
                            sort_attribute_values=True, 
                            disable_infer_genes=True, 
                            disable_infer_transcripts=True)
    
if os.path.exists('/seq/vgb/swofford/ref/canFam3.104.db'):
    canFam3db = gffutils.FeatureDB('/seq/vgb/swofford/ref/canFam3.104.db')
else:
    #canFam = gzip.open('/seq/vgb/swofford/ref/Canis_familiaris.CanFam3.1.99.gtf.gz', 'rb')
    canFam3db = gffutils.create_db('/seq/vgb/swofford/ref/Canis_familiaris.CanFam3.1.99.gtf', 
                            dbfn='/seq/vgb/swofford/ref/canFam3.99.db', 
                            force = True, keep_order = True, 
                            merge_strategy='merge', 
                            sort_attribute_values=True, 
                            disable_infer_genes=True, 
                            disable_infer_transcripts=True)
'''
if os.path.exists('/seq/vgb/swofford/forShirley/dogZoonomiaHigh'):
    dogZoodb = gffutils.FeatureDB('/seq/vgb/swofford/forShirley/dogZoonomiaHigh')
else:
    #canFam = gzip.open('/seq/vgb/swofford/ref/Canis_familiaris.CanFam3.1.99.gtf.gz', 'rb')
    dogZoodb = gffutils.create_db('/seq/vgb/swofford/forShirley/Canis_lupus_familiaris_highQual.gtf', 
                            dbfn='/seq/vgb/swofford/forShirley/dogZoonomiaHigh', 
                            force = True, keep_order = True, 
                            merge_strategy='merge', 
                            sort_attribute_values=True, 
                            disable_infer_genes=True, 
                            disable_infer_transcripts=True)
'''

for gene in humanIDs:
    try:
        name = hg38db[gene]['gene_name'][0]
    except (gffutils.exceptions.FeatureNotFoundError, KeyError) as e:
        try:
            name = hg19db[gene]['gene_name'][0]
        except (gffutils.exceptions.FeatureNotFoundError, KeyError) as e:
            continue
        #name = gene
    try:
        for i in hg38db.children(gene, featuretype='transcript'):
            if i['transcript_biotype'][0] == 'protein_coding':
                cdsStarts = []
                cdsEnds = []
                cdsExonNumber = []
                chrom = str(i.chrom)
                for j in hg38db.children(i, featuretype='CDS', order_by='start'):
                    #print(j)
                    cdsStarts.append(j.start -1)
                    cdsEnds.append(j.end)
                
                    
                with open('hg38CDS.bed', 'a') as outBed:
                    for pair in zip(cdsStarts, cdsEnds):
                        outBed.write('chr' + chrom + '\t' + str(pair[0]) + '\t' 
                              + str(pair[1]) + '\t' + gene + ';' + name + '\n')
    except gffutils.exceptions.FeatureNotFoundError:
        with open('hg38NotFound.txt', 'a') as failOut:
            failOut.write(gene + '\n')
    
    try:        
        for i in hg19db.children(gene, featuretype='transcript'):
            if i['transcript_biotype'][0] == 'protein_coding':
                cdsStarts = []
                cdsEnds = []
                cdsExonNumber = []
                chrom = str(i.chrom)
                for j in hg19db.children(i, featuretype='CDS', order_by='start'):
                    #print(j)
                    cdsStarts.append(j.start -1)
                    cdsEnds.append(j.end)
                
                    
                with open('hg19CDS.bed', 'a') as outBed:
                    for pair in zip(cdsStarts, cdsEnds):
                        outBed.write('chr' + chrom + '\t' + str(pair[0]) + '\t' 
                              + str(pair[1]) + '\t' + gene + ';' + name + '\n')
    except gffutils.exceptions.FeatureNotFoundError:
        with open('hg19NotFound.txt', 'a') as failOut:
            failOut.write(gene + '\n')
            
for gene in dogIDs:
    try:
        name = canFam3db[gene]['gene_name'][0]
    except:
        name = gene
    try:
        for i in canFam3db.children(gene, featuretype='transcript'):
            if i['transcript_biotype'][0] == 'protein_coding':
                cdsStarts = []
                cdsEnds = []
                cdsExonNumber = []
                chrom = str(i.chrom)
                for j in canFam3db.children(i, featuretype='CDS', order_by='start'):
                    #print(j)
                    cdsStarts.append(j.start -1)
                    cdsEnds.append(j.end)
                
                    
                with open('canFam3CDS.bed', 'a') as outBed:
                    for pair in zip(cdsStarts, cdsEnds):
                        outBed.write('chr' + chrom + '\t' + str(pair[0]) + '\t' 
                              + str(pair[1]) + '\t' + gene + ';' + name + '\n')
    
    except:
        with open('canFam3NotFound.txt', 'a') as failOut:
            failOut.write(gene + '\n')
            
os.system('sort -k1,1 -k2,2n canFam3CDS.bed | bedtools merge -i stdin > canFam3CDSMerged.bed')
os.system('sort -k1,1 -k2,2n hg38CDS.bed | bedtools merge -i stdin > hg38CDSMerged.bed')
os.system('sort -k1,1 -k2,2n hg19CDS.bed | bedtools merge -i stdin > hg19CDSMerged.bed')

'''           
for gene in zooIDs:
    geneID = gene + '_Canis_lupus_familiaris_zoonomia'
    with open('zooDbIds.txt', 'a') as zooOut:
        zooOut.write(geneID + '\n')
    
    try:
        name = dogZoodb[geneID]['gene_name'][0]
    except:
        name = geneID
    try:
        for i in dogZoodb.children(geneID, featuretype='transcript'):
            if i['transcript_biotype'][0] == 'protein_coding':
                cdsStarts = []
                cdsEnds = []
                cdsExonNumber = []
                chrom = str(i.chrom)
                for j in dogZoodb.children(i, featuretype='CDS', order_by='start'):
                    #print(j)
                    cdsStarts.append(j.start -1)
                    cdsEnds.append(j.end)
                
                    
                with open('zooDogCDS.bed', 'a') as outBed:
                    for pair in zip(cdsStarts, cdsEnds):
                        outBed.write(chrom + '\t' + str(pair[0]) + '\t' 
                              + str(pair[1]) + '\n')
    
    except:
        with open('zooDogNotFound.txt', 'a') as failOut:
            failOut.write(geneID + '\n')
'''


