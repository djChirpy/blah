'''
To do:
take the top species in each order as a representative. Cycle through and find all genes present in n-1 of the species
as well as all of the species as a whole. Finding a minimum gene set.

Take the minimum gene set and perform go enrichment analysis between it and the full gene set. 
Full gene set can mean:
all input genes
all input genes that succeeded in human
all input genes that have an annotation in any species
all input genes that have an annotation in primates
all input genes that have an annotation in the primate representative

set of genes shared between all primates, but missing from the rest - I may already have this
'''



import argparse
import os
#from goatools.go_enrichment import GOEnrichmentStudy
from goatools.goea.go_enrichment_ns import GOEnrichmentStudyNS
from goatools.base import download_go_basic_obo
from goatools.base import download_ncbi_associations
from goatools.obo_parser import GODag
from goatools.anno.genetogo_reader import Gene2GoReader
from genes_ncbi_9606_proteincoding import GENEID2NT as GeneID2nt_hum
from goatools.godag_plot import plot_gos, plot_results, plot_goid2goobj

#from goatools import obo_parser
import wget
import Bio.UniProt.GOA as GOA
import gzip
import statsmodels


parser = argparse.ArgumentParser()
parser.add_argument('-s', '--statsFile', default = 'conservedGenesGenomeStatsNew.txt')
parser.add_argument('-t', '--tax', default = '/seq/vgb/swofford/ref/kraken2/coronaVirusProject/taxonomy/names.dmp')

args = parser.parse_args()

class species:
    def __init__(self, name, orderID, familyID, genusID, clade, contigN50, scaffoldN50, HQGeneCount, branchLength):
        
        self.name = name
        #self.tribe = taxonLookup(tribeID)
        self.clade = clade
        #self.superorder = taxonLookup(superorderID)
        self.order = taxonLookup(orderID)
        #self.infraorder = taxonLookup(infraorderID)
        #self.suborder = taxonLookup(suborderID)
        self.family = taxonLookup(familyID)
        self.genus = taxonLookup(genusID)
        self.contigN50 = contigN50
        self.scaffoldN50 = scaffoldN50
        self.HQGeneCount = HQGeneCount
        self.branchLength = branchLength
        self.HQGenes = []

def taxonLookup(taxon):
    return(taxons[taxon])

cladeSpecies = {}
orderSpecies = {}
familySpecies = {}
genusSpecies = {}


taxons = {}
speciesList = {}
geneSet = set()
geneNames = {}
topSpecies = {}
ncbiToEns ={}
ensToNCBI = {}
ensToUni = {}

with open('/seq/vgb/swofford/zoonomia/conservedGenes/bestOfOrder.txt', 'r') as inFile:
    next(inFile)
    for line in inFile:
        fields = line.split('\t')
        topSpecies[fields[0].strip()] = fields[1].strip()

'''
url = 'http://purl.obolibrary.org/obo/go/go-basic.obo'
goFolder = '/seq/vgb/swofford/ref/GO/data'

if not os.path.exists('/seq/vgb/swofford/ref/GO'):
    os.mkdir('/seq/vgb/swofford/ref/GO')
if not os.path.exists(goFolder):
    os.mkdir(goFolder)

if not os.path.exists(goFolder + '/go-basic.obo'):
    go_obo = wget.download(url, goFolder + '/go-basic.obo')
else:
    go_obo = goFolder + '/go-basic.obo'
    
go = obo_parser.GODag(go_obo)
'''
obo_fname = download_go_basic_obo()
fin_gene2go = download_ncbi_associations()
obodag = GODag("go-basic.obo")
objanno = Gene2GoReader(fin_gene2go, taxids=[9606])

ns2assoc = objanno.get_ns2assc()

for nspc, id2gos in ns2assoc.items():
    print("{NS} {N:,} annotated human genes".format(NS=nspc, N=len(id2gos)))
    
print(len(GeneID2nt_hum))

#print(GeneID2nt_hum.keys())


'''
with gzip.open('/seq/vgb/swofford/ref/GO/goa_human.gaf.gz', 'rt') as gafFile:
    hum_funcs = {}  # Initialise the dictionary of functions
    
    # Iterate on each function using Bio.UniProt.GOA library.
    for entry in GOA.gafiterator(gafFile):
        uniprot_id = entry.pop('DB_Object_ID')
        hum_funcs[uniprot_id] = entry
        
pop = hum_funcs.keys()

assoc = {}

for x in hum_funcs:
    if x not in assoc:
        assoc[x] = set()
    assoc[x].add(str(hum_funcs[x]['GO_ID']))
print(assoc)
'''   

with open('ensemblToNCBI.txt', 'r') as inFile:
    
    for line in inFile:
        fields = line.split('\t')
        try:
            ensID = fields[0].strip()
            ncbi = int(fields[5].strip())
        except (IndexError, ValueError):
            print('no ncbi entry for ' + ensID + '\n')
            continue
        try:
            ensToNCBI[ensID].add(ncbi)
            ncbiToEns[ncbi] = ensID
        except KeyError:
            ensToNCBI[ensID] = set()
            ensToNCBI[ensID].add(ncbi)
            ncbiToEns[ncbi] = ensID
            
with gzip.open('/seq/vgb/swofford/ref/ensemblToUniprotWithGOTermsAllHumanGenesFrom99.gz', 'rt') as inFile:
    
    #assoc = {}
    for line in inFile:
        fields = line.split('\t')
        uniProt = fields[0].strip()
        ensID = fields[12].strip()
        #GO = fields[10].split(';')
        #goSet = set()
        #for i in GO:
        #    goSet.add(i.strip())
        ensToUni[ensID] = uniProt
        #assoc[uniProt] = goSet

#methods = ["bonferroni", "sidak", "holm", "fdr"]

#g = GOEnrichmentStudy(pop, assoc, go,
#                         propagate_counts=True,
#                         alpha=0.05,
#                         methods=methods)        



with open('geneNames.txt', 'r') as inGene:
    for line in inGene:
        fields = line.split('\t')
        name = fields[1].strip()
        id = fields[0]
        geneNames[id] = name
        geneSet.add(id)

with open(args.tax, 'r') as inTax:
    for line in inTax:
        fields = line.split('|')
        taxon = fields[0].strip()
        name = fields[1].strip()
        type = fields[3].strip()
        if type == 'scientific name':
            taxons[taxon] = name

              
n = 0
with open(args.statsFile, 'r') as inStat, open(args.statsFile + '.names.txt', 'w') as outFile:
    for line in inStat:
        if n == 0:
            outFile.write(line)
            n = 1
            continue
        fields = line.split('\t')
        name = fields[0].strip()
        orderID = fields[1].strip()
        familyID = fields[2].strip()
        genusID = fields[3].strip()
        clade = fields[4].strip()
        contigN50 = fields[5].strip()
        scaffoldN50 = fields[6].strip()
        HQGenes = fields[7].strip()
        branchLength = fields[8].strip()
        newSpecies = species(name = name, 
                                orderID = orderID, 
                                familyID = familyID, 
                                genusID = genusID, 
                                clade = clade, 
                                contigN50 = contigN50, 
                                scaffoldN50 = scaffoldN50, 
                                HQGeneCount = HQGenes, 
                                branchLength = branchLength)
        speciesList[name] = newSpecies
        outLine = [name, taxons[orderID], taxons[familyID], taxons[genusID]]
        outLine = outLine + fields[4:]
        outFile.write('\t'.join(outLine))
        
### creating lists of species per taxonomic group to be able
### to assign levels later           
for species in speciesList.keys():
    clade = speciesList[species].clade
    order = speciesList[species].order
    genus = speciesList[species].genus
    family = speciesList[species].family
    fileName = 'speciesGenes/' + species + '.geneList'
    with open(fileName, 'r') as inFile:
        for line in inFile:
            speciesList[species].HQGenes.append(line.strip())
            #geneSet.add(line.strip())
    try:
        cladeSpecies[clade].append(species)
    except KeyError:
        cladeSpecies[clade] = [species]
    try:
        orderSpecies[order].append(species)
    except KeyError:
        orderSpecies[order] = [species]
    try:
        genusSpecies[genus].append(species)
    except KeyError:
        genusSpecies[genus] = [species]
    try:
        familySpecies[family].append(species)
    except KeyError:
        familySpecies[family] = [species]
    
genesWithAnn = set()
genesWithoutAnn = set()
goGeneSet = set()

for species in speciesList.keys():
    for gene in speciesList[species].HQGenes:
        genesWithAnn.add(gene)
        #print(gene)
        #print(list(ensToNCBI[gene])[0])
        try:
            goGeneSet.add(list(ensToNCBI[gene])[0])
        except KeyError:
            continue
goGeneSet = list(goGeneSet)

with open('genesPresentInAllSpecies.txt' ,'w') as outFile:
    for gene in geneSet:
        inAll = True
        for species in speciesList.keys():
            if inAll:
                if gene not in speciesList[species].HQGenes:
                    inAll = False
        if inAll:
            outFile.write(gene + '\n')      
    


with open('fullGeneSetNames.txt', 'w') as outFile:
    for i in list(genesWithAnn):
        outFile.write(geneNames[i] + '\n')


print("subset: " + str(len(goGeneSet)))
print("total: " + str(len(GeneID2nt_hum.keys())))

for gene in geneSet:
    if gene not in genesWithAnn:
        genesWithoutAnn.add(gene) 
        

        
goeaobj = GOEnrichmentStudyNS(
        #GeneID2nt_hum.keys(), # List of human protein-coding genes
        goGeneSet,
        ns2assoc, # geneid/GO associations
        obodag, # Ontologies
        propagate_counts = False,
        alpha = 0.05, # default significance cut-off
        methods = ['fdr_bh']) # defult multipletest correction method               
    
if not os.path.exists('output'):
    os.mkdir('output')
if not os.path.exists('output/clade'):
    os.mkdir('output/clade')
if not os.path.exists('output/genus'):
    os.mkdir('output/genus')
if not os.path.exists('output/order'):
    os.mkdir('output/order')
if not os.path.exists('output/family'):
    os.mkdir('output/family')


### creating a list of genes annotated in at least one species
### in a given taxonomic level. For use in questions like
### "is this gene present in this clade"
cladeGeneSets = {}   
orderGeneSets = {}
genusGeneSets = {}
orderGenesEns = {}


'''Start adding in the stuff referenced at the top here. Things requested from Kerstin/Elinor meeting'''
bestOfOrder = {}

for order in orderSpecies.keys():
    if order != 'root':
        best = ''
        for critter in orderSpecies[order]:
            if best == '':
                best = speciesList[critter]
            elif len(speciesList[critter].HQGenes) > len(best.HQGenes):
                best = speciesList[critter]
        bestOfOrder[order] = speciesList[critter]

with open('output/bestOfOrderSummary.txt', 'w') as outFile:
    orders = list(bestOfOrder.keys())
    genesPresent = {}
    goList = []
    for i in range(len(orders)):
        outcast = orders[i]
        os.system('rm output/orderGenesMinus_' + outcast)
        os.system('rm without_' + outcast + '.png')
        for species in bestOfOrder.values():
            if species.order != outcast:
                for gene in species.HQGenes:
                    try:
                        genesPresent[gene] += 1
                    except KeyError:
                        genesPresent[gene] = 1
        for gene in genesPresent.keys():
            if genesPresent[gene] == 17:
                #outFile.write(outcast + '\t' + str(len(genesPresent.keys())) + '\n')
                with open('output/orderGenesMinus_' + outcast, 'a') as genesOut:
                    #for gene in genesPresent:
                    genesOut.write(gene + '\n')
                    try:
                        goList.append(list(ensToNCBI[gene])[0])
                    except KeyError:
                        continue
                    
        outcastpng = 'output/without_' + outcast + '.png'            
        g_pres = goeaobj.run_study(goList)
        g_presSig = [r for r in g_pres if r.p_fdr_bh < 0.05]
        plot_results(outcastpng, g_presSig)
        genesPresent = {}
        goList = []
        
#genes present in all of the best representatives
with open('output/bestOfOrderSharedGenes.txt', 'w') as outFile:
    orders = list(bestOfOrder.keys())
    genesPresent = {}
    goList = []
    for species in bestOfOrder.values():
        for gene in species.HQGenes:
            try:
                genesPresent[gene] += 1
            except KeyError:
                genesPresent[gene] = 1
    for gene in genesPresent.keys():
            if genesPresent[gene] == len(orders):
                #with open('output/bestOfOrderAllGenes.txt', 'a') as genesOut:
                    #for gene in genesPresent:
                outFile.write(gene + '\n')
                try:
                    goList.append(list(ensToNCBI[gene])[0])
                except KeyError:
                    continue
    bestOfOrderpng = 'output/bestOfOrderAllGenes.png'            
    g_pres = goeaobj.run_study(goList)
    g_presSig = [r for r in g_pres if r.p_fdr_bh < 0.05]
    with open('output/bestOfOrderAllGenesGO.txt', 'w') as goOut:
        for i in g_presSig:
            print(i)
    plot_results(bestOfOrderpng, g_presSig)
    genesPresent = {}
    goList = []
    
    
for i in range(len(orders)):
    outcast = orders[i]
    os.system('rm output/orderGenesMinus_' + outcast)
    os.system('rm without_' + outcast + '.png')
    for species in bestOfOrder.values():
        if species.order != outcast:
            for gene in species.HQGenes:
                try:
                    genesPresent[gene] += 1
                except KeyError:
                    genesPresent[gene] = 1
    for gene in genesPresent.keys():
        if genesPresent[gene] == 17:
            #outFile.write(outcast + '\t' + str(len(genesPresent.keys())) + '\n')
            with open('output/orderGenesMinus_' + outcast, 'a') as genesOut:
                #for gene in genesPresent:
                genesOut.write(gene + '\n')
                try:
                    goList.append(list(ensToNCBI[gene])[0])
                except KeyError:
                    continue
                
    outcastpng = 'output/without_' + outcast + '.png'            
    g_pres = goeaobj.run_study(goList)
    g_presSig = [r for r in g_pres if r.p_fdr_bh < 0.05]
    plot_results(outcastpng, g_presSig)
    genesPresent = {}
    goList = []
        
'''end of new stuff, I think I need to move the go declarations around so that I can do enrichment on these gene sets, which are small'''
            

for species in speciesList.keys():
    for gene in speciesList[species].HQGenes:
        try:
            orderGenesEns[speciesList[species].order].add(gene)
        except KeyError:
            orderGenesEns[speciesList[species].order] = set()
            orderGenesEns[speciesList[species].order].add(gene)
        try: 
            ncbiGene = list(ensToNCBI[gene])[0]
        except KeyError:
            continue
        try:
            cladeGeneSets[speciesList[species].clade].add(ncbiGene)
        except KeyError:
            try:
                cladeGeneSets[speciesList[species].clade] = set()
                cladeGeneSets[speciesList[species].clade].add(ncbiGene)
            except KeyError:
                continue
        try:
            orderGeneSets[speciesList[species].order].add(ncbiGene)
        except KeyError:
            try:
                orderGeneSets[speciesList[species].order] = set()
                orderGeneSets[speciesList[species].order].add(ncbiGene)
            except KeyError:
                continue
        try:
            genusGeneSets[speciesList[species].genus].add(ncbiGene)
        except KeyError:
            try:
                genusGeneSets[speciesList[species].genus] = set()
                genusGeneSets[speciesList[species].genus].add(ncbiGene)
            except KeyError:
                continue
for i in cladeGeneSets.keys():
    #print(len(cladeGeneSets[i]))
    #print(cladeGeneSets[i])
    cladeGeneSets[i] = list(cladeGeneSets[i])
for i in orderGeneSets.keys():
    orderGeneSets[i] = list(orderGeneSets[i])    
for i in genusGeneSets.keys():
    genusGeneSets[i] = list(genusGeneSets[i])    

orderGeneCounts = {}
for order in  orderGenesEns.keys():
    if order != 'root':
        for gene in list(orderGenesEns[order]):
            try:
                orderGeneCounts[gene] += 1
            except KeyError:
                orderGeneCounts[gene] = 1

with open('output/anyOfOrderSharedGenes.txt', 'w') as outFile:
    #print(str(len(orderGeneCounts.keys()) - 1) + ' orders')
    for gene in orderGeneCounts.keys():
        print(gene)
        if orderGeneCounts[gene] == len(orders):
            outFile.write(gene + '\n')
        else:
            print(orderGeneCounts[gene])
    
    

    
'''
for clade in cladeSpecies.keys():
    if not os.path.exists('output/clade/' + clade):
        os.mkdir('output/clade/' + clade)

    speciesCount = len(cladeSpecies[clade])
    genePercentage = {}
    for i in genesWithAnn:
        genePercentage[i] = 0
    for species in cladeSpecies[clade]:
        for gene in speciesList[species].HQGenes:
            genePercentage[gene] += 1
    
    percentFile = 'output/clade/' + clade + '/percentPresent.txt'
    missingFile = 'output/clade/' + clade + '/missingInGroup.txt'
    presentFile = 'output/clade/' + clade + '/presentInGroup.txt'
    present100File = 'output/clade/' + clade + '/presentIn100Group.txt'
    present95File = 'output/clade/' + clade + '/presentIn95Group.txt'
    present90File = 'output/clade/' + clade + '/presentIn90Group.txt'
    missing95File = 'output/clade/' + clade + '/missingIn95Group.txt'
    missing90File = 'output/clade/' + clade + '/missingIn90Group.txt'
    missing90png = 'output/clade/' + clade + '/missingIn90GO.png'
    present90png = 'output/clade/' + clade + '/presentIn90Go.png'
    missing100png = 'output/clade/' + clade + '/missingIn100GO.png'
    present100png = 'output/clade/' + clade + '/presentIn100Go.png'
    
    #create list of uniprot IDs for missing/present in 90+ for GO enricment
    missingIn90 = []
    presentIn90 = []
    missingIn100 = []
    presentIn100 = []
    
    with open(percentFile, 'w') as percentOut, \
        open(missingFile, 'w') as missingOut, \
        open(presentFile, 'w') as presentOut, \
        open(present100File, 'w') as present100Out, \
        open(present95File, 'w') as present95Out, \
        open(present90File, 'w') as present90Out, \
        open(missing95File, 'w') as missing95Out, \
        open(missing90File, 'w') as missing90Out:
            
        for i in genePercentage.keys():
            percentage = genePercentage[i]/speciesCount
            percentOut.write(i + '\t' + geneNames[i] + '\t' + str(percentage) + '\n')
            if percentage > 0:
                presentOut.write(i + '\t' + geneNames[i] + '\n')
            if percentage == 1:
                present100Out.write(i + '\t' + geneNames[i] + '\n')
                try:
                    for ncbi in ensToNCBI[i]:
                        presentIn100.append(ncbi)
                except KeyError:
                    pass
            if percentage == 0:
                missingOut.write(i + '\t' + geneNames[i] + '\n')
                try:
                    for ncbi in ensToNCBI[i]:
                        missingIn100.append(ncbi)
                except KeyError:
                    pass
            if percentage >= .9:
                try:
                    #presentIn90.append(i)
                    for ncbi in ensToNCBI[i]:
                        presentIn90.append(ncbi)
                except KeyError:
                    pass
                    #print(i + ' is not found in uniprot\n')
                present90Out.write(i + '\t' + geneNames[i] + '\n')
            if percentage >= .95:
                present95Out.write(i + '\t' + geneNames[i] + '\n')
            if percentage <= .1:
                try:
                    for ncbi in ensToNCBI[i]:
                        missingIn90.append(ncbi)
                    #missingIn90.append(x for x in ensToNCBI[i])
                except KeyError:
                    pass
                    #print(i + ' is not found in uniprot\n')
                missing90Out.write(i + '\t' + geneNames[i] + '\n')
            if percentage <= .05:
                missing95Out.write(i + '\t' + geneNames[i] + '\n')

    print(clade + '\n')
    print('missing in 90% of species:\n')
    print(missingIn90)
    g_mis90 = goeaobj.run_study(missingIn90)
    g_mis90Sig = [r for r in g_mis90 if r.p_fdr_bh < 0.05]
    plot_results(missing90png, g_mis90Sig)
    print('present in 90% of species:\n')
    g_pres90 = goeaobj.run_study(presentIn90)
    g_pres90Sig = [r for r in g_pres90 if r.p_fdr_bh < 0.05]
    plot_results(present90png, g_pres90Sig)


    print('missing in 100% of species:\n')
    print(missingIn100)
    g_mis100 = goeaobj.run_study(missingIn100)
    g_mis100Sig = [r for r in g_mis100 if r.p_fdr_bh < 0.05]
    plot_results(missing100png, g_mis100Sig)
    print('present in 90% of species:\n')
    g_pres100 = goeaobj.run_study(presentIn100)
    g_pres100Sig = [r for r in g_pres100 if r.p_fdr_bh < 0.05]
    plot_results(present100png, g_pres100Sig)
'''

for order in orderSpecies.keys():
    parent = speciesList[orderSpecies[order][0]].clade
    privateGoGeneSet = cladeGeneSets[parent]
    
    goeaobjPrivate = GOEnrichmentStudyNS(
        #GeneID2nt_hum.keys(), # List of human protein-coding genes
        privateGoGeneSet,
        ns2assoc, # geneid/GO associations
        obodag, # Ontologies
        propagate_counts = False,
        alpha = 0.05, # default significance cut-off
        methods = ['fdr_bh']) # defult multipletest correction method     
    
    '''
    orderGeneSets = {}    
    for species in speciesList.keys():
        if species.order = order:
            for gene in speciesList[species].HQGenes:
                genesWithAnn.add(gene)
        #print(gene)
        #print(list(ensToNCBI[gene])[0])
            try:
                goGeneSet.add(list(ensToNCBI[gene])[0])
            except KeyError:
                continue
    goGeneSet = list(goGeneSet)        
    '''
    
    if not os.path.exists('output/order/' + order):
        os.mkdir('output/order/' + order)

    speciesCount = len(orderSpecies[order])
    genePercentage = {}
    for i in genesWithAnn:
        genePercentage[i] = 0
    for species in orderSpecies[order]:
        for gene in speciesList[species].HQGenes:
            genePercentage[gene] += 1
    
    percentFile = 'output/order/' + order + '/percentPresent.txt'
    missingFile = 'output/order/' + order + '/missingInGroup.txt'
    presentFile = 'output/order/' + order + '/presentInGroup.txt'
    present100File = 'output/order/' + order + '/presentIn100Group.txt'
    present95File = 'output/order/' + order + '/presentIn95Group.txt'
    present90File = 'output/order/' + order + '/presentIn90Group.txt'
    missing95File = 'output/order/' + order + '/missingIn95Group.txt'
    missing90File = 'output/order/' + order + '/missingIn90Group.txt'
    missing90png = 'output/order/' + order + '/missingIn90GO.png'
    present90png = 'output/order/' + order + '/presentIn90Go.png'
    missing100png = 'output/order/' + order + '/missingIn100GO.png'
    present100png = 'output/order/' + order + '/presentIn100Go.png'
    missing100pngPrivate = 'output/order/' + order + '/missingIn100PrivGo.png'
    present100pngPrivate = 'output/order/' + order + '/presentIn100PrivGo.png'
    #missingFromParent = 'output/order/' + order + '/missingFromParent.txt'
    
    #create list of uniprot IDs for missing/present in 90+ for GO enricment
    missingIn90 = []
    presentIn90 = []
    missingIn100 = []
    presentIn100 = []
    
    with open(percentFile, 'w') as percentOut, \
        open(missingFile, 'w') as missingOut, \
        open(presentFile, 'w') as presentOut, \
        open(present100File, 'w') as present100Out, \
        open(present95File, 'w') as present95Out, \
        open(present90File, 'w') as present90Out, \
        open(missing95File, 'w') as missing95Out, \
        open(missing90File, 'w') as missing90Out:
            
        for i in genePercentage.keys():
            percentage = genePercentage[i]/speciesCount
            percentOut.write(i + '\t' + geneNames[i] + '\t' + str(percentage) + '\n')
            if percentage > 0:
                presentOut.write(i + '\t' + geneNames[i] + '\n')
            if percentage == 1:
                present100Out.write(i + '\t' + geneNames[i] + '\n')
                try:
                    for ncbi in ensToNCBI[i]:
                        presentIn100.append(ncbi)
                except KeyError:
                    pass
            if percentage == 0:
                missingOut.write(i + '\t' + geneNames[i] + '\n')
                try:
                    for ncbi in ensToNCBI[i]:
                        missingIn100.append(ncbi)
                except KeyError:
                    pass
            if percentage >= .9:
                try:
                    #presentIn90.append(i)
                    for ncbi in ensToNCBI[i]:
                        presentIn90.append(ncbi)
                except KeyError:
                    pass
                    #print(i + ' is not found in uniprot\n')
                present90Out.write(i + '\t' + geneNames[i] + '\n')
            if percentage >= .95:
                present95Out.write(i + '\t' + geneNames[i] + '\n')
            if percentage <= .1:
                try:
                    for ncbi in ensToNCBI[i]:
                        missingIn90.append(ncbi)
                    #missingIn90.append(x for x in ensToNCBI[i])
                except KeyError:
                    pass
                    #print(i + ' is not found in uniprot\n')
                missing90Out.write(i + '\t' + geneNames[i] + '\n')
            if percentage <= .05:
                missing95Out.write(i + '\t' + geneNames[i] + '\n')

    print(order + '\n')
    print('missing in 90% of species:\n')
    #print(missingIn90)
    g_mis90 = goeaobj.run_study(missingIn90)
    g_mis90Sig = [r for r in g_mis90 if r.p_fdr_bh < 0.05]
    plot_results(missing90png, g_mis90Sig)
    '''
    for i in g_misRes:
        if i.p_bonferroni <= 0.05:
            print(i.goterm.id, i.p_bonferroni)
    
    try:
        g.print_results(g_misRes, min_ratio=None, pval=0.01)
    except:
        print('arg')
    '''
    print('present in 90% of species:\n')
    g_pres90 = goeaobj.run_study(presentIn90)
    g_pres90Sig = [r for r in g_pres90 if r.p_fdr_bh < 0.05]
    plot_results(present90png, g_pres90Sig)


    print('missing in 100% of species:\n')
    #print(missingIn100)
    g_mis100 = goeaobj.run_study(missingIn100)
    g_mis100Sig = [r for r in g_mis100 if r.p_fdr_bh < 0.05]
    plot_results(missing100png, g_mis100Sig)
    '''
    for i in g_misRes:
        if i.p_bonferroni <= 0.05:
            print(i.goterm.id, i.p_bonferroni)
    
    try:
        g.print_results(g_misRes, min_ratio=None, pval=0.01)
    except:
        print('arg')
    '''
    print('present in 100% of species:\n')
    g_pres100 = goeaobj.run_study(presentIn100)
    g_pres100Sig = [r for r in g_pres100 if r.p_fdr_bh < 0.05]
    plot_results(present100png, g_pres100Sig)

    print('fraction from parent:\n')
    g_mis100Priv = goeaobjPrivate.run_study(missingIn100)
    g_mis100SigPriv = [r for r in g_mis100Priv if r.p_fdr_bh < 0.05]
    plot_results(missing100pngPrivate, g_mis100SigPriv)
    


'''
for family in familySpecies.keys():
    parent = speciesList[familySpecies[family][0]].genus
    privateGoGeneSet = genusGeneSets[parent]
    
    
    goeaobjPrivate = GOEnrichmentStudyNS(
        #GeneID2nt_hum.keys(), # List of human protein-coding genes
        privateGoGeneSet,
        ns2assoc, # geneid/GO associations
        obodag, # Ontologies
        propagate_counts = False,
        alpha = 0.05, # default significance cut-off
        methods = ['fdr_bh']) # defult multipletest correction method   
    
    if not os.path.exists('output/family/' + family):
        os.mkdir('output/family/' + family)

    speciesCount = len(familySpecies[family])
    genePercentage = {}
    for i in genesWithAnn:
        genePercentage[i] = 0
    for species in familySpecies[family]:
        for gene in speciesList[species].HQGenes:
            genePercentage[gene] += 1
    
    percentFile = 'output/family/' + family + '/percentPresent.txt'
    missingFile = 'output/family/' + family + '/missingInGroup.txt'
    presentFile = 'output/family/' + family + '/presentInGroup.txt'
    present100File = 'output/family/' + family + '/presentIn100Group.txt'
    present95File = 'output/family/' + family + '/presentIn95Group.txt'
    present90File = 'output/family/' + family + '/presentIn90Group.txt'
    missing95File = 'output/family/' + family + '/missingIn95Group.txt'
    missing90File = 'output/family/' + family + '/missingIn90Group.txt'
    missing90png = 'output/family/' + family + '/missingIn90GO.png'
    present90png = 'output/family/' + family + '/presentIn90Go.png'
    missing100png = 'output/family/' + family + '/missingIn100GO.png'
    present100png = 'output/family/' + family + '/presentIn100Go.png'
    missing100pngPrivate = 'output/family/' + family + '/missingIn100PrivGo.png'
    present100pngPrivate = 'output/family/' + family + '/presentIn100PrivGo.png'
    missingFromParent = ''
    
    
    #create list of uniprot IDs for missing/present in 90+ for GO enricment
    missingIn90 = []
    presentIn90 = []
    missingIn100 = []
    presentIn100 = []
    
    with open(percentFile, 'w') as percentOut, \
        open(missingFile, 'w') as missingOut, \
        open(presentFile, 'w') as presentOut, \
        open(present100File, 'w') as present100Out, \
        open(present95File, 'w') as present95Out, \
        open(present90File, 'w') as present90Out, \
        open(missing95File, 'w') as missing95Out, \
        open(missing90File, 'w') as missing90Out:
            
        for i in genePercentage.keys():
            percentage = genePercentage[i]/speciesCount
            percentOut.write(i + '\t' + geneNames[i] + '\t' + str(percentage) + '\n')
            if percentage > 0:
                presentOut.write(i + '\t' + geneNames[i] + '\n')
            if percentage == 1:
                present100Out.write(i + '\t' + geneNames[i] + '\n')
                try:
                    for ncbi in ensToNCBI[i]:
                        presentIn100.append(ncbi)
                except KeyError:
                    pass
            if percentage == 0:
                missingOut.write(i + '\t' + geneNames[i] + '\n')
                try:
                    for ncbi in ensToNCBI[i]:
                        missingIn100.append(ncbi)
                except KeyError:
                    pass
            if percentage >= .9:
                try:
                    #presentIn90.append(i)
                    for ncbi in ensToNCBI[i]:
                        presentIn90.append(ncbi)
                except KeyError:
                    pass
                    #print(i + ' is not found in uniprot\n')
                present90Out.write(i + '\t' + geneNames[i] + '\n')
            if percentage >= .95:
                present95Out.write(i + '\t' + geneNames[i] + '\n')
            if percentage <= .1:
                try:
                    for ncbi in ensToNCBI[i]:
                        missingIn90.append(ncbi)
                    #missingIn90.append(x for x in ensToNCBI[i])
                except KeyError:
                    pass
                    #print(i + ' is not found in uniprot\n')
                missing90Out.write(i + '\t' + geneNames[i] + '\n')
            if percentage <= .05:
                missing95Out.write(i + '\t' + geneNames[i] + '\n')

    print(family + '\n')
    print('missing in 90% of species:\n')
    #print(missingIn90)
    g_mis90 = goeaobj.run_study(missingIn90)
    g_mis90Sig = [r for r in g_mis90 if r.p_fdr_bh < 0.05]
    plot_results(missing90png, g_mis90Sig)

    print('present in 90% of species:\n')
    g_pres90 = goeaobj.run_study(presentIn90)
    g_pres90Sig = [r for r in g_pres90 if r.p_fdr_bh < 0.05]
    plot_results(present90png, g_pres90Sig)


    print('missing in 100% of species:\n')
    #print(missingIn100)
    g_mis100 = goeaobj.run_study(missingIn100)
    g_mis100Sig = [r for r in g_mis100 if r.p_fdr_bh < 0.05]
    plot_results(missing100png, g_mis100Sig)

    print('present in 90% of species:\n')
    g_pres100 = goeaobj.run_study(presentIn100)
    g_pres100Sig = [r for r in g_pres100 if r.p_fdr_bh < 0.05]
    plot_results(present100png, g_pres100Sig)
    
    
    print('fraction from parent:\n')
    g_mis100Priv = goeaobjPrivate.run_study(missingIn100)
    g_mis100SigPriv = [r for r in g_mis100Priv if r.p_fdr_bh < 0.05]
    plot_results(missing100pngPrivate, g_mis100SigPriv)



for genus in genusSpecies.keys():
    
    parent = speciesList[genusSpecies[genus][0]].order
    privateGoGeneSet = orderGeneSets[parent]
    
    goeaobjPrivate = GOEnrichmentStudyNS(
        #GeneID2nt_hum.keys(), # List of human protein-coding genes
        privateGoGeneSet,
        ns2assoc, # geneid/GO associations
        obodag, # Ontologies
        propagate_counts = False,
        alpha = 0.05, # default significance cut-off
        methods = ['fdr_bh']) # defult multipletest correction method   
    
    if not os.path.exists('output/genus/' + genus):
        os.mkdir('output/genus/' + genus)

    speciesCount = len(genusSpecies[genus])
    genePercentage = {}
    for i in genesWithAnn:
        genePercentage[i] = 0
    for species in genusSpecies[genus]:
        for gene in speciesList[species].HQGenes:
            genePercentage[gene] += 1
    
    percentFile = 'output/genus/' + genus + '/percentPresent.txt'
    missingFile = 'output/genus/' + genus + '/missingInGroup.txt'
    presentFile = 'output/genus/' + genus + '/presentInGroup.txt'
    present100File = 'output/genus/' + genus + '/presentIn100Group.txt'
    present95File = 'output/genus/' + genus + '/presentIn95Group.txt'
    present90File = 'output/genus/' + genus + '/presentIn90Group.txt'
    missing95File = 'output/genus/' + genus + '/missingIn95Group.txt'
    missing90File = 'output/genus/' + genus + '/missingIn90Group.txt'
    missing90png = 'output/genus/' + genus + '/missingIn90GO.png'
    present90png = 'output/genus/' + genus + '/presentIn90Go.png'
    missing100png = 'output/genus/' + genus + '/missingIn100GO.png'
    present100png = 'output/genus/' + genus + '/presentIn100Go.png'
    missing100pngPrivate = 'output/genus/' + genus + '/missingIn100PrivGo.png'
    present100pngPrivate = 'output/genus/' + genus + '/presentIn100PrivGo.png'
    
    #create list of uniprot IDs for missing/present in 90+ for GO enricment
    missingIn90 = []
    presentIn90 = []
    missingIn100 = []
    presentIn100 = []
    
    with open(percentFile, 'w') as percentOut, \
        open(missingFile, 'w') as missingOut, \
        open(presentFile, 'w') as presentOut, \
        open(present100File, 'w') as present100Out, \
        open(present95File, 'w') as present95Out, \
        open(present90File, 'w') as present90Out, \
        open(missing95File, 'w') as missing95Out, \
        open(missing90File, 'w') as missing90Out:
            
        for i in genePercentage.keys():
            percentage = genePercentage[i]/speciesCount
            percentOut.write(i + '\t' + geneNames[i] + '\t' + str(percentage) + '\n')
            if percentage > 0:
                presentOut.write(i + '\t' + geneNames[i] + '\n')
            if percentage == 1:
                present100Out.write(i + '\t' + geneNames[i] + '\n')
                try:
                    for ncbi in ensToNCBI[i]:
                        presentIn100.append(ncbi)
                except KeyError:
                    pass
            if percentage == 0:
                missingOut.write(i + '\t' + geneNames[i] + '\n')
                try:
                    for ncbi in ensToNCBI[i]:
                        missingIn100.append(ncbi)
                except KeyError:
                    pass
            if percentage >= .9:
                try:
                    #presentIn90.append(i)
                    for ncbi in ensToNCBI[i]:
                        presentIn90.append(ncbi)
                except KeyError:
                    pass
                    #print(i + ' is not found in uniprot\n')
                present90Out.write(i + '\t' + geneNames[i] + '\n')
            if percentage >= .95:
                present95Out.write(i + '\t' + geneNames[i] + '\n')
            if percentage <= .1:
                try:
                    for ncbi in ensToNCBI[i]:
                        missingIn90.append(ncbi)
                    #missingIn90.append(x for x in ensToNCBI[i])
                except KeyError:
                    pass
                    #print(i + ' is not found in uniprot\n')
                missing90Out.write(i + '\t' + geneNames[i] + '\n')
            if percentage <= .05:
                missing95Out.write(i + '\t' + geneNames[i] + '\n')

    print(genus + '\n')
    print('missing in 90% of species:\n')
    #print(missingIn90)
    g_mis90 = goeaobj.run_study(missingIn90)
    g_mis90Sig = [r for r in g_mis90 if r.p_fdr_bh < 0.05]
    plot_results(missing90png, g_mis90Sig)

    print('present in 90% of species:\n')
    g_pres90 = goeaobj.run_study(presentIn90)
    g_pres90Sig = [r for r in g_pres90 if r.p_fdr_bh < 0.05]
    plot_results(present90png, g_pres90Sig)


    print('missing in 100% of species:\n')
    #print(missingIn100)
    g_mis100 = goeaobj.run_study(missingIn100)
    g_mis100Sig = [r for r in g_mis100 if r.p_fdr_bh < 0.05]
    plot_results(missing100png, g_mis100Sig)

    print('present in 90% of species:\n')
    g_pres100 = goeaobj.run_study(presentIn100)
    g_pres100Sig = [r for r in g_pres100 if r.p_fdr_bh < 0.05]
    plot_results(present100png, g_pres100Sig)
    
    print('fraction from parent:\n')
    g_mis100Priv = goeaobjPrivate.run_study(missingIn100)
    g_mis100SigPriv = [r for r in g_mis100Priv if r.p_fdr_bh < 0.05]
    plot_results(missing100pngPrivate, g_mis100SigPriv)
    




for order in orderSpecies.keys():
    if not os.path.exists('output/order/' + order):
        os.mkdir('output/order/' + order)

    speciesCount = len(orderSpecies[order])
    genePercentage = {}
    for i in genesWithAnn:
        genePercentage[i] = 0
    for species in orderSpecies[order]:
        for gene in speciesList[species].HQGenes:
            genePercentage[gene] += 1
    
    percentFile = 'output/order/' + order + '/percentPresent.txt'
    missingFile = 'output/order/' + order + '/missingInGroup.txt'
    presentFile = 'output/order/' + order + '/presentInGroup.txt'
    present100File = 'output/order/' + order + '/presentIn100Group.txt'
    present95File = 'output/order/' + order + '/presentIn95Group.txt'
    present90File = 'output/order/' + order + '/presentIn90Group.txt'
    missing95File = 'output/order/' + order + '/missingIn95Group.txt'
    missing90File = 'output/order/' + order + '/missingIn90Group.txt'
    
    
    
    with open(percentFile, 'w') as percentOut, \
        open(missingFile, 'w') as missingOut, \
        open(presentFile, 'w') as presentOut, \
        open(present100File, 'w') as present100Out, \
        open(present95File, 'w') as present95Out, \
        open(present90File, 'w') as present90Out, \
        open(missing95File, 'w') as missing95Out, \
        open(missing90File, 'w') as missing90Out:
            
        for i in genePercentage.keys():
            percentage = genePercentage[i]/speciesCount
            percentOut.write(i + '\t' + geneNames[i] + '\t' + str(percentage) + '\n')
            if percentage > 0:
                presentOut.write(i + '\t' + geneNames[i] + '\n')
            if percentage == 1:
                present100Out.write(i + '\t' + geneNames[i] + '\n')
            if percentage == 0:
                missingOut.write(i + '\t' + geneNames[i] + '\n')
            if percentage >= .9:
                present90Out.write(i + '\t' + geneNames[i] + '\n')
            if percentage >= .95:
                present95Out.write(i + '\t' + geneNames[i] + '\n')
            if percentage <= .1:
                missing90Out.write(i + '\t' + geneNames[i] + '\n')
            if percentage <= .05:
                missing95Out.write(i + '\t' + geneNames[i] + '\n')
                
for genus in genusSpecies.keys():
    if not os.path.exists('output/genus/' + genus):
        os.mkdir('output/genus/' + genus)

    speciesCount = len(genusSpecies[genus])
    genePercentage = {}
    for i in genesWithAnn:
        genePercentage[i] = 0
    for species in genusSpecies[genus]:
        for gene in speciesList[species].HQGenes:
            genePercentage[gene] += 1
    
    percentFile = 'output/genus/' + genus + '/percentPresent.txt'
    missingFile = 'output/genus/' + genus + '/missingInGroup.txt'
    presentFile = 'output/genus/' + genus + '/presentInGroup.txt'
    present100File = 'output/genus/' + genus + '/presentIn100Group.txt'
    present95File = 'output/genus/' + genus + '/presentIn95Group.txt'
    present90File = 'output/genus/' + genus + '/presentIn90Group.txt'
    missing95File = 'output/genus/' + genus + '/missingIn95Group.txt'
    missing90File = 'output/genus/' + genus + '/missingIn90Group.txt'
    
    
    
    with open(percentFile, 'w') as percentOut, \
        open(missingFile, 'w') as missingOut, \
        open(presentFile, 'w') as presentOut, \
        open(present100File, 'w') as present100Out, \
        open(present95File, 'w') as present95Out, \
        open(present90File, 'w') as present90Out, \
        open(missing95File, 'w') as missing95Out, \
        open(missing90File, 'w') as missing90Out:
            
        for i in genePercentage.keys():
            percentage = genePercentage[i]/speciesCount
            percentOut.write(i + '\t' + geneNames[i] + '\t' + str(percentage) + '\n')
            if percentage > 0:
                presentOut.write(i + '\t' + geneNames[i] + '\n')
            if percentage == 1:
                present100Out.write(i + '\t' + geneNames[i] + '\n')
            if percentage == 0:
                missingOut.write(i + '\t' + geneNames[i] + '\n')
            if percentage >= .9:
                present90Out.write(i + '\t' + geneNames[i] + '\n')
            if percentage >= .95:
                present95Out.write(i + '\t' + geneNames[i] + '\n')
            if percentage <= .1:
                missing90Out.write(i + '\t' + geneNames[i] + '\n')
            if percentage <= .05:
                missing95Out.write(i + '\t' + geneNames[i] + '\n')
                
for family in familySpecies.keys():
    if not os.path.exists('output/family/' + family):
        os.mkdir('output/family/' + family)

    speciesCount = len(familySpecies[family])
    genePercentage = {}
    for i in genesWithAnn:
        genePercentage[i] = 0
    for species in familySpecies[family]:
        for gene in speciesList[species].HQGenes:
            genePercentage[gene] += 1
    
    percentFile = 'output/family/' + family + '/percentPresent.txt'
    missingFile = 'output/family/' + family + '/missingInGroup.txt'
    presentFile = 'output/family/' + family + '/presentInGroup.txt'
    present100File = 'output/family/' + family + '/presentIn100Group.txt'
    present95File = 'output/family/' + family + '/presentIn95Group.txt'
    present90File = 'output/family/' + family + '/presentIn90Group.txt'
    missing95File = 'output/family/' + family + '/missingIn95Group.txt'
    missing90File = 'output/family/' + family + '/missingIn90Group.txt'
    
    
    
    with open(percentFile, 'w') as percentOut, \
        open(missingFile, 'w') as missingOut, \
        open(presentFile, 'w') as presentOut, \
        open(present100File, 'w') as present100Out, \
        open(present95File, 'w') as present95Out, \
        open(present90File, 'w') as present90Out, \
        open(missing95File, 'w') as missing95Out, \
        open(missing90File, 'w') as missing90Out:
            
        for i in genePercentage.keys():
            percentage = genePercentage[i]/speciesCount
            percentOut.write(i + '\t' + geneNames[i] + '\t' + str(percentage) + '\n')
            if percentage > 0:
                presentOut.write(i + '\t' + geneNames[i] + '\n')
            if percentage == 1:
                present100Out.write(i + '\t' + geneNames[i] + '\n')
            if percentage == 0:
                missingOut.write(i + '\t' + geneNames[i] + '\n')
            if percentage >= .9:
                present90Out.write(i + '\t' + geneNames[i] + '\n')
            if percentage >= .95:
                present95Out.write(i + '\t' + geneNames[i] + '\n')
            if percentage <= .1:
                missing90Out.write(i + '\t' + geneNames[i] + '\n')
            if percentage <= .05:
                missing95Out.write(i + '\t' + geneNames[i] + '\n')

'''
                
            
#find the genes that always fail

    
            
    
            
        
            
        
        