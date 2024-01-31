import pysam
from pysam import VariantFile, FastaFile
import re
import csv
import argparse
import gzip



##Example usage:
# python parseVCF.py -i ../SNPeff_output/LymphomaForTesting/*.vcf /seq/vgb/cancer_r01/cbioportal/angiosarcoma/data_mutations_mskcc.txt /seq/vgb/cancer_r01/icgc/*.tsv.gz /seq/vgb/cancer_r01/tcga_data/maf_files/*/*.somatic.maf.gz -r /seq/vgb/swofford/ref/EnsemblHumanDogOrthologs.txt -p /seq/vgb/swofford/ref/pathways/C2CuratedGeneSetsGeneNames.gmt.txt -o allTest -hf /seq/vgb/swofford/ref/hg38.fa -df /seq/vgb/swofford/ref/Canis_lupus_familiaris_assembly3.fasta -hg19 /seq/vgb/swofford/ref/hg19plusChrM.fa -bp -tri -m -gc

'''Need to modify vcf parsing so as not to count variants more than once in the mutation signatures stuff. Check other
filetypes, but I think it's just the VCFs that has this issue. -done for vcfs, need to add for everything, for saftey'''


parser = argparse.ArgumentParser()
parser.add_argument('-i', '--inFile', nargs='+')
parser.add_argument('-sp', '--sparsityMin', default = 0)
parser.add_argument('-o', '--outFile')
parser.add_argument('-s', '--species')
parser.add_argument('-r', '--refFile')
parser.add_argument('-p', '--pathways')
parser.add_argument('-bp', '--brokenPathways', action='store_true', default=False)
parser.add_argument('-g', '--brokenGene', action='store_true', default=False)
parser.add_argument('-c', '--cores', default = 1)
parser.add_argument('-maxVar', '--maximumVariantsSample')
parser.add_argument('-m', '--mutationCount', action='store_true', default=False)
parser.add_argument('-hf', '--humanFasta')
parser.add_argument('-hg19', '--hg19Fasta')
parser.add_argument('-df', '--dogFasta')
parser.add_argument('-tri', '--triNucCon', action = 'store_true', default=False)
parser.add_argument('-gc', '--geneCount', action = 'store_true' , default=False)
parser.add_argument('-pc', '--pathwayCount', action = 'store_true', default=False)
parser.add_argument('-gr', '--granularity', type = int, default = 1)

args = parser.parse_args()



'''Samples are stored as an instance of the sample class. Each sample has a label which
will be used for learning at a later stage, a dictionary of variants, a species (currently poorly implemented)
a list of broken genes to speed lookup, a count of broken genes, a list of broken pathways, a count of variants,
and a list of features to be written to an output file. feature list is added after parsing based on options.
'''
 
'''added set{} test on 8/20 but have not tested. if breaks, remove or start here'''
class sample:
    def __init__(self, variants, name, species, label = 'NA'):
        pathBrokenInSample = {}
        genesBrokenInSample = {}
        ####
        variantSet = set()
        codingVariants = []
        
        for variant in variants:
            ###
            if variant.variantID not in variantSet and variant.geneID != 'NA':
            ###
                variantSet.add(variant.variantID)
                codingVariants.append(variant)
                if variant.breaksGene:
                    try:
                        genesBrokenInSample[variant.geneID] += 1                    
                    except KeyError:
                        genesBrokenInSample[variant.geneID] = 1
                    if variant.geneID != 'NA':
                        for pathway in reference[variant.geneID].pathways:
                            try:
                                pathBrokenInSample[pathway[0]] += 1
                            except KeyError:
                                pathBrokenInSample[pathway[0]] = 1
        
        self.species = species
        self.name = name
        self.variants = codingVariants                        
        self.brokenGenes = genesBrokenInSample.keys()
        self.brokenGeneCount = genesBrokenInSample
        self.brokenPathways = pathBrokenInSample.keys()
        self.brokenPathwayCount = pathBrokenInSample
        self.variantCount = len(variantSet)
        self.toWrite = []
        #Just for testing purposes, found rando vcfs in one of Ginger's folders, 
        #ran snpeff on them and am calling them "Lymphoma" (they all have 'PONfilt' in their filenames)
        '''
        if label == 'NA' and 'PONfilt' in name:
            self.label = 'Lymphoma'
        elif label == 'NA':
            self.label = 'HSA'
        else:
        '''
        self.label = label
        labels.add(label)


''' Will need to modify the below to better work wth the icgc data. Probably should get rid of the "fields" input from vcf and just pass values one at a time, so can be consistent coming from any source.
so, move the converting from fields array into the vcf portion, rather than doing here at the variant creation level. Need comprehensive list of types of info '''

class variant:
    
    def __init__(self, Effect_Impact = 'NA', Functional_Class = 'NA', Codon_change = 'NA', Amino_Acid_Change = 'NA', 
                 Amino_Acid_length = 'NA', Gene_Name = 'NA', Transcript_Biotype = 'NA', Gene_Coding = 'NA', 
                 Transcript_ID = 'NA', Exon_Rank = 'NA', Genotype = 'NA', Effect = 'NA', chrom = 'NA', pos = 'NA', ref = 'NA', FASTA='NA', changeTri='NA'):
         
        geneID = getHumanGeneID(Gene_Name)
        
        self.chrom = chrom
        self.pos = pos
        self.Effect = Effect
        self.Effect_Impact = Effect_Impact
        self.Functional_Class = Functional_Class
        self.Codon_change = Codon_change
        self.Amino_Acid_Change = Amino_Acid_Change
        self.Amino_Acid_length = Amino_Acid_length
        self.Gene_Name = Gene_Name
        self.Transcript_BioType = Transcript_Biotype
        self.Gene_Coding = Gene_Coding
        self.Transcript_ID = Transcript_ID
        self.Exon_Rank = Exon_Rank
        self.Genotype = Genotype.upper()
        self.ref = ref.upper()
        self.geneID = geneID
        self.dogGeneID = getDogGeneID(Gene_Name)
        self.variantID = chrom+str(pos)+Genotype.upper()
        
        if ('missense' in Effect 
            or 'frameshift_variant' in Effect 
            or 'disruptive_inframe_deletion' in Effect 
            or 'stop_gained' in Effect 
            or 'stop_lost' in Effect
            or 'MISSENSE' in Effect
            or 'NONSENSE' in Effect
            or Effect in damagingEffects):
                self.breaksGene = True
                try:
                    brokenGenes[geneID] += 1
                except KeyError:
                    brokenGenes[geneID] = 1
        else:
            self.breaksGene = False
                    
        if len(ref) == 1 and len(Genotype) == 1 and ref != '-' and Genotype != '-':
            
            refTrinucleotide = FASTA.fetch(chrom,pos-2,pos+1).upper() 
            
            ref = ref.upper()
            
            if ref == 'A' or ref == 'G':
                refTrinucleotide = reverse(complement(refTrinucleotide))
                ref = complement(ref)
                Genotype = complement(Genotype)
                
            changeTri = refTrinucleotide[:1] + '[' + ref + '>' + Genotype + ']' + refTrinucleotide[2:]
            
            if refTrinucleotide[1:2] != ref:
                print('AAAAAHHHHH WHAT DO YOU DO TO THE BRAIN')
                print(refTrinucleotide)
                print(ref)
                print(Genotype)################################################            
                
                     
        self.trinucleotide = changeTri
       
        
        


#SNPeff output has the following in the EFF info section:
#Type=String,Description="Predicted effects for this variant.Format: 'Effect ( Effect_Impact | Functional_Class | Codon_Change | Amino_Acid_Change| Amino_Acid_length | Gene_Name | Transcript_BioType | Gene_Coding | Transcript_ID | Exon_Rank  | Genotype [ | ERRORS | WARNINGS ] )'">

#icgc output guide chrom = chromosome, pos = chromosome_start, Effect = consequence_type, ref = reference_genome_allele, Genotype = mutated_to_allele

''' Class defining a single gene ortholog relationship between dog and human, to facilitate bi-directional lookups later'''
class ortholog:
    def __init__(self, dogGeneName, humanGeneName, dogGeneID, geneID, dogChr, humanChr, dogStart, dogEnd, humanStart, humanEnd, conservationScore, orthologConfidence, pathways):
        self.dogGeneName = dogGeneName
        self.humanGeneName = humanGeneName
        self.dogGeneID = dogGeneID
        self.geneID = geneID
        self.dogChr = dogChr
        self.humanChr = humanChr
        self.dogStart = dogStart
        self.dogEnd = dogEnd
        self.humanStart = humanStart
        self.humanEnd = humanEnd
        self.conservationScore = conservationScore
        self.orthologConfidence = orthologConfidence
        self.pathways = pathways


'''builds a reference table lookup between human/dog orthologous genes, parsed from an Ensembl table. File must contain the below columns, extra columns are fine and will be ignored.'''
#include pathway stuff here? 

def buildRef(refFile, pathways):
    with open(refFile, 'r') as refFile, open(pathways, 'r') as pathRef:
        ref = {}
        pathwayList = {}
        genePath = {}
        for line in pathRef:
            genes = line.split('\t')
            pathGeneList = []
            pathName = genes[0]
            print(pathName)
            source = genes[1]
            print(source)
            
            for gene in genes[2:]:
                try:
                    genePath[gene].append([pathName])
                except KeyError:
                    genePath[gene] = [pathName]

        print(str(len(genePath)) + ' genes with pathways found.')
            
        
        refReader =  csv.DictReader(refFile, delimiter='\t')
        for line in refReader:
            if line['Dog homology type'] == 'ortholog_one2one' and line['Dog gene name'] != '':
                
                dogGeneName = line['Dog gene name']
                humanGeneName = line['Gene name']
                dogGeneID = line['Dog gene stable ID']
                geneID = line['Gene stable ID']
                dogChr = line['Dog chromosome/scaffold name']
                humanChr = line['Chromosome/scaffold name']
                dogStart = line['Dog chromosome/scaffold start (bp)']
                dogEnd = line['Dog chromosome/scaffold end (bp)']
                humanStart = line['Gene start (bp)']
                humanEnd = line['Gene end (bp)']
                conservationScore = line['Dog Gene-order conservation score']
                orthologConfidence = line['Dog orthology confidence [0 low, 1 high]']
                
                try:
                    pathwaysEffected = genePath[humanGeneName]
                except KeyError:
                    pathwaysEffected = ''
                        
                
                ref[geneID] = ortholog(dogGeneName, humanGeneName, dogGeneID, geneID, dogChr, humanChr, dogStart, dogEnd, humanStart, humanEnd, conservationScore, orthologConfidence, pathwaysEffected) 
                #following: might be useful to have an easy lookup from either directiion (lookup from dog id, vs human id) possible solution is just to store a version for each - cludgy, and implemented!
                ref[dogGeneID] = ortholog(dogGeneName, humanGeneName, dogGeneID, geneID, dogChr, humanChr, dogStart, dogEnd, humanStart, humanEnd, conservationScore, orthologConfidence, pathwaysEffected)
                ref[dogGeneName] = ortholog(dogGeneName, humanGeneName, dogGeneID, geneID, dogChr, humanChr, dogStart, dogEnd, humanStart, humanEnd, conservationScore, orthologConfidence, pathwaysEffected)

        #print('the length of ref is ' + str(len(ref)))
        print('Finished building Reference. ' + str(len(ref)) + ' orthologs recorded.')
        return ref

'''get the corresponding dog stableGeneID from human stable geneID'''
def getDogGeneID(geneID):
    try:
        return reference[geneID].dogGeneID
    except (KeyError):
        return 'NA'

'''get the corresponding Human stable gene ID from the stable dog geneID'''
def getHumanGeneID(geneID):
    try:
        return reference[geneID].geneID
    except (KeyError):
        return 'NA'

'''get ortholog object for geneID'''
def getOrtholog(geneID):
    try:
        return reference[geneID]
    except (KeyError):
        return 'NA'

'''find list of genes that are "broken" in any of the samples. May need to modify what snpeff effects truly break genes'''            
def getBrokenGenes(samples):
    genes = set()
    for sample in samples:
        for gene in sample.brokenGenes:
            genes.add(gene)
        
        '''
        for variant in sample.variants:
            if ('missense' in variant.Effect 
            or 'frameshift_variant' in variant.Effect 
            or 'disruptive_inframe_deletion' in variant.Effect 
            or 'stop_gained' in variant.Effect 
            or 'stop_lost' in variant.Effect):
                genes.add(variant.geneID)
        '''
    return(genes)


'''queries whether or not a given gene is broken in an individual sample, only considers genes with human/dog ortholog pairs'''
def isGeneBroken(sample, geneID):
    
    
    if geneID in sample.brokenGenes:
        return True
    return False

    
def isPathwayBroken(sample, pathway):
    if pathway in sample.brokenPathways:
        return True
    return False      

def complement(seq):
    '''Returns the complement of a DNA sequence passed to it. Accepts a
    string as input, and returns a string. Valid bases for input are
    "ATCGWN"'''
    compSeq = ''
    codes = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'W': 'W', 'N':'N'}
    for letter in seq:
        compSeq += codes[letter]
    return(compSeq)


def reverse(seq):
    '''Returns the reverse of a string passed to it. Input is a string,
    output is a string.'''
    return(seq[::-1])



brokenGenes = {}

try:
    reference = buildRef(args.refFile, args.pathways)
except:
    print('Failed to load reference files, dying. Because I can\'t even.')

samples = []
labels = set()
featureList = []
#genes = set()
if args.maximumVariantsSample:
    maxVar = int(args.maximumVariantsSample)
else:
    maxVar = 1000000000

if args.granularity == 1:    
    studyClass = {
        'ALL-US':'ALL',
        'B-cell acute lymphoblastic leukemia, non-hypodiploid':'ALL',
        'B-cell acute lymphoblastic leukemia, hypodiploid':'ALL',
        'T-cell acute lymphoblastic leukemia':'ALL',
        'AML-US':'AML',
        'LAML-CN':'Remove',
        'LAML-KR':'AML',
        'Acute myeloid leukemias':'AML',
        'LAML':'AML',
        'ACC':'ACC',
        'Adrenocortical carcinoma':'ACC',
        'BLCA-CN':'BLCA',
        'BLCA':'BLCA',
        'BOCA-FR':'Remove',
        'BOCA-UK':'Remove',
        'Ewing\'s sarcoma':'Ewing\'s',
        'Osteosarcoma':'OSA',
        'BPLL-FR':'BPLL',
        'BRCA-EU':'BRCA',
        'BRCA-FR':'BRCA',
        'BRCA-KR':'BRCA',
        'BRCA-UK':'BRCA',
        'BRCA':'BRCA',
        'BTCA-JP':'BTCA',
        'BTCA-SG':'BTCA',
        'CHOL':'BTCA',
        'CESC':'CESC',
        'CLLE-ES':'CLLE',
        'CMDI-UK':'Remove',
        'COCA-CN':'CORE',
        'COAD':'CORE',
        'READ':'CORE',
        'Embryonal tumor with multilayered rosettes':'ETMR',
        'Ependymoma supratentorial':'Ependymoma',
        'Ependymoma infratentorial':'Ependymoma',
        'ESAD-UK':'ESCA',
        'ESCA-CN':'ESCA',
        'ESCA':'ESCA',
        'GACA-CN':'GACA',
        'GACA-JP':'GACA',
        'STAD':'GACA',
        'GBM':'Glioma',
        'High-grade glioma K27Mmut':'Glioma',
        'High-grade glioma K27Mwt':'Glioma',
        'LGG':'Glioma',
        'HNSC':'HNSC',
        'Hepatoblastoma':'Hepatoblastoma',
        'KIRP':'KIRP',
        'KICH':'KICH',
        'KIRC':'KIRC',
        'RECA-CN':'KIRC',
        'RECA-EU':'KIRC',
        'Wilms\' tumors':'WT',
        'WT-US':'WT',
        'LICA-CN':'LICA',
        'LICA-FR':'LICA',
        'LIHC':'LICA',
        'LINC-JP':'LICA',
        'LIRI-JP':'LICA',
        'LIHM-FR':'Remove',
        'LIAD-FR':'Remove',
        'LMS-FR':'LMS',
        'LUAD':'LUAD',
        'LUSC':'LUSC',
        'LUSC-CN':'LUSC',
        'LUSC-KR':'Remove',
        'DLBC':'LSA',
        'MALY-DE':'LSA',
        'Burkitt\'s lymphoma':'LSA',
        'NKTL-SG':'LSA',
        'Medulloblastoma WNT':'PEME',
        'Medulloblastoma SHH':'PEME',
        'Medulloblastoma Group3':'PEME',
        'Medulloblastoma Group4':'PEME',
        'PEME-CA':'PEME',
        'PBCA-DE':'Remove',
        'SKCM':'MELA',
        'MELA-AU':'MELA',
        'UVM':'MELA',
        'MESO':'MESO',
        'NACA-CN':'NACA',
        'NBL-US':'NBL',
        'Neuroblastoma':'NBL',
        'ORCA-IN':'ORCA',
        'OV':'OV',
        'OV-AU':'OV',
        'PACA-AU':'PACA',
        'PACA-CA':'PACA',
        'PAAD':'PACA',
        'PAEN-AU':'PAEN',
        'PAEN-IT':'PAEN',
        'PBCA-US':'Remove',
        'PCPG':'PCPG',
        'PBCA-DE':'Remove',
        'Pilocytic astrocytoma':'PAST',
        'EOPC-DE':'PRAD',
        'PRAD-CA':'PRAD',
        'PRAD-CN':'PRAD',
        'PRAD-FR':'PRAD',
        'PRAD-UK':'PRAD',
        'PRAD':'PRAD',
        'Retinoblastoma':'RB',
        'RT-US':'RT',
        'Atypial teratoid/rhabdoid tumor':'RT',
        'Rhabdomyosarcoma':'RMS',
        'SARC':'Remove',
        'SKCA-BR':'SKAD',
        'TGCT':'TGCT',
        'THYM':'THYM',
        'THCA-CN':'THCA',
        'THCA-SA':'THCA',
        'THCA':'THCA',
        'UTCA-FR':'UCS',
        'UCS':'UCS',
        'UCEC':'UCEC',
        'Soft Tissue Sarcoma':'AS'
        }

elif args.granularity == 2:

    studyClass= {
        'ALL-US':'Remove',
        'B-cell acute lymphoblastic leukemia, non-hypodiploid':'BALL',
        'B-cell acute lymphoblastic leukemia, hypodiploid':'BALL',
        'T-cell acute lymphoblastic leukemia':'TALL',
        'AML-US':'AML',
        'LAML-CN':'Remove',
        'LAML-KR':'AML',
        'Acute myeloid leukemias':'AML',
        'LAML':'AML',
        'ACC':'ACC',
        'Adrenocortical carcinoma':'ACC',
        'BLCA-CN':'BLCA',
        'BLCA':'BLCA',
        'BOCA-FR':'Remove',
        'BOCA-UK':'Remove',
        'Ewing\'s sarcoma':'Ewing\'s',
        'Osteosarcoma':'OSA',
        'BPLL-FR':'BPLL',
        'BRCA-EU':'BRCA',
        'BRCA-FR':'BRCA',
        'BRCA-KR':'BRCA',
        'BRCA-UK':'BRCA',
        'BRCA':'BRCA',
        'BTCA-JP':'Remove',
        'BTCA-SG':'CHOL',
        'CHOL':'CHOL',
        'CESC':'CESC',
        'CLLE-ES':'CLLE',
        'CMDI-UK':'Remove',
        'COCA-CN':'Remove',
        'COAD':'COAD',
        'READ':'READ',
        'Embryonal tumor with multilayered rosettes':'ETMR',
        'Ependymoma supratentorial':'Ependymoma supratentorial',
        'Ependymoma infratentorial':'Ependymoma infratentorial',
        'ESAD-UK':'ESAD',
        'ESCA-CN':'ESCA',
        'ESCA':'Remove',
        'GACA-CN':'',
        'GACA-JP':'Remove',
        'STAD':'STAD',
        'GBM':'GBM',
        'High-grade glioma K27Mmut':'HGG',
        'High-grade glioma K27Mwt':'HGG',
        'LGG':'LGG',
        'HNSC':'HNSC',
        'Hepatoblastoma':'Hepatoblastoma',
        'KIRP':'KIRP',
        'KICH':'KICH',
        'KIRC':'KIRC',
        'RECA-CN':'KIRC',
        'RECA-EU':'KIRC',
        'Wilms\' tumors':'WT',
        'WT-US':'WT',
        'LICA-CN':'VLICA',
        'LICA-FR':'LICA',
        'LIHC':'LICA',
        'LINC-JP':'VLICA',
        'LIRI-JP':'VLICA',
        'LIHM-FR':'Remove',
        'LIAD-FR':'Remove',
        'LMS-FR':'LMS',
        'LUAD':'LUAD',
        'LUSC':'LUSC',
        'LUSC-CN':'LUSC',
        'LUSC-KR':'Remove',
        'DLBC':'DLBCL',
        'MALY-DE':'Remove',
        'Burkitt\'s lymphoma':'BULY',
        'NKTL-SG':'NKTL',
        'Medulloblastoma WNT':'PEME',
        'Medulloblastoma SHH':'PEME',
        'Medulloblastoma Group3':'PEME',
        'Medulloblastoma Group4':'PEME',
        'PEME-CA':'PEME',
        'PBCA-DE':'Remove',
        'SKCM':'MELA',
        'MELA-AU':'MELA',
        'UVM':'UVM',
        'MESO':'MESO',
        'NACA-CN':'NACA',
        'NBL-US':'NBL',
        'Neuroblastoma':'NBL',
        'ORCA-IN':'ORCA',
        'OV':'OV',
        'OV-AU':'OV',
        'PACA-AU':'PACA',
        'PACA-CA':'PACA',
        'PAAD':'PACA',
        'PAEN-AU':'PAEN',
        'PAEN-IT':'PAEN',
        'PBCA-US':'Remove',
        'PCPG':'PCPG',
        'PBCA-DE':'Remove',
        'Pilocytic astrocytoma':'PAST',
        'EOPC-DE':'PRAD',
        'PRAD-CA':'PRAD',
        'PRAD-CN':'PRAD',
        'PRAD-FR':'PRAD',
        'PRAD-UK':'PRAD',
        'PRAD':'PRAD',
        'Retinoblastoma':'RB',
        'RT-US':'RT',
        'Atypial teratoid/rhabdoid tumor':'ATRT',
        'Rhabdomyosarcoma':'RMS',
        'SARC':'Remove',
        'SKCA-BR':'SKAD',
        'TGCT':'TGCT',
        'THYM':'THYM',
        'THCA-CN':'THCA',
        'THCA-SA':'THCA',
        'THCA':'THCA',
        'UTCA-FR':'UCS',
        'UCS':'UCS',
        'UCEC':'UCEC',
        'Soft Tissue Sarcoma':'AS'
        }

elif int(args.granularity) > 2 or int(args.granularity) < 1:
    print('Level of granularity selected is out of bounds, I die now.')
    quit()
        
ASKeep = ['Angio-ASCProject_03IpIBIR-Tumor-SM-DAE1F',
    'Angio-ASCProject_09iyC6Cp-Tumor-SM-DACME',
    'Angio-ASCProject_3NflfGHo-Tumor-SM-DADBW',
    'Angio-ASCProject_5DFyF0HA-Tumor-SM-DAD34',
    'Angio-ASCProject_5NsKslTb-Tumor-SM-DACSO',
    'Angio-ASCProject_5Ph4hetv-Tumor-SM-DAE5T',
    'Angio-ASCProject_bJi9iyCb-Tumor-SM-DAE73',
    'Angio-ASCProject_blh1uyhK-Tumor-SM-DADU5',
    'Angio-ASCProject_bRSGSlCG-Tumor-SM-DACJV',
    'Angio-ASCProject_EAuOu4uD-Tumor-SM-DADQD',
    'Angio-ASCProject_EdCVC0um-Tumor-SM-DADA1',
    'Angio-ASCProject_eVIAsEIy-Tumor-SM-DADLC',
    'Angio-ASCProject_GXTMTxU7-Tumor-SM-DACLR',
    'Angio-ASCProject_JkivIBSV-Tumor-SM-DACIM',
    'Angio-ASCProject_kJFaIXcN-Tumor-SM-DAE4K',
    'Angio-ASCProject_KxFGsofW-Tumor-SM-DAD4D',
    'Angio-ASCProject_loC8UQUk-Tumor-SM-DACX3',
    'Angio-ASCProject_lYIMuWu5-Tumor-SM-DACJ9',
    'Angio-ASCProject_mpUNibFx-Tumor-SM-DAD3Q',
    'Angio-ASCProject_NdUxUwCM-Tumor-SM-DAE57',
    'Angio-ASCProject_nvfdfZCz-Tumor-SM-DADGX',
    'Angio-ASCProject_oQi7i3U5-Tumor-SM-DADHK',
    'Angio-ASCProject_QYsAsrTO-Tumor-SM-DAE3X',
    'Angio-ASCProject_XbieiohK-Tumor-SM-DADIT',
    'Angio-ASCProject_xdigi2t7-Tumor-SM-DADLY',
    'Angio-ASCProject_XjHAU9uR-Tumor-SM-DADEF',
    'Angio-ASCProject_Y4SjSnIg-Tumor-SM-DADSQ',
    'Angio-ASCProject_YLCRINC8-Tumor-SM-DACNN',
    'Angio-ASCProject_YzTlTytj-Tumor-SM-DAD69',
    'Angio-ASCProject_YzTxcouV-Tumor-SM-DAD9E',
    'Angio-ASCProject_zQtVt5tr-Tumor-SM-E76O4',
    'Angio-ASCProject_7bt2t8IJ-Tumor-SM-DAD18',
    'Angio-ASCProject_aPskieIL-Tumor-SM-DAD5M',
    'Angio-ASCProject_dyhLT8sG-Tumor-SM-DACWG',
    'Angio-ASCProject_RvtOtjtj-Tumor-SM-DACOW',
    'Angio-ASCProject_ZXi6FXFo-Tumor-SM-DADXW'
    ]

damagingEffects = ['chromosome','CHROMOSOME_LARGE DELETION',
                   'duplication','CHROMOSOME_LARGE_DUPLICATION',
                   'inversion','CHROMOSOME_LARGE_INVERSION',
                   'exon_loss_variant','EXON_DELETED',
                   'exon_loss_variant','EXON_DELETED_PARTIAL',
                   'duplication','EXON_DUPLICATION',
                   'duplication','EXON_DUPLICATION_PARTIAL',
                   'inversion','EXON_INVERSION',
                   'inversion','EXON_INVERSION_PARTIAL',
                   'frameshift_variant','FRAME_SHIFT',
                   'feature_ablation','GENE_DELETED',
                   'gene_fusion','GENE_FUSION',
                   'gene_fusion','GENE_FUSION_HALF',
                   'bidirectional_gene_fusion','GENE_FUSION_REVERESE',
                   'rearranged_at_DNA_level','GENE_REARRANGEMENT',
                   'protein_protein_contact','PROTEIN_PROTEIN_INTERACTION_LOCUS',
                   'structural_interaction_variant','PROTEIN_STRUCTURAL_INTERACTION_LOCUS',
                   'rare_amino_acid_variant','RARE_AMINO_ACID',
                   'splice_acceptor_variant','SPLICE_SITE_ACCEPTOR',
                   'splice_donor_variant','SPLICE_SITE_DONOR',
                   'stop_lost','STOP_LOST',
                   'start_lost','START_LOST',
                   'stop_gained','STOP_GAINED',
                   'feature_ablation','TRANSCRIPT_DELETED',
                   'inframe_insertion','CODON_INSERTION',
                   'disruptive_inframe_insertion','CODON_CHANGE_PLUS CODON_INSERTION',
                   'inframe_deletion','CODON_DELETION',
                   'disruptive_inframe_deletion','CODON_CHANGE_PLUS CODON_DELETION',
                   'duplication','GENE_DUPLICATION',
                   'missense_variant','NON_SYNONYMOUS_CODING',
                   'splice_region_variant','SPLICE_SITE_BRANCH_U12',
                   '3_prime_UTR_truncation + exon_loss','UTR_3_DELETED',
                   '5_prime_UTR_truncation + exon_loss_variant','UTR_5_DELETED'
                   ]
                   
    

    
# Go through individual files and collect variant info for each sample

with pysam.Fastafile(args.humanFasta) as humFast, pysam.Fastafile(args.dogFasta) as dogFast, pysam.Fastafile(args.hg19Fasta) as hg19Fast:
        
    #currently requiring fasta reference file, can fix later if necessary. Will need separate dog and human fasta references.
    #with pysam.Fastafile(args.fasta) as fastaRef:        
        
    for fileName in args.inFile:
        print(fileName)
        
        if 'Heather' in fileName:
            with open(fileName, 'rt') as inFile:
                sampleVariants = {}
                label = 'dogOSA'
                reader = csv.DictReader(inFile, delimiter = '\t')
                for line in reader:
                    name = line['Dog ID']
                    chrom = 'chr' + line['Chromosome']
                    pos = int(line['Position'])
                    ref = line['Reference Allele']
                    Genotype = line['Alternate Allele']
                    Gene_Name = line['Gene Name']
                    Gene_ID = line['Ensembl Gene ID']
                    Transcript_ID = line['Ensembl Transcript ID']
                    Effect = line['Effect']
                    
                    newVariant = variant(Gene_Name = Gene_Name, Transcript_ID = Transcript_ID, 
                                         Genotype = Genotype, Effect = Effect, chrom = chrom, pos = pos, ref = ref, FASTA = dogFast)
                     
                    try:
                        sampleVariants[name].append(newVariant)
                    except KeyError:
                        sampleVariants[name] = [newVariant]

                for name in sampleVariants:
                    if len(sampleVariants[name]) < maxVar:
                        samples.append(sample(sampleVariants[name], name, 'Dog', label))
        '''
        if '.final.txt' in fileName:
            with open(fileName, 'rt') as inFile:
                names = ['chrom', 'pos', 'Gene_Name', 'Effect', 'Amino_Acid_Change', 'unk', 'ref', 'genotype', 'Effect_Impact', 'Transcript_ID']
                reader = csv.DictReader(inFile, delimiter = '\t', fieldnames=names)
                sampleVariants = {}
                
                for line in reader:
                    name = line['Tumor_Sam']      
                    Amino_Acid_Change = line['Amino_acids']
                    Exon_Rank = line['Exon_Number']
                    Gene_Name = line['Gene']
                    Amino_Acid_Change = line['HGVSp']
                    Transcript_ID = line['Feature']
                    Genotype = line['Tumor_Seq_Allele2']
                    Effect = line['One_Consequence']
                    chrom = line['Chromosome']
                    pos = int(line['Start_Position'])
                    ref = line['Reference_Allele']
                    
                    newVariant = variant(Amino_Acid_Change = Amino_Acid_Change, 
                                         Gene_Name = Gene_Name, Exon_Rank = Exon_Rank, Transcript_ID = Transcript_ID, 
                                         Genotype = Genotype, Effect = Effect, chrom = chrom, pos = pos, ref = ref, FASTA = humFast)
                     
                    try:
                        sampleVariants[name].append(newVariant)
                    except KeyError:
                        sampleVariants[name] = [newVariant]

                for name in sampleVariants:
                    if len(sampleVariants[name]) < maxVar:
                        samples.append(sample(sampleVariants[name], name, 'Human', label))
        '''            
                
        
        if 'data_mutations' in fileName:
            basePath = ''
            fullPath = fileName.split('/')
            for i in range(len(fullPath) - 1):
                basePath = basePath + fullPath[i] + '/'
            phenotypeFile = basePath + 'data_clinical_sample.txt'
            

            with open(fileName, 'rt') as maf_in, open(phenotypeFile, 'rt') as pheno_in:
                sampleLabel = {}
                sampleVariants = {}
                
                # skip header, if there is one, and get column names
                for line in pheno_in:
                    if line[:1] != '#':
                        line = line.strip()
                        names = line.split('\t')
                        break
                    
                phenoReader = csv.DictReader(pheno_in, delimiter = '\t', fieldnames = names)
                for line in phenoReader:
                    #pediatric pan-cancer sample names end in '-R' when sample is from a relapse, rather than primary
                    #We're going to ignore the relapse samples
                    if studyClass[line['CANCER_TYPE']] != 'AS':
                        if '-R' not in line['SAMPLE_ID']:
                            sampleID = line['SAMPLE_ID']
                            study = line['CANCER_TYPE']
                            sampleLabel[sampleID] = studyClass[study]
                    else:
                        if line['SAMPLE_ID'] in ASKeep:
                            sampleID = line['SAMPLE_ID']
                            study = line['CANCER_TYPE']
                            sampleLabel[sampleID] = studyClass[study]
                    
                for line in maf_in:
                    if line[:1] != '#':
                        names = line.split('\t')
                        break
                
                mafReader = csv.DictReader(maf_in, delimiter = '\t', fieldnames=names)

                
                
                for line in mafReader:
                    name = line['Tumor_Sample_Barcode']
                    Exon_Rank = line['Protein_position']
                    Gene_Name = line['Hugo_Symbol']
                    Amino_Acid_Change = line['HGVSp']
                    Transcript_ID = line['Transcript_ID']
                    Genotype = line['Tumor_Seq_Allele2']
                    Effect = line['Consequence']
                    chrom = 'chr' + str(line['Chromosome'])
                    try:
                        pos = int(line['Start_Position'])
                    except ValueError:
                        pos = 'NA'
                    ref = line['Reference_Allele']
                    
                    newVariant = variant(Amino_Acid_Change = Amino_Acid_Change, 
                                         Gene_Name = Gene_Name, Exon_Rank = Exon_Rank, Transcript_ID = Transcript_ID, 
                                         Genotype = Genotype, Effect = Effect, chrom = chrom, pos = pos, ref = ref, FASTA = hg19Fast)
                    try:
                        sampleVariants[name].append(newVariant)
                    except:
                        sampleVariants[name] = [newVariant]

                for name in sampleVariants:

                    try:
                        #print(name)
                        samples.append(sample(sampleVariants[name], name, 'Human', sampleLabel[name]))
                    except KeyError:
                        print('Sample name not found in phenotype file: ' + name)
                        
                        
        elif 'TCGA' in fileName:
            #assumes that TCGA filenames all follow the same naming convention
            #species = 'human'
            sampleVariants = {}
            
            label = studyClass[fileName.split('.')[1]]
            if label == 'Remove':
                continue
            
            print(label)
            
            with gzip.open(fileName, 'rt') as maf_in:
                for line in maf_in:
                    if line[:1] != '#':
                        names = line.split('\t')
                        break
                    
                mafReader = csv.DictReader(maf_in, delimiter = '\t', fieldnames=names)
                variantList = []
                oldName = ''
                
                
                for line in mafReader:
                    name = line['Tumor_Sample_UUID']
                    if oldName == '':
                        oldName = name
                                 
                    Amino_Acid_Change = line['Amino_acids']
                    Exon_Rank = line['Exon_Number']
                    Gene_Name = line['Gene']
                    Amino_Acid_Change = line['HGVSp']
                    Transcript_ID = line['Feature']
                    Genotype = line['Tumor_Seq_Allele2']
                    Effect = line['One_Consequence']
                    chrom = line['Chromosome']
                    pos = int(line['Start_Position'])
                    ref = line['Reference_Allele']
                    
                    newVariant = variant(Amino_Acid_Change = Amino_Acid_Change, 
                                         Gene_Name = Gene_Name, Exon_Rank = Exon_Rank, Transcript_ID = Transcript_ID, 
                                         Genotype = Genotype, Effect = Effect, chrom = chrom, pos = pos, ref = ref, FASTA = humFast)
                     
                    try:
                        sampleVariants[name].append(newVariant)
                    except KeyError:
                        sampleVariants[name] = [newVariant]

                for name in sampleVariants:
                    if len(sampleVariants[name]) < maxVar:
                        samples.append(sample(sampleVariants[name], name, 'Human', label))
        
        elif 'vcf' in fileName:
            if 'LSA' in fileName:
                label = 'dogLSA'
            elif 'MELANOMA' in fileName:
                label = 'dogMELANOMA'
            elif 'HSA' in fileName:
                label = 'dogHSA'
            elif 'OSA' in fileName:
                label = 'dogOSA'
                
            with VariantFile(fileName) as vcf_in:
                variants = []
            #split each record into its various transcripts (each record can list
            #effects for multiple transcripts
    
                for rec in vcf_in.fetch():
                
                    ref = rec.ref
                    chrom = rec.contig
                    pos = int(rec.pos)
                    
                    
                    #snpEFF = re.split('[,\'()]', str(rec.info['EFF']))
                    
                    
                    
                    transcripts = re.split('[,\'()]', str(rec.info['EFF']))
                    #print(str(transcripts) + '\n\n\n')
                    
                #split out fields from each transcript. Push to variants list
                    for transcript in transcripts:
                        fields = re.split('[|]', transcript)
                        if len(fields) == 1 and len(fields[0].strip()) > 0:
                            #print('effect: ' + transcript)
                            Effect = fields[0]
                        elif len(fields) > 1:
                            #print('multi-field: ' + transcript)
                            Effect_Impact = fields[0]
                            Functional_Class = fields[1]
                            Codon_change = fields[2]
                            Amino_Acid_Change = fields[3]
                            Amino_Acid_length = fields[4]
                            Gene_Name = fields[5]
                            Transcript_Biotype = fields[6]
                            Gene_Coding = fields[7]
                            Transcript_ID = fields[8]
                            Exon_Rank = fields[9]
                            Genotype = fields[10]
                            if getHumanGeneID(Gene_Name) != 'NA':
                                variants.append(variant(Effect_Impact, Functional_Class, Codon_change, Amino_Acid_Change, Amino_Acid_length, Gene_Name, Transcript_Biotype, Gene_Coding, Transcript_ID, Exon_Rank, Genotype, Effect, chrom, pos, ref, FASTA = dogFast))
                            #print(Genotype)
                            #if len(fields[5]) > 0:
                            #    genes.add(fields[5])
                        if (len(variants) > maxVar):
                            break
            if (len(variants) < maxVar):
                samples.append(sample(variants,fileName,args.species,label))
    
        elif 'tsv' in fileName:
            sampleVariants = {}
            #print(fileName)
            with gzip.open(fileName, 'rt') as ssm:
                label = studyClass[fileName.split('.')[2]]
                if label == 'Remove':
                    continue
                #label = label.split('-')[0]
                print(label)
                

                for line in ssm:
                    if line[:1] != '#':
                        names = line.split('\t')
                        break
                
                ssmreader = csv.DictReader(ssm, delimiter='\t', fieldnames=names)
                ssmSamples = []
                tsvList = []
                oldName = ''
                #mutationList = set()
                mutationID = ''
                
                for row in ssmreader:
                    if row['icgc_mutation_id'] == mutationID:
                        continue 
                    mutationID = row['icgc_mutation_id']
                    #mutationList.add(mutationID)
                    name = row['icgc_donor_id']
                    Effect = row['consequence_type']
                    Gene_Name = row['gene_affected']
                    Transcript_ID = row['transcript_affected']
                    Amino_Acid_Change = row['aa_mutation']
                    Codon_change = row['cds_mutation']
                    chrom = 'chr' + row['chromosome']
                    pos = int(row['chromosome_start'])
                    ref = row['reference_genome_allele']
                    Genotype = row['mutated_to_allele']
                    #mutation type is something I don't think we have for dog, but seems important (ex. 'single base substitution') should look into getting for dog
                    #also need "Gene_Coding, Exon_Rank, functional class
                    mutationType = row['mutation_type']
                    newVariant = variant(Codon_change = Codon_change, Amino_Acid_Change = Amino_Acid_Change, Gene_Name = Gene_Name, Transcript_ID = Transcript_ID, Genotype = Genotype, Effect = Effect, chrom = chrom, pos = pos, ref = ref, FASTA = hg19Fast)
                    
                    #newVariant = variant(Codon_change = Codon_change, Amino_Acid_Change = Amino_Acid_Change, Amino_Acid_length = Amino_Acid_length, Gene_Name = Gene_Name, Transcript_ID = Transcript_ID, Genotype = Genotype, Effect = Effect, chrom = chrom, pos = pos, ref = ref)
                    #Exon_Rank, Gene_Coding, Transcript_Biotype, Amino_Acid_length, functional_Class, Effect_Impact - all missing from ssm files. Need to make default entries for variant class
                    '''
                    if oldName == '':
                        oldName = name
                    if name == oldName:
                    '''
                    try:
                        sampleVariants[name].append(newVariant)
                    except KeyError:
                        sampleVariants[name] = [newVariant]
                        
                for name in sampleVariants:
                    if len(sampleVariants[name]) < maxVar:
                        samples.append(sample(sampleVariants[name], name, 'Human', label))
                        #if len(tsvList) > maxVar:
                        #    break
                        #add variant to running list for sample
                        #tsvList.append(newVariant)
                    #else:
                        #push completed sample with variant list to samples as a new sample
                        #clear variant list, and start new
                        #newSample = sample(tsvList, oldName, 'Human', label)
                        #print(oldName + '\t' + label)
                        #print(newSample)
                        #print(oldName)
                        #if (len(tsvList) < maxVar):
                            #samples.append(sample(tsvList, oldName, 'Human', label))
                        #tsvList = [newVariant]
                        #print('donor; ' + oldName + ' number of mutations: ' + str(mutationCount))
                        #oldName = name
                        #mutationCount = 0
        
    print('Number of broken genes in dataset: ' + str(len(brokenGenes.keys())))

#####add a keys list with the following entries. Will give an order for looping through and outputting later. Not adding now, because I'm offline and want to just copy from the existing one.
contextKeys = [
            'A[C>A]A','A[C>A]C','A[C>A]G','A[C>A]T','C[C>A]A','C[C>A]C','C[C>A]G','C[C>A]T',
            'G[C>A]A','G[C>A]C','G[C>A]G','G[C>A]T','T[C>A]A','T[C>A]C','T[C>A]G','T[C>A]T',
            'A[C>G]A','A[C>G]C','A[C>G]G','A[C>G]T','C[C>G]A','C[C>G]C','C[C>G]G','C[C>G]T',
            'G[C>G]A','G[C>G]C','G[C>G]G','G[C>G]T','T[C>G]A','T[C>G]C','T[C>G]G','T[C>G]T',
            'A[C>T]A','A[C>T]C','A[C>T]G','A[C>T]T','C[C>T]A','C[C>T]C','C[C>T]G','C[C>T]T',
            'G[C>T]A','G[C>T]C','G[C>T]G','G[C>T]T','T[C>T]A','T[C>T]C','T[C>T]G','T[C>T]T',
            'A[T>A]A','A[T>A]C','A[T>A]G','A[T>A]T','C[T>A]A','C[T>A]C','C[T>A]G','C[T>A]T',
            'G[T>A]A','G[T>A]C','G[T>A]G','G[T>A]T','T[T>A]A','T[T>A]C','T[T>A]G','T[T>A]T',
            'A[T>G]A','A[T>G]C','A[T>G]G','A[T>G]T','C[T>G]A','C[T>G]C','C[T>G]G','C[T>G]T',
            'G[T>G]A','G[T>G]C','G[T>G]G','G[T>G]T','T[T>G]A','T[T>G]C','T[T>G]G','T[T>G]T',
            'A[T>C]A','A[T>C]C','A[T>C]G','A[T>C]T','C[T>C]A','C[T>C]C','C[T>C]G','C[T>C]T',
            'G[T>C]A','G[T>C]C','G[T>C]G','G[T>C]T','T[T>C]A','T[T>C]C','T[T>C]G','T[T>C]T',
            ]
for context in contextKeys:
    featureList.append(context)
    
if args.triNucCon:
    with open(args.outFile + '_triCon.txt', 'w') as contextOutFile:
        
        contextOutFile.write('disease')
        
        for context in contextKeys:
            contextOutFile.write('\t' + context)
        contextOutFile.write('\n')
        
        for sample in samples:
            contextOutFile.write(sample.label)
            triContextDict = {
            'A[C>A]A':0,'A[C>A]C':0,'A[C>A]G':0,'A[C>A]T':0,'C[C>A]A':0,'C[C>A]C':0,'C[C>A]G':0,'C[C>A]T':0,
            'G[C>A]A':0,'G[C>A]C':0,'G[C>A]G':0,'G[C>A]T':0,'T[C>A]A':0,'T[C>A]C':0,'T[C>A]G':0,'T[C>A]T':0,
            'A[C>G]A':0,'A[C>G]C':0,'A[C>G]G':0,'A[C>G]T':0,'C[C>G]A':0,'C[C>G]C':0,'C[C>G]G':0,'C[C>G]T':0,
            'G[C>G]A':0,'G[C>G]C':0,'G[C>G]G':0,'G[C>G]T':0,'T[C>G]A':0,'T[C>G]C':0,'T[C>G]G':0,'T[C>G]T':0,
            'A[C>T]A':0,'A[C>T]C':0,'A[C>T]G':0,'A[C>T]T':0,'C[C>T]A':0,'C[C>T]C':0,'C[C>T]G':0,'C[C>T]T':0,
            'G[C>T]A':0,'G[C>T]C':0,'G[C>T]G':0,'G[C>T]T':0,'T[C>T]A':0,'T[C>T]C':0,'T[C>T]G':0,'T[C>T]T':0,
            'A[T>A]A':0,'A[T>A]C':0,'A[T>A]G':0,'A[T>A]T':0,'C[T>A]A':0,'C[T>A]C':0,'C[T>A]G':0,'C[T>A]T':0,
            'G[T>A]A':0,'G[T>A]C':0,'G[T>A]G':0,'G[T>A]T':0,'T[T>A]A':0,'T[T>A]C':0,'T[T>A]G':0,'T[T>A]T':0,
            'A[T>G]A':0,'A[T>G]C':0,'A[T>G]G':0,'A[T>G]T':0,'C[T>G]A':0,'C[T>G]C':0,'C[T>G]G':0,'C[T>G]T':0,
            'G[T>G]A':0,'G[T>G]C':0,'G[T>G]G':0,'G[T>G]T':0,'T[T>G]A':0,'T[T>G]C':0,'T[T>G]G':0,'T[T>G]T':0,
            'A[T>C]A':0,'A[T>C]C':0,'A[T>C]G':0,'A[T>C]T':0,'C[T>C]A':0,'C[T>C]C':0,'C[T>C]G':0,'C[T>C]T':0,
            'G[T>C]A':0,'G[T>C]C':0,'G[T>C]G':0,'G[T>C]T':0,'T[T>C]A':0,'T[T>C]C':0,'T[T>C]G':0,'T[T>C]T':0,
            }
            
            snpCount = 0
            
            
            for variant in  sample.variants:
                #if variant.variantID not in variantSet:
                try:
                    triContextDict[variant.trinucleotide] += 1
                    #variantSet.add(variant.variantID)
                        #snpCount += 1
                except KeyError:
                    if variant.trinucleotide != 'NA':
                        print('This ain\'t got no fitting in gud!')
                        print(variant.trinucleotide)
            #sample.VariantCount = len(variantSet)
            #print(sample.variants)
            mutCount = 0
            for context in contextKeys:
                mutCount += triContextDict[context]
                #print(mutCount)
            for context in contextKeys:
                try:
                    #print(str(triContextDict[context]/mutCount))               
                    contextOutFile.write('\t' + str(triContextDict[context]/mutCount))
                    sample.toWrite.append(triContextDict[context]/mutCount)
                except ZeroDivisionError:
                    contextOutFile.write('\t' + '0')
                    sample.toWrite.append(0)
            contextOutFile.write('\n')
                
        
        
            
        

###Report binary broken/not broken for genes    
if args.brokenGene:
    #broked = list(getBrokenGenes(samples))
    #p = multiprocessing.Pool(4)
    #broked = brokenGenes
    print('broked: '+ str(len(brokenGenes)))
    if args.sparsityMin > 0:
        
        brokenAtLeastMin = []
        for gene in brokenGenes:
            mutCount = 0
            for sample in samples:
                if isGeneBroken(sample, gene):
                    mutCount +=1
                if mutCount >= args.sparsityMin:
                    brokenAtLeastMin.append(gene)
                    break
    else: 
        brokenAtLeastMin = brokenGenes
        #if mutCount > 5:
        #    brokenAtLeastMin.append(gene)
        #    continue
    for geneID in brokenAtLeastMin:
        featureList.append(geneID)
    
    
                
    #print(len(broked))
    with open(args.outFile + 'testplus.out', 'w') as testOut:
        print('genes broken in at least ' + str(args.sparsityMin) + ' samples: ' + str(len(brokenAtLeastMin)))
        testOut.write('disease')
        for geneID in brokenAtLeastMin:
            try:
                testOut.write('\t' + getOrtholog(geneID).dogGeneName)
            except: 
                testOut.write('\t' + getOrtholog(geneID))
        testOut.write('\n')
        for sample in samples:
            brokenInSample = []
            for geneID in brokenAtLeastMin:
                
                if isGeneBroken(sample, geneID):
                    brokenInSample.append(1)
                    sample.toWrite.append(1)
                else:
                    brokenInSample.append(0)
                    sample.toWrite.append(0)
                #geneCount = 0
                #for gene in brokenInSample:
                #    geneCount += gene
            #print('sample Name: ' + str(sample.name))
            testOut.write(sample.label)
            for gene in brokenInSample:
                testOut.write('\t' + str(gene))
            testOut.write('\n')

if args.brokenPathways:
    broked = brokenGenes
    
    brokenPathways = set()
    for gene in broked:
        try:
            for pathway in reference[gene].pathways:
                brokenPathways.add(pathway[0])
        except:
            print(gene)
            
    print('Number of broken pathways: ' + str(len(brokenPathways)))
    
    brokenPathways = list(brokenPathways)
    for pathway in brokenPathways:
        featureList.append(pathway)
    
    with open(args.outFile + 'testPathways.out', 'w') as testPathOut:
        testPathOut.write('disease')
        for pathway in brokenPathways:
            testPathOut.write('\t' + pathway)
        testPathOut.write('\n')
        for sample in samples:
            brokenPathInSample = []
            for pathway in brokenPathways:
                
                if pathway in sample.brokenPathways:
                #if isPathwayBroken(sample, pathway):
                    brokenPathInSample.append(1)
                    sample.toWrite.append(1)
                else:
                    brokenPathInSample.append(0)
                    sample.toWrite.append(0)
            testPathOut.write(sample.label)
            for pathway in brokenPathInSample:
                testPathOut.write('\t' + str(pathway))
            testPathOut.write('\n')
            
if args.mutationCount:
    featureList.append('mutCount')
    with open(args.outFile + 'VariantCounts.out', 'w') as varCountOut:
        for sample in samples:           
            variantSet = set()
            for variant in sample.variants:
                variantSet.add(variant.variantID)
            sample.variantCount = len(variantSet)
            varCountOut.write(sample.label + '\t' + str(sample.variantCount) + '\n')
            sample.toWrite.append(sample.variantCount)
                
if args.geneCount:
    print('broked genes to count: '+ str(len(brokenGenes)))
    if args.sparsityMin > 0:
        
        brokenAtLeastMin = []
        for gene in brokenGenes:
            mutCount = 0
            for sample in samples:
                if isGeneBroken(sample, gene):
                    mutCount +=1
                if mutCount >= args.sparsityMin:
                    brokenAtLeastMin.append(gene)
                    break
    else: 
        brokenAtLeastMin = brokenGenes
        #if mutCount > 5:
        #    brokenAtLeastMin.append(gene)
        #    continue
    for geneID in brokenAtLeastMin:
        featureList.append(geneID + '_count')
    
    
                
    #print(len(broked))
    with open(args.outFile + 'brokenGenesCount.out', 'w') as brokenCountOut:
        print("Number of genes testing: " + str(len(brokenAtLeastMin)))
        brokenCountOut.write('disease')
        for geneID in brokenAtLeastMin:
            try:
                brokenCountOut.write('\t' + getOrtholog(geneID).dogGeneName)
            except: 
                brokenCountOut.write('\t' + getOrtholog(geneID))
        brokenCountOut.write('\n')
        for sample in samples:
            brokenInSample = []
            for geneID in brokenAtLeastMin:
                
                if isGeneBroken(sample, geneID):
                    
                    brokenInSample.append(sample.brokenGeneCount[geneID])
                    sample.toWrite.append(sample.brokenGeneCount[geneID])
                else:
                    brokenInSample.append(0)
                    sample.toWrite.append(0)
                #geneCount = 0
                #for gene in brokenInSample:
                #    geneCount += gene
            #print('sample Name: ' + str(sample.name))
            brokenCountOut.write(sample.label)
            for gene in brokenInSample:
                brokenCountOut.write('\t' + str(gene))
            brokenCountOut.write('\n')
            
if args.pathwayCount:
    brokenPathways = set()
    for gene in brokenGenes:
        try:
            for pathway in reference[gene].pathways:
                brokenPathways.add(pathway[0])
        except:
            print(gene)
            
    print('Number of broken pathways: ' + str(len(brokenPathways)))
    
    brokenPathways = list(brokenPathways)
    for pathway in brokenPathways:
        featureList.append(pathway + '_count')
    
    with open(args.outFile + 'PathwaysCount.out', 'w') as pathCountOut:
        pathCountOut.write('disease')
        for pathway in brokenPathways:
            pathCountOut.write('\t' + pathway)
        pathCountOut.write('\n')
        for sample in samples:
            brokenPathInSample = []
            for pathway in brokenPathways:
                
                if pathway in sample.brokenPathways:
                #if isPathwayBroken(sample, pathway):
                    brokenPathInSample.append(sample.brokenPathwayCount[pathway])
                    sample.toWrite.append(sample.brokenPathwayCount[pathway])
                else:
                    brokenPathInSample.append(0)
                    sample.toWrite.append(0)
            pathCountOut.write(sample.label)
            for pathway in brokenPathInSample:
                pathCountOut.write('\t' + str(pathway))
            pathCountOut.write('\n')
    
                
            
            
with open(args.outFile + 'CumulativeForLearning.out', 'w') as learnOut:
    #labels = list(labels)
    learnOut.write('label')
    for feature in featureList:
        learnOut.write('\t' + feature)
    learnOut.write('\n')
    for sample in samples:
        learnOut.write(sample.label)
        for feature in sample.toWrite:
            learnOut.write('\t' + str(feature))
        learnOut.write('\n')
    
    

