import pysam
from pysam import VariantFile, FastaFile
import re
import csv
import argparse
import gzip
import os
import sys
import gffutils
import subprocess
import pandas as pd




##Example usage:
#python parseVCFminimizingUpdate061621_updatingForMultipleOrthologs.py -r /seq/vgb/swofford/dogHumanOrthologs/old/feb2/ncbiGeneNameAdded.bed -i /seq/vgb/cancer_r01/cbioportal/metastatic_prostate_cancer/MPCP_provisional_nov_2019/prad_mpcproject_2018/data_mutations_mskcc.txt /seq/vgb/cancer_r01/cbioportal/metastatic_breast_cancer/MBCP_provisional_feb_2020/brca_mbcproject_wagle_2017/data_mutations_mskcc.txt /seq/vgb/cancer_r01/cbioportal/metastatic_solid_umich_2017/metastatic_solid_tumors_mich_2017/data_mutations_mskcc.txt /seq/vgb/cancer_r01/cbioportal/pediatric_columbia_2019/mixed_pipseq_2017/data_mutations_mskcc.txt /seq/vgb/cancer_r01/cbioportal/angiosarcoma_project_july_2020/angs_painter_2020/data_mutations_extended.txt /seq/vgb/cancer_r01/cbioportal/pediatric_dkfz/data_mutations_mskcc.txt /seq/vgb/cancer_r01/icgc/release_28/simple_somatic_mutation_files/*.tsv.gz /seq/vgb/cancer_r01/tcga_data/simple_somatic_mutations/*.somatic.maf.gz /seq/vgb/cancer_r01/cbioportal/dfci_dlbcl_2018/dlbcl_dfci_2018/data_mutations_mskcc.txt /seq/vgb/cancer_r01/cbioportal/broad_dlbcl_2012/dlbc_broad_2012/data_mutations.txt /seq/vgb/cancer_r01/icgc/release_28/simple_somatic_mutation_files/simple_somatic_mutation.open.BOCA-UK.tsv.gz /seq/vgb/swofford/R01/convertingToSKLInput/testVCFs/BCL/snpEFFrerun/*.vcf /seq/vgb/swofford/R01/convertingToSKLInput/testVCFs/CMT/*.vcf /seq/vgb/cancer_r01/canine_glioma_verhaak/maf_nogistic_maftools_weirdPosRemoved.maf /seq/vgb/swofford/R01/convertingToSKLInput/testVCFs/HSA/redo/*.vcf /seq/vgb/swofford/R01/convertingToSKLInput/testVCFs/LSA/redo/comboTest/*.vcf /seq/vgb/swofford/R01/convertingToSKLInput/testVCFs/OSA/rename/*.vcf /seq/vgb/swofford/R01/convertingToSKLInput/testVCFs/OSA/HeatherOSA.txt /seq/vgb/swofford/R01/convertingToSKLInput/testVCFs/MELANOMA/*.vcf -o test080322 -bp -tri -g -m -gc -pc
#python parseVCF.py -i /seq/vgb/cancer_r01/cbioportal/metastatic_prostate_cancer/MPCP_provisional_nov_2019/prad_mpcproject_2018/data_mutations_mskcc.txt /seq/vgb/cancer_r01/cbioportal/metastatic_breast_cancer/MBCP_provisional_feb_2020/brca_mbcproject_wagle_2017/data_mutations_mskcc.txt /seq/vgb/cancer_r01/cbioportal/metastatic_solid_umich_2017/metastatic_solid_tumors_mich_2017/data_mutations_mskcc.txt /seq/vgb/cancer_r01/cbioportal/pediatric_columbia_2019/mixed_pipseq_2017/data_mutations_mskcc.txt /seq/vgb/cancer_r01/cbioportal/angiosarcoma_project_july_2020/angs_painter_2020/data_mutations_extended.txt /seq/vgb/swofford/R01/convertingToSKLInput/testVCFs/BCL/snpEFFrerun/*.vcf /seq/vgb/swofford/R01/convertingToSKLInput/testVCFs/CMT/*.vcf /seq/vgb/cancer_r01/canine_glioma_verhaak/maf_nogistic_maftools_weirdPosRemoved.maf ./HSA/*.vcf ./LSA/*.vcf ./OSA/*fixed.vcf ./OSA/HeatherOSA.txt ./MELANOMA/*.vcf /seq/vgb/cancer_r01/cbioportal/angiosarcoma/data_mutations_mskcc.txt /seq/vgb/cancer_r01/cbioportal/pediatric_dkfz/data_mutations_mskcc.txt /seq/vgb/cancer_r01/icgc/release_28/simple_somatic_mutation_files/*.tsv.gz /seq/vgb/cancer_r01/tcga_data/simple_somatic_mutations/*.somatic.maf.gz -r /seq/vgb/swofford/ref/EnsemblHumanDogOrthologs.txt -p /seq/vgb/swofford/ref/pathways/C2CuratedGeneSetsGeneNames.gmt.txt -o allTestNoCounts -hf /seq/vgb/swofford/ref/hg38.fa -df /seq/vgb/swofford/ref/Canis_lupus_familiaris_assembly3.fasta -hg19 /seq/vgb/swofford/ref/hg19plusChrM.fa -bp -tri -g
# python parseVCF.py -i ../SNPeff_output/LymphomaForTesting/*.vcf /seq/vgb/cancer_r01/cbioportal/angiosarcoma/data_mutations_mskcc.txt /seq/vgb/cancer_r01/icgc/*.tsv.gz /seq/vgb/cancer_r01/tcga_data/maf_files/*/*.somatic.maf.gz -r /seq/vgb/swofford/ref/EnsemblHumanDogOrthologs.txt -p /seq/vgb/swofford/ref/pathways/C2CuratedGeneSetsGeneNames.gmt.txt -o allTest -hf /seq/vgb/swofford/ref/hg38.fa -df /seq/vgb/swofford/ref/Canis_lupus_familiaris_assembly3.fasta -hg19 /seq/vgb/swofford/ref/hg19plusChrM.fa -bp -tri -m -gc

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--inFile', nargs='+')
parser.add_argument('-sp', '--sparsityMin', default = 1)
parser.add_argument('-o', '--outFile')
#parser.add_argument('-s', '--species')
parser.add_argument('-r', '--refFile', default = '/seq/vgb/swofford/dogHumanOrthologs/old/feb2/ncbiGeneNameAdded_zeroRemoved_030623.bed')
#parser.add_argument('-r', '--refFile', default = '/seq/vgb/swofford/dogHumanOrthologs/old/feb2/ncbiGeneNameAdded_chrMAndUnanchoredRemoved.bed')
parser.add_argument('-p', '--pathways', default = '/seq/vgb/swofford/ref/pathways/C2CuratedGeneSetsGeneNames.gmt.txt')
parser.add_argument('-bp', '--brokenPathways', action='store_true', default=False)
parser.add_argument('-g', '--brokenGene', action='store_true', default=False)
parser.add_argument('-nd', '--nonDamage', action='store_true', default=False)
parser.add_argument('-ndc', '--nonDamageGeneCount', action='store_true', default=False)
parser.add_argument('-maxVar', '--maximumVariantsSample', default = 10000000)
parser.add_argument('-minVar', '--minimumVariantsSample', default = 5)
parser.add_argument('-m', '--mutationCount', action='store_true', default=False)
parser.add_argument('-hf', '--humanFasta', default = '/seq/vgb/swofford/ref/hg38.fa')
parser.add_argument('-hg19', '--hg19Fasta', default = '/seq/vgb/swofford/ref/hg19plusChrM.fa')
parser.add_argument('-df', '--dogFasta', default = '/seq/vgb/swofford/ref/Canis_lupus_familiaris_assembly3.fasta')
parser.add_argument('-tri', '--triNucCon', action = 'store_true', default=False)
parser.add_argument('-gc', '--geneCount', action = 'store_true' , default=False)
parser.add_argument('-pc', '--pathwayCount', action = 'store_true', default=False)
parser.add_argument('-gr', '--granularity', type = int, default = 4)
#parser.add_argument('-hg38bed', '--hg38bed', default = '/seq/vgb/swofford/R01/convertingToSKLInput/testVCFs/hg38CDS_chrM_unanch_removed_merged.bed')
#parser.add_argument('-hg19bed', '--hg19bed', default = '/seq/vgb/swofford/R01/convertingToSKLInput/testVCFs/hg19CDS_chrM_unanch_removed_merged.bed')
#parser.add_argument('-canFam3bed', '--canfam3bed', default = '/seq/vgb/swofford/R01/convertingToSKLInput/testVCFs/canFam3CDS_chrM_unanch_removed_merged.bed')
parser.add_argument('-hg38bed', '--hg38bed', default = '/seq/vgb/swofford/R01/convertingToSKLInput/testVCFs/hg38CDSMerged.bed')
parser.add_argument('-hg19bed', '--hg19bed', default = '/seq/vgb/swofford/R01/convertingToSKLInput/testVCFs/hg19CDSMerged.bed')
parser.add_argument('-canFam3bed', '--canfam3bed', default = '/seq/vgb/swofford/R01/convertingToSKLInput/testVCFs/canFam3CDSMerged.bed')
parser.add_argument('-noFly', '--noFlyList', default = '/seq/vgb/jason/cancer_mutations/cancer_machine_learning/duplicate_issue/samples.remove')
parser.add_argument('-db', '--dogDb', default = '/seq/vgb/swofford/ref/canFam3.104.db')
parser.add_argument('-hg38db', '--hg38Db', default = '/seq/vgb/swofford/ref/hg38.104.db')
parser.add_argument('-hg19db', '--hg19Db', default = '/seq/vgb/swofford/ref/hg19.87.db')
parser.add_argument('-dm', '--driverMutCount', action = 'store_true', default=False)
parser.add_argument('-s', '--source', action = 'store_true', default=False)


args = parser.parse_args()



'''Samples are stored as an instance of the sample class. Each sample has a label which
will be used for learning at a later stage, a dictionary of variants, a species (currently poorly implemented)
a list of broken genes to speed lookup, a count of broken genes, a list of broken pathways, a count of variants,
and a list of features to be written to an output file. feature list is added after parsing based on options.
'''
 
class sample:
    def __init__(self, variants, name, species, label = 'NA', assembly = 'NA', source = 'NA'):
        pathBrokenInSample = {}
        genesBrokenInSample = {}
        genesWithNonDamagingMuts = {}
        ####
        variantSet = set()
        codingVariants = []
        outVariants = {}
        outCodons = {}
        driverMutCount = 0
        
        
        for variant in variants:
            ###
            if variant.variantID not in variantSet:# and variant.geneID != 'NA':
            ###
                if any(id in variant.geneID for id in driverIDs):
                    driverMutCount += 1
                '''    
                for id in list(variant.geneID):
                    if id in driverIDs:
                        driverMutCount += 1
                        break
                '''
                #if variant.geneID in driverIDs:
                #    driverMutCount += 1
                #else:
                #    print(variant.geneID)
                    
                #print('variant ID = ' + variant.variantID)
                if '-1' not in variant.outVariantID:
                    outVariants[variant.outVariantID] = [variant.ref, variant.Genotype]
                    outVariantSet.add(variant.outVariantID)
                    '''
                    try:
                        for codon in variantToCodon[variant.outVariantID]:
                            outCodonSet.add(codon)
                            outCodons[codon] = 1
                    except KeyError:
                        pass
                    '''
                if assembly == 'NA':
                    assembly = variant.assembly
                variantSet.add(variant.variantID)
                codingVariants.append(variant)
                
                if variant.breaksGene:
                    if variant.geneID != 'NA':
                        for id in variant.geneID:
                            try:
                                genesBrokenInSample[id] += 1                    
                            except KeyError:
                                genesBrokenInSample[id] = 1
                                #print('do I ever get here?')
                            for pathway in reference[id].pathways:
                                try:
                                    pathBrokenInSample[pathway] += 1
                                except KeyError:
                                    pathBrokenInSample[pathway] = 1
                else:
                    if variant.geneID != 'NA':
                        for id in variant.geneID:
                            try:
                                genesWithNonDamagingMuts[id] += 1                    
                            except KeyError:
                                genesWithNonDamagingMuts[id] = 1
                                                           
        
        self.species = species
        self.name = name
        self.variants = codingVariants                        
        self.brokenGenes = genesBrokenInSample.keys()
        self.nonDamageGenes = genesWithNonDamagingMuts.keys()
        self.nonDamageGeneCount = genesWithNonDamagingMuts
        self.brokenGeneCount = genesBrokenInSample
        self.brokenPathways = pathBrokenInSample.keys()
        self.brokenPathwayCount = pathBrokenInSample
        self.variantCount = len(variantSet)
        self.toWrite = []
        self.label = label
        self.assembly = assembly
        self.outVariants = outVariants
        self.outCodons = outCodons
        self.driverMutCount = driverMutCount
        self.source = source
        labels.add(label)



class variant:
    
    def __init__(self, Effect_Impact = 'NA', Functional_Class = 'NA', Codon_change = 'NA', Amino_Acid_Change = 'NA', 
                 Amino_Acid_length = 'NA', Gene_Name = 'NA', Transcript_Biotype = 'NA', Gene_Coding = 'NA', 
                 Transcript_ID = 'NA', Exon_Rank = 'NA', Genotype = 'NA', Effect = 'NA', chrom = 'NA', 
                 pos = 'NA', ref = 'NA', assembly='NA', changeTri='NA', hg38Chr = 'NA', hg38Pos = 'Na'):
        
        if 'ENS' in Transcript_ID:
            #print(Transcript_ID)
            if assembly == 'hg19':
                db = hg19db
                '''
                try:
                    chr, start, end = getLiftover(chrom, pos - 1, pos, hg19ToHg38)
                    self.hg38Chr = chr
                    self.hg38Pos = end
                except:
                    self.hg38Chr = 'NA'
                    self.hg38Pos = -1
                '''
            elif assembly == 'hg38':
                db = hg38db
                '''
                self.hg38Chr = chrom
                self.hg38Pos = pos
                '''
            else:
                db = dogdb
                '''
                try:
                    chr, start, end = getLiftover(chrom, pos - 1, pos, canFam3ToHg38)
                    self.hg38Chr = chr
                    self.hg38Pos = end
                except:
                    self.hg38Chr = 'NA'
                    self.hg38Pos = -1
                '''
            try:
                Gene_Name = db[Transcript_ID.split('.')[0]]['gene_id'][0]
            except:
                Gene_Name = Gene_Name
        geneID = getHumanGeneID(Gene_Name)
        
        self.chrom = chrom
        self.pos = pos
        self.Effect = Effect
        self.Effect_Impact = Effect_Impact #don'tNeed
        self.Functional_Class = Functional_Class #don'tNeed
        self.Codon_change = Codon_change #don'tNeed
        self.Amino_Acid_Change = Amino_Acid_Change #don'tNeed
        self.Amino_Acid_length = Amino_Acid_length #don'tNeed
        self.Gene_Name = Gene_Name
        self.Transcript_BioType = Transcript_Biotype #don'tNeed 
        self.Gene_Coding = Gene_Coding #don'tNeed
        self.Transcript_ID = Transcript_ID
        self.Exon_Rank = Exon_Rank #don'tNeed
        self.Genotype = Genotype.upper()
        self.ref = ref.upper()
        self.geneID = geneID
        self.dogGeneID = getDogGeneID(Gene_Name)
        self.variantID = chrom+str(pos)+Genotype.upper()
        try:
            self.outVariantID = variantIDLookup['@'.join([chrom,str(pos),assembly])]
        except KeyError:
            self.outVariantID = 'NA'
        #self.outVariantID = '@'.join([chrom,str(pos),assembly])
        #self.outVariantID = hg38Chr+str(hg38Pos)+Genotype.upper()
        self.assembly = assembly
        self.refTrinucleotide = ''
        if assembly == 'hg38':
            FASTA = humFast
        elif assembly == 'hg19':
            FASTA = hg19Fast
        elif assembly == 'canFam3':
            FASTA = dogFast
        
        if ('missense' in Effect 
            or 'frameshift_variant' in Effect 
            or 'disruptive_inframe_deletion' in Effect 
            or 'stop_gained' in Effect 
            or 'stop_lost' in Effect
            or 'MISSENSE' in Effect
            or 'NONSENSE' in Effect
            or 'nonsense' in Effect
            or Effect in damagingEffects
            and geneID != 'NA'):
                #print(Effect)
                self.breaksGene = True            
                for id in geneID:
                    try:
                        brokenGenes[id] += 1
                    except KeyError:
                        brokenGenes[id] = 1
        else:
            #print('effect = ' + Effect)
            self.breaksGene = False
            if geneID != 'NA':
                for id in geneID:
                    try:
                        nonDamageGenes[id] += 1
                    except:
                        nonDamageGenes[id] = 1
                    
        if len(ref) == 1 and len(Genotype) == 1 and ref != '-' and Genotype != '-':
            
            refTrinucleotide = FASTA.fetch(chrom,pos-2,pos+1).upper() 
            
            ref = ref.upper()
            
            if ref == 'A' or ref == 'G':
                self.refTrinucleotide = reverse(complement(refTrinucleotide))
                refTrinucleotid = reverse(complement(refTrinucleotide))
                self.ref = complement(ref)
                ref = complement(ref)
                self.Genotype = complement(Genotype)
                Genotype = complement(Genotype)
                
            else:
                self.refTrinucleotide = refTrinucleotide
    
                
            changeTri = refTrinucleotide[:1] + '[' + ref + '>' + Genotype + ']' + refTrinucleotide[2:]
            
            if refTrinucleotide[1:2] != ref:
                if reverse(complement(refTrinucleotide))[1:2] == ref:
                    self.refTrinucleotide = reverse(complement(refTrinucleotide))
                    refTrinucleotid = reverse(complement(refTrinucleotide))
                    changeTri = refTrinucleotide[:1] + '[' + ref + '>' + Genotype + ']' + refTrinucleotide[2:]
                else:
                    print('Ok, weirdness is still happening')
                    print('derived ref : ' + refTrinucleotide + ' ref according to input file: ' + ref + ' genotype from file: ' + Genotype)
                    print(self.chrom + ': ' + str(self.pos))
                    #print('huh, reverse compliment worked')
                #print('AAAAAHHHHH WHAT DO YOU DO TO THE BRAIN')
                #print('derived ref : ' + refTrinucleotide + ' ref according to input file: ' + ref + ' genotype from file: ' + Genotype)
                #print(ref)
                #print(Genotype)
            #else: print('good')           
        #else:
            #print('ref: ' + ref + ' genotype: ' + Genotype)        
        #print(changeTri)             
        self.trinucleotide = changeTri
       
        
        


#SNPeff output has the following in the EFF info section:
#Type=String,Description="Predicted effects for this variant.Format: 'Effect ( Effect_Impact | Functional_Class | Codon_Change | Amino_Acid_Change| Amino_Acid_length | Gene_Name | Transcript_BioType | Gene_Coding | Transcript_ID | Exon_Rank  | Genotype [ | ERRORS | WARNINGS ] )'">

#icgc output guide chrom = chromosome, pos = chromosome_start, Effect = consequence_type, ref = reference_genome_allele, Genotype = mutated_to_allele

''' Class defining a single gene ortholog relationship between dog and human, to facilitate bi-directional lookups later'''
class ortholog:
    def __init__(self, dogGeneName = '', humanGeneName = set(), dogGeneID = set(), geneID = set(), pathways = set()):
        self.dogGeneName = dogGeneName
        self.humanGeneName = humanGeneName
        self.dogGeneID = dogGeneID
        self.geneID = geneID
        self.pathways = pathways

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

def hammingDist(seq1, seq2):
    '''Determines the Hamming distance between two strings. Equivalent to the
    minimum number of character substitutions required to change from one
    string to the other. Takes two strings as input, and returns an integer
    value for the Hamming distance'''
    return sum(ch1 != ch2 for ch1, ch2 in zip(seq1, seq2))

def getLiftover(chr, start, end, chain):
    outChr = ''
    outStart = ''
    outEnd = ''
    with open('tempToLift.bed', 'w') as outFile:
            outFile.write('\t'.join([chr, str(start), str(end)]) + '\n') #str(pos - 1), str(pos)]) + '\n')
    subprocess.run(['liftOver', 'tempToLift.bed', chain, 'tempLifted.bed', 'tempUnmap.bed'])
    with open('tempLifted.bed', 'r') as inFile:
        try:
            line = next(inFile)
            bed = line.strip().split('\t')
            outChr = bed[0]
            outStart = bed[1]
            outEnd = bed[2]
        except StopIteration: #error thrown when file is empty/liftover failed
            outChr = 'NA'
            outStart = 0
            outEnd = 0
    subprocess.run(['rm', 'tempToLift.bed', 'tempLifted.bed', 'tempUnmap.bed'])
    return(outChr, int(outStart), int(outEnd))

def getFasta(chr, start, end, fasta):
    try:
        return fasta.fetch(chr, start, end).upper()
    except:
        return 'NA'



'''builds a reference table lookup between human/dog orthologous genes, parsed from an Ensembl table. File must contain the below columns, extra columns are fine and will be ignored.'''
def buildRef(refFile, pathways):
    with open(refFile, 'r') as refFile, open(pathways, 'r') as pathRef:
        ref = {}
        pathwayList = {}
        genePath = {}
        for line in pathRef:
            genes = line.strip().split('\t')
            pathName = genes[0]
            for gene in genes[2:]:
                try:
                    if 'http' in gene:
                        print(gene)
                    genePath[gene].append(pathName)
                except KeyError:
                    genePath[gene] = [pathName]
        print(str(len(genePath)) + ' genes with pathways found.')                 
        for line in refFile:
            humEns = set()
            dogEns = set()
            humanGeneName = set()
            pathwaysEffected = set()
            if 'Ensembl' in line or 'NCBI' in line:
                fields = line.strip().split()
                meta = fields[3].split(':')
                if meta[1] == '' or meta[3] == '':
                    continue
                dogEnsembl = meta[1].strip()
                if dogEnsembl != '':
                    dogEns.add(dogEnsembl)
                geneID = meta[3].strip()
                if geneID != '':
                    humEns.add(geneID)
                geneName = meta[2].strip()
                if geneName != '':
                    humanGeneName.add(geneName)
                dogGeneName = meta[0].strip()
                ncbiDogID = meta[6].strip()
                if ncbiDogID != '':
                    dogEns.add(ncbiDogID)
                ncbiHumID = meta[5].strip()
                if ncbiHumID != '':
                    humEns.add(ncbiHumID)
                ncbiHumName = meta[9].strip()
                if ncbiHumName != '':
                    humanGeneName.add(ncbiHumName)
                pathwaysEffected = set()
                for i in humanGeneName:
                    try:
                        pathwaysEffected.update(set(genePath[i]))
                    except KeyError:
                        continue                
                orth = ortholog(dogGeneName, humanGeneName, dogEns, humEns, pathways = pathwaysEffected)
                for id in humEns:
                    try:
                        ref[id].pathways.update(pathwaysEffected)
                        ref[id].geneID.update(humEns)
                        ref[id].humanGeneName.update(humanGeneName)
                    except KeyError:
                        ref[id] = orth
                for name in humanGeneName:
                    if name != '' and name != 'NA':
                        try:
                            ref[name].pathways.update(pathwaysEffected)
                            ref[name].geneID.update(humEns)
                            ref[name].humanGeneName.update(humanGeneName)
                        except KeyError:
                            ref[name] = orth 
                for id in dogEns:
                    try:
                        ref[id].pathways.update(pathwaysEffected)
                        ref[id].geneID.update(humEns)
                        ref[id].humanGeneName.update(humanGeneName)
                    except KeyError:
                        ref[id] = orth    
        print('Finished building Reference. ' + str(len(ref)/3) + ' orthologs recorded.')
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
    return(genes)


'''queries whether or not a given gene is broken in an individual sample, only considers genes with human/dog ortholog pairs'''
def isGeneBroken(sample, geneID):
    if geneID in sample.brokenGenes:
        return True
    return False

def isGeneNonDamage(sample, geneID):
    if geneID in sample.nonDamageGenes:
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

def inCDS(assembly, chrom, pos):
    if isNumber(pos):
        pos = int(pos)
    else: 
        return False
    if chrom == 'chrNA':
        return False
    #print(assembly + '\t' + chrom + '\t' + str(pos))
    if assembly == 'dogFast':
        ref = canFam3CDS
    elif assembly == 'hg19Fast':
        ref = hg19CDS
    elif assembly == 'humFast':
        ref = hg38CDS        
    if chrom == 'chrM' or chrom == 'chrMT':
        return False
    try:    
        for interval in ref[chrom]:
            if pos > interval[1]:
                #print(str(pos) + '\t' + str(interval[0]))
                continue
            else:
                if pos <= interval[1] and pos > interval[0]:
                    return True
                else:
                    return False
    except KeyError:
        return False   

def isNumber(testNum):
    try:
        float(testNum)
        return True
    except ValueError:
        return False

driverIDs = []
driverAnno = pd.read_csv('/seq/vgb/cancer_r01/ML_annotations/orthologousGeneSizesNonZeroGenesUpdatedToMatchGiantFile_091223.txt', sep = '\t')
driverAnno = driverAnno[driverAnno['driver'] == 'yes']
for i in list(driverAnno['binaryFeatureNames']):
    driverIDs.append(i.strip().split('_')[1])
    

brokenGenes = {}
nonDamageGenes = {}
noFlyList = []
with open(args.noFlyList, 'r') as inFile:
    for line in inFile:
        noFlyList.append(line.strip())

#try:
reference = buildRef(args.refFile, args.pathways)
#except:
#    print('Failed to load reference files, dying. Because I can\'t even.')
    
hg38CDS = {}
hg19CDS = {}
canFam3CDS = {}
hg19ToHg38 = '/seq/vgb/swofford/software/liftover/hg19ToHg38.over.chain'
canFam3ToHg38 = '/seq/vgb/swofford/software/liftover/canFam3ToHg38.over.chain'
outVariantSet = set()
outCodonSet = set()

rejectGenes = []
with open('/seq/vgb/swofford/R01/convertingToSKLInput/testVCFs/genesToRemoveFromParsing.txt', 'r') as inFile:
    for line in inFile:
        rejectGenes.append(line.strip())

with open(args.canfam3bed, 'r') as cf3, open(args.hg38bed, 'r') as hg38, open(args.hg19bed, 'r') as hg19:
    for line in cf3:
        fields = line.strip().split()
        contig = fields[0].strip()
        start = int(fields[1].strip())
        end = int(fields[2].strip())
        try:
            canFam3CDS[contig].append((start,end))
        except KeyError:
            canFam3CDS[contig] = [(start,end)]
    for i in canFam3CDS.keys():
        canFam3CDS[i].sort()
    for line in hg19:
        fields = line.strip().split()
        contig = fields[0].strip()
        start = int(fields[1].strip())
        end = int(fields[2].strip())
        try:
            hg19CDS[contig].append((start,end))
        except KeyError:
            hg19CDS[contig] = [(start,end)]
    for i in hg19CDS.keys():
        hg19CDS[i].sort()
    for line in hg38:
        fields = line.strip().split()
        contig = fields[0].strip()
        start = int(fields[1].strip())
        end = int(fields[2].strip())
        try:
            hg38CDS[contig].append((start,end))
        except KeyError:
            hg38CDS[contig] = [(start,end)]
    for i in hg38CDS.keys():
        hg38CDS[i].sort()
        
        
dogdb = gffutils.FeatureDB(args.dogDb)
hg19db = gffutils.FeatureDB(args.hg19Db)
hg38db = gffutils.FeatureDB(args.hg38Db)
samples = []
labels = set()
featureList = []
#genes = set()
#if args.maximumVariantsSample:
maxVar = int(args.maximumVariantsSample)
minVar = int(args.minimumVariantsSample)
#else:
#    maxVar = 1000000000

if args.granularity == 1:    
    studyClass = {
        'ALL-US':'ALL',
        'B-cell acute lymphoblastic leukemia, non-hypodiploid':'ALL',
        'B-cell acute lymphoblastic leukemia, hypodiploid':'ALL',
        'T-cell acute lymphoblastic leukemia':'ALL',
        'AML-US':'Remove',
        'LAML-CN':'Remove',
        'LAML-KR':'AML',
        'Acute myeloid leukemias':'AML',
        'LAML':'AML',
        'ACC':'ACC',
        'Adrenocortical carcinoma':'ACC',
        'BLCA-CN':'BLCA',
        'BLCA':'BLCA',
        'BOCA-FR':'Ewing',
        'BOCA-UK':'Remove',
        'Ewing\'s sarcoma':'Ewing',
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
        'Soft Tissue Sarcoma':'AS',
        'Breast Cancer': 'MBCP',
        'Prostate Cancer': 'MPCP'
        }
elif args.granularity == 2:
    studyClass= {
        'ALL-US':'Remove',
        'B-cell acute lymphoblastic leukemia, non-hypodiploid':'BALL',
        'B-cell acute lymphoblastic leukemia, hypodiploid':'BALL',
        'T-cell acute lymphoblastic leukemia':'TALL',
        'AML-US':'Remove',
        'LAML-CN':'Remove',
        'LAML-KR':'AML',
        'Acute myeloid leukemias':'AML',
        'LAML':'AML',
        'ACC':'ACC',
        'Adrenocortical carcinoma':'ACC',
        'BLCA-CN':'BLCA',
        'BLCA':'BLCA',
        'BOCA-FR':'Ewing',
        'BOCA-UK':'Remove',
        'Ewing\'s sarcoma':'Ewing',
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
        'Soft Tissue Sarcoma':'AS',
        'Breast Cancer': 'MBCP',
        'Prostate Cancer': 'MPCP'
        }
elif args.granularity == 3:
    studyClass= {
        'ALL-US':'ALL',
        'Soft Tissue Sarcoma':'AS',
        'B-cell acute lymphoblastic leukemia, non-hypodiploid':'ALL',
        'B-cell acute lymphoblastic leukemia, hypodiploid':'ALL',
        'T-cell acute lymphoblastic leukemia':'ALL',
        'AML-US':'Remove',
        'LAML-CN':'Remove',
        'LAML-KR':'AML',
        'Acute myeloid leukemias':'AML',
        'LAML':'AML',
        'ACC':'ACC',
        'Adrenocortical carcinoma':'ACC',
        'BLCA-CN':'BLCA',
        'BLCA':'BLCA',
        'BOCA-FR':'Ewing',
        'BOCA-UK':'Remove',
        'Ewing\'s sarcoma':'Ewing',
        'Osteosarcoma':'OSA',
        'Bone Cancer':'OSA',
        'BPLL-FR':'BPLL',
        'BRCA-EU':'BRCA',
        'BRCA-FR':'BRCA',
        'BRCA-KR':'BRCA',
        'BRCA-UK':'BRCA',
        'BRCA':'BRCA',
        'MBCP':'BRCA',
        'Breast Cancer':'BRCA',
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
        'Wilms\` tumors':'WT',
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
        'Burkitt\`s lymphoma':'LSA',
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
        'MPCP':'PRAD',
        'Prostate Cancer':'PRAD',
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
        'Wilms\' tumors':'Remove',
        'Burkitt\'s lymphoma':'Remove',
        'Soft Tissue Sarcoma':'AS',
        'Bone Cancer':'OSA'
        
        }  
elif args.granularity == 4:
    studyClass = {
        'ALL-US':'Remove',
        'B-cell acute lymphoblastic leukemia, non-hypodiploid':'BALL',
        'B-cell acute lymphoblastic leukemia, hypodiploid':'BALL',
        'T-cell acute lymphoblastic leukemia':'TALL',
        'AML-US':'Remove',
        'LAML-CN':'Remove',
        'LAML-KR':'AML',
        'Acute myeloid leukemias':'AML',
        'LAML':'AML',
        'ACC':'ACC',
        'Adrenocortical carcinoma':'ACC',
        'AS':'AS',
        'BLCA-CN':'BLCA',
        'BLCA':'BLCA',
        'BOCA-FR':'Ewings',
        'BOCA-UK':'OSA',
        'Ewing\'s sarcoma':'Ewings',
        'Osteosarcoma':'OSA',
        'BPLL-FR':'BPLL',
        'BRCA-EU':'BRCA',
        'BRCA-FR':'BRCA',
        'BRCA-KR':'BRCA',
        'BRCA-UK':'BRCA',
        'BRCA':'BRCA',
        'MBCP':'Remove',
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
        'Ependymoma supratentorial':'Ependymoma',
        'Ependymoma infratentorial':'Ependymoma',
        'ESAD-UK':'ESAD',
        'ESCA-CN':'ESCA',
        'ESCA':'ESAD',
        'ESCA':'ESCA',
        'GACA-CN':'Remove',
        'GACA-JP':'Remove',
        'STAD':'STAD',
        'GBM':'HGG',
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
        'SARC':'LMS',
        'LUAD':'LUAD',
        'LUSC':'LUSC',
        'LUSC-CN':'LUSC',
        'LUSC-KR':'Remove',
        'DLBC':'DLBCL',
        'MALY-DE':'Remove',
        'Burkitt\'s lymphoma':'BULY',
        'NKTL-SG':'Remove',
        'Medulloblastoma WNT':'PEME',
        'Medulloblastoma SHH':'PEME',
        'Medulloblastoma Group3':'PEME',
        'Medulloblastoma Group4':'PEME',
        'PEME-CA':'PEME',
        'PBCA-DE':'Remove',
        'SKCM':'MELA',
        'MELA-AU':'MELA',
        'SKCA-BR':'MELA',
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
        'PAEN-IT':'Remove',
        'PBCA-US':'Remove',
        'PCPG':'PCPG',
        'PBCA-DE':'Remove',
        'Pilocytic astrocytoma':'PAST',
        'EOPC-DE':'Remove',
        'PRAD-CA':'PRAD',
        'PRAD-CN':'Remove',
        'PRAD-FR':'PRAD',
        'PRAD-UK':'PRAD',
        'PRAD':'PRAD',
        'MPCP':'Remove',
        'Retinoblastoma':'RB',
        'RT-US':'RT',
        'Atypial teratoid/rhabdoid tumor':'ATRT',
        'Rhabdomyosarcoma':'RMS',
        'SARC':'Remove',
        'TGCT':'TGCT',
        'THYM':'THYM',
        'THCA-CN':'THCA',
        'THCA-SA':'THCA',
        'THCA':'THCA',
        'UTCA-FR':'UCS',
        'UCS':'UCS',
        'UCEC':'UCEC',
        'Prostate Cancer':'Remove',
        'Breast Cancer':'Remove',
        'Bone Cancer':'Remove',
        'Soft Tissue Sarcoma':'AS',
        'DLBCLNOS':'DLBCL',
        'Mature B-Cell Neoplasms':'DLBCL'
        }
elif int(args.granularity) > 4 or int(args.granularity) < 1:
    print('Level of granularity selected is out of bounds, I die now.')
    quit()
'''Old Angio list: keeping for reference but should be duplicated in new samples
    despite different IDs
    'Angio-ASCProject_03IpIBIR-Tumor-SM-DAE1F',
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
    'Angio-ASCProject_ZXi6FXFo-Tumor-SM-DADXW',
'''
    
'''Hand-curated list of AS samples to keep in analysis'''
ASKeep = [
    'RP-1447_ASCProject_03IpIBIR_T1_v2_Exome',
    'RP-1447_ASCProject_09iyC6Cp_T2_v3_Exome',
    'RP-1447_ASCProject_0RU8CnuV_T1_v4_Exome',
    'RP-1447_ASCProject_0xI6UVSy_T2_v3_Exome',
    'RP-1447_ASCProject_3NflfGHo_T1_v3_Exome',
    'RP-1447_ASCProject_44uvIDfp_BLOOD_P_v1_Exome',
    'RP-1447_ASCProject_47SztnTm_T1_v3_Exome',
    'RP-1447_ASCProject_5DFyF0HA_T1_v3_Exome',
    'RP-1447_ASCProject_5NsKslTb_T1_v3_Exome',
    'RP-1447_ASCProject_5Ph4hetv_T1_v3_Exome',
    'RP-1447_ASCProject_64i5tauQ_T2_v2_Exome',
    'RP-1447_ASCProject_7bt2t8IJ_T2_v2_Exome',
    'RP-1447_ASCProject_87uAC2t8_T1_v2_Exome',
    'RP-1447_ASCProject_A8UvuRTN_BLOOD_P_v1_Exome',
    'RP-1447_ASCProject_aPskieIL_T1_v3_Exome',
    'RP-1447_ASCProject_B0FQhzI9_T1_v4_Exome',
    'RP-1447_ASCProject_bJi9iyCb_T1_v3_Exome',
    'RP-1447_ASCProject_blh1uyhK_T1_v1_Exome',
    'RP-1447_ASCProject_bRSGSlCG_T2_v1_Exome',
    'RP-1447_ASCProject_dvH7FdsQ_T2_v3_Exome',
    'RP-1447_ASCProject_dyhLT8sG_T1_v2_Exome',
    'RP-1447_ASCProject_EAuOu4uD_T1_v1_Exome',
    'RP-1447_ASCProject_EdCVC0um_T1_v2_Exome',
    'RP-1447_ASCProject_eVIAsEIy_T1_v1_Exome',
    'RP-1447_ASCProject_EzU1Tmcp_T1A_v4_Exome',
    'RP-1447_ASCProject_GXTMTxU7_T1_v3_Exome',
    'RP-1447_ASCProject_JkivIBSV_T1_v1_Exome',
    'RP-1447_ASCProject_JltRf4Sd_T1_v2_Exome',
    'RP-1447_ASCProject_JvCxfgU1_T1_v2_Exome',
    'RP-1447_ASCProject_jWhlFLi3_BLOOD_P_v1_Exome',
    'RP-1447_ASCProject_Jyiwc4t6_T2A_v4_Exome',
    'RP-1447_ASCProject_kJFaIXcN_T2_v2_Exome',
    'RP-1447_ASCProject_kWCRuOFB_T1_v3_Exome',
    'RP-1447_ASCProject_KxFGsofW_T2_v3_Exome',
    'RP-1447_ASCProject_loC8UQUk_T1_v3_Exome',
    'RP-1447_ASCProject_LRS0tkCV_BLOOD_P_v1_Exome',
    'RP-1447_ASCProject_lYIMuWu5_T1_v1_Exome',
    'RP-1447_ASCProject_mpUNibFx_T1_v1_Exome',
    'RP-1447_ASCProject_NdUxUwCM_T1_v2_Exome',
    'RP-1447_ASCProject_nvfdfZCz_T1_v2_Exome',
    'RP-1447_ASCProject_oQi7i3U5_T2_v1_Exome',
    'RP-1447_ASCProject_OwHXUEs2_BLOOD_P_v1_Exome',
    'RP-1447_ASCProject_QYsAsrTO_T1_v2_Exome',
    'RP-1447_ASCProject_rnhwT3Hk_T1_v4_Exome',
    'RP-1447_ASCProject_RRFGCkUD_BLOOD_P_v1_Exome',
    'RP-1447_ASCProject_RvtOtjtj_T3_v2_Exome',
    'RP-1447_ASCProject_XbieiohK_T1_v1_Exome',
    'RP-1447_ASCProject_xdigi2t7_T1_v1_Exome',
    'RP-1447_ASCProject_xEFYtWU3_T1_v2_Exome',
    'RP-1447_ASCProject_XjHAU9uR_T2_v3_Exome',
    'RP-1532_ASCProject_XOHgH4CR_T1_v1_Exome',
    'RP-1447_ASCProject_Y4SjSnIg_T1_v1_Exome',
    'RP-1447_ASCProject_y5TafLf8_BLOOD_P_v1_Exome',
    'RP-1447_ASCProject_ybi8CAUN_T1_v4_Exome',
    'RP-1447_ASCProject_YLCRINC8_T1_v3_Exome',
    'RP-1447_ASCProject_YqsAUeFb_T1_v2_Exome',
    'RP-1447_ASCProject_YzTlTytj_T1_v1_Exome',
    'RP-1447_ASCProject_YzTxcouV_T2_v2_Exome',
    'RP-1447_ASCProject_ZbSaSACL_T1_v2_Exome',
    'RP-1447_ASCProject_zgCLCGco_T1_v1_Exome',
    'RP-1447_ASCProject_zQtVt5tr_BLOOD_v1_Exome',
    'RP-1447_ASCProject_ZXi6FXFo_T1_v1_Exome'
    ]

canineGliomaKeep = ['i_E7AB',
                    'i_51A5',
                    'i_F8E7',
                    'i_1166',
                    'i_8743',
                    'i_2EC9',
                    'i_3F8C',
                    'i_B023',
                    'i_0FF0',
                    'i_D756',
                    'i_22C7',
                    'i_92AC',
                    'i_BF76',
                    'i_E271',
                    'i_D026',
                    'i_607E',
                    'i_4D0C',
                    'i_C3C0',
                    'i_5CE5',
                    'i_157E',
                    'i_42D9',
                    'i_B02B',
                    'i_6561',
                    'i_4990',
                    'i_A71E',
                    'i_99AF',
                    'i_B2DC',
                    'i_DCD0',
                    'i_AB3E',
                    'i_FECA',
                    'i_D7EC',
                    'i_2C25',
                    'i_81E2',
                    'i_4BAC',
                    'i_A974',
                    'i_8A0A',
                    'i_5E9A',
                    'i_6E45',
                    'i_FC65',
                    'i_2C4F',
                    'i_03A6',
                    'i_6638',
                    'i_E952',
                    'i_C04D',
                    'i_750B',
                    'i_4AAB',
                    'i_B70F',
                    'i_F840',
                    'i_63FE',
                    'i_3688',
                    'i_B3CE',
                    'i_D0EE',
                    'i_66A9',
                    'i_6D5C',
                    'i_49E6',
                    'i_1793',
                    'i_1165',
                    'i_E2CD',
                    'i_05CA',
                    'i_6454',
                    'i_C561',
                    'i_B4F5',
                    'i_ED99',
                    'i_6254',
                    'i_8228',
                    'i_B813'
                    ]
cmtLookup = ['CMT_TN_26_final_ponmin1_nlod_0_0_PASS_IDOG_PON_EX.vcf.snpeff.vcf',
              'CMT_TN_29_final_ponmin1_nlod_0_0_PASS_IDOG_PON_EX.vcf.snpeff.vcf',
                'CMT_TN_31_final_ponmin1_nlod_0_0_PASS_IDOG_PON_EX.vcf.snpeff.vcf',
                'CMT_TN_33_final_ponmin1_nlod_0_0_PASS_IDOG_PON_EX.vcf.snpeff.vcf',
                'CMT_TN_32_final_ponmin1_nlod_0_0_PASS_IDOG_PON_EX.vcf.snpeff.vcf',
                'CMT_TN_9_final_ponmin1_nlod_0_0_PASS_IDOG_PON_EX.vcf.snpeff.vcf',
                'CMT_TN_16_final_ponmin1_nlod_0_0_PASS_IDOG_PON_EX.vcf.snpeff.vcf',
                'CMT_TN_17_final_ponmin1_nlod_0_0_PASS_IDOG_PON_EX.vcf.snpeff.vcf',
                'CMT_TN_78_final_ponmin1_nlod_0_0_PASS_IDOG_PON_EX.vcf.snpeff.vcf',
                'CMT_TN_23_final_ponmin1_nlod_0_0_PASS_IDOG_PON_EX.vcf.snpeff.vcf',
                'CMT_TN_24_final_ponmin1_nlod_0_0_PASS_IDOG_PON_EX.vcf.snpeff.vcf',
                'CMT_TN_51_final_ponmin1_nlod_0_0_PASS_IDOG_PON_EX.vcf.snpeff.vcf',
                'CMT_TN_7_final_ponmin1_nlod_0_0_PASS_IDOG_PON_EX.vcf.snpeff.vcf',
                'CMT_TN_13_final_ponmin1_nlod_0_0_PASS_IDOG_PON_EX.vcf.snpeff.vcf',
                'CMT_TN_36_final_ponmin1_nlod_0_0_PASS_IDOG_PON_EX.vcf.snpeff.vcf',
                'CMT_TN_38_final_ponmin1_nlod_0_0_PASS_IDOG_PON_EX.vcf.snpeff.vcf',
                'CMT_TN_39_final_ponmin1_nlod_0_0_PASS_IDOG_PON_EX.vcf.snpeff.vcf',
                'CMT_TN_40_final_ponmin1_nlod_0_0_PASS_IDOG_PON_EX.vcf.snpeff.vcf',
                'CMT_TN_43_final_ponmin1_nlod_0_0_PASS_IDOG_PON_EX.vcf.snpeff.vcf',
                'CMT_TN_54_final_ponmin1_nlod_0_0_PASS_IDOG_PON_EX.vcf.snpeff.vcf',
                'CMT_TN_55_final_ponmin1_nlod_0_0_PASS_IDOG_PON_EX.vcf.snpeff.vcf',
                'CMT_TN_59_final_ponmin1_nlod_0_0_PASS_IDOG_PON_EX.vcf.snpeff.vcf',
                'CMT_TN_75_final_ponmin1_nlod_0_0_PASS_IDOG_PON_EX.vcf.snpeff.vcf']
'''
cmtLookup = {'CMT.mutect2.somatic.CMT-002-tumor.snpeff.vcf':'dog_CMT_Carcinoma',
            'CMT.mutect2.somatic.CMT-003-tumor.snpeff.vcf':'dog_CMT_Carcinoma',
            'CMT.mutect2.somatic.CMT-013-tumor.snpeff.vcf':'dog_CMT_Carcinoma',
            'CMT.mutect2.somatic.CMT-021-tumor.snpeff.vcf':'dog_CMT_Carcinoma',
            'CMT.mutect2.somatic.CMT-026-tumor.snpeff.vcf':'dog_CMT_Carcinoma',
            'CMT.mutect2.somatic.CMT-030-tumor.snpeff.vcf':'dog_CMT_Carcinoma',
            'CMT.mutect2.somatic.CMT-043-tumor.snpeff.vcf':'dog_CMT_Carcinoma',
            'CMT.mutect2.somatic.CMT-048-tumor.snpeff.vcf':'dog_CMT_Carcinoma',
            'CMT.mutect2.somatic.CMT-053-tumor.snpeff.vcf':'dog_CMT_Carcinoma',
            'CMT.mutect2.somatic.CMT-054-tumor.snpeff.vcf':'dog_CMT_Carcinoma',
            'CMT.mutect2.somatic.CMT-057-tumor.snpeff.vcf':'dog_CMT_Carcinoma',
            'CMT.mutect2.somatic.CMT-059-tumor.snpeff.vcf':'dog_CMT_Carcinoma',
            'CMT.mutect2.somatic.CMT-063-tumor.snpeff.vcf':'dog_CMT_Carcinoma',
            'CMT.mutect2.somatic.CMT-072-tumor.snpeff.vcf':'dog_CMT_Carcinoma',
            'CMT.mutect2.somatic.CMT-077-tumor.snpeff.vcf':'dog_CMT_Carcinoma',
            'CMT.mutect2.somatic.CMT-078-tumor.snpeff.vcf':'dog_CMT_Carcinoma',
            'CMT.mutect2.somatic.CMT-084-tumor.snpeff.vcf':'dog_CMT_Carcinoma',
            'CMT.mutect2.somatic.CMT-099-tumor.snpeff.vcf':'dog_CMT_Carcinoma',
            'CMT.mutect2.somatic.CMT-100-tumor.snpeff.vcf':'dog_CMT_Carcinoma',
            'CMT.mutect2.somatic.CMT-105-tumor.snpeff.vcf':'dog_CMT_Carcinoma',
            'CMT.mutect2.somatic.CMT-107-tumor.snpeff.vcf':'dog_CMT_Carcinoma',
            'CMT.mutect2.somatic.CMT-109-tumor.snpeff.vcf':'dog_CMT_Carcinoma',
            'CMT.mutect2.somatic.CMT-111-tumor.snpeff.vcf':'dog_CMT_Carcinoma',
            'CMT.mutect2.somatic.CMT-118-tumor.snpeff.vcf':'dog_CMT_Carcinoma',
            'CMT.mutect2.somatic.CMT-119-tumor.snpeff.vcf':'dog_CMT_Carcinoma',
            'CMT.mutect2.somatic.CMT-121-tumor.snpeff.vcf':'dog_CMT_Carcinoma',
            'CMT.mutect2.somatic.CMT-129-tumor.snpeff.vcf':'dog_CMT_Carcinoma',
            'CMT.mutect2.somatic.CMT-142-tumor.snpeff.vcf':'dog_CMT_Carcinoma',
            'CMT.mutect2.somatic.CMT-147-tumor.snpeff.vcf':'dog_CMT_Carcinoma',
            'CMT.mutect2.somatic.CMT-150-tumor.snpeff.vcf':'dog_CMT_Carcinoma',
            'CMT.mutect2.somatic.CMT-159-tumor.snpeff.vcf':'dog_CMT_Carcinoma',
            'CMT.mutect2.somatic.CMT-168-tumor.snpeff.vcf':'dog_CMT_Carcinoma',
            'CMT.mutect2.somatic.CMT-169-tumor.snpeff.vcf':'dog_CMT_Carcinoma',
            'CMT.mutect2.somatic.CMT-183-tumor.snpeff.vcf':'dog_CMT_Carcinoma',
            'CMT.mutect2.somatic.CMT-184-tumor.snpeff.vcf':'dog_CMT_Carcinoma',
            'CMT.mutect2.somatic.CMT-187-tumor.snpeff.vcf':'dog_CMT_Carcinoma',
            'CMT.mutect2.somatic.CMT-188-tumor.snpeff.vcf':'dog_CMT_Carcinoma',
            'CMT.mutect2.somatic.CMT-190-tumor.snpeff.vcf':'dog_CMT_Carcinoma',
            'CMT.mutect2.somatic.CMT-193-tumor.snpeff.vcf':'dog_CMT_Carcinoma',
            'CMT.mutect2.somatic.CMT-197-tumor.snpeff.vcf':'dog_CMT_Carcinoma',
            'CMT.mutect2.somatic.CMT-201-tumor.snpeff.vcf':'dog_CMT_Carcinoma',
            'CMT.mutect2.somatic.CMT-206-tumor.snpeff.vcf':'dog_CMT_Carcinoma',
            'CMT.mutect2.somatic.CMT-210-tumor.snpeff.vcf':'dog_CMT_Carcinoma',
            'CMT.mutect2.somatic.CMT-217-tumor.snpeff.vcf':'dog_CMT_Carcinoma',
            'CMT.mutect2.somatic.CMT-227-tumor.snpeff.vcf':'dog_CMT_Carcinoma',
            'CMT.mutect2.somatic.CMT-231-tumor.snpeff.vcf':'dog_CMT_Carcinoma',
            'CMT.mutect2.somatic.CMT-233-tumor.snpeff.vcf':'dog_CMT_Carcinoma',
            'CMT.mutect2.somatic.CMT-234-tumor.snpeff.vcf':'dog_CMT_Carcinoma',
            'CMT.mutect2.somatic.CMT-236-tumor.snpeff.vcf':'dog_CMT_Carcinoma',
            'CMT.mutect2.somatic.CMT-239-tumor.snpeff.vcf':'dog_CMT_Carcinoma',
            'CMT.mutect2.somatic.CMT-241-tumor.snpeff.vcf':'dog_CMT_Carcinoma',
            'CMT.mutect2.somatic.CMT-242-tumor.snpeff.vcf':'dog_CMT_Carcinoma',
            'CMT.mutect2.somatic.CMT-246-tumor.snpeff.vcf':'dog_CMT_Carcinoma',
            'CMT.mutect2.somatic.CMT-247-tumor.snpeff.vcf':'dog_CMT_Carcinoma',
            'CMT.mutect2.somatic.CMT-250-tumor.snpeff.vcf':'dog_CMT_Carcinoma',
            'CMT.mutect2.somatic.CMT-252-tumor.snpeff.vcf':'dog_CMT_Carcinoma',
            'CMT.mutect2.somatic.CMT-259-tumor.snpeff.vcf':'dog_CMT_Carcinoma',
            'CMT.mutect2.somatic.CMT-262-tumor.snpeff.vcf':'dog_CMT_Carcinoma',
            'CMT.mutect2.somatic.CMT-263-tumor.snpeff.vcf':'dog_CMT_Carcinoma',
            'CMT.mutect2.somatic.CMT-271-tumor.snpeff.vcf':'dog_CMT_Carcinoma',
            'CMT.mutect2.somatic.CMT-278-tumor.snpeff.vcf':'dog_CMT_Carcinoma',
            'CMT.mutect2.somatic.CMT-283-tumor.snpeff.vcf':'dog_CMT_Carcinoma',
            'CMT.mutect2.somatic.CMT-292-tumor.snpeff.vcf':'dog_CMT_Carcinoma',
            'CMT.mutect2.somatic.CMT-301-tumor.snpeff.vcf':'dog_CMT_Carcinoma',
            'CMT.mutect2.somatic.CMT-306-tumor.snpeff.vcf':'dog_CMT_Carcinoma',
            'CMT.mutect2.somatic.CMT-307-tumor.snpeff.vcf':'dog_CMT_Carcinoma',
            'CMT.mutect2.somatic.CMT-308-tumor.snpeff.vcf':'dog_CMT_Carcinoma',
            'CMT.mutect2.somatic.CMT-318-tumor.snpeff.vcf':'dog_CMT_Carcinoma',
            'CMT.mutect2.somatic.CMT-319-tumor.snpeff.vcf':'dog_CMT_Carcinoma',
            'CMT.mutect2.somatic.CMT-321-tumor.snpeff.vcf':'dog_CMT_Carcinoma',
            'CMT.mutect2.somatic.CMT-336-tumor.snpeff.vcf':'dog_CMT_Carcinoma',
            'CMT.mutect2.somatic.CMT-340-tumor.snpeff.vcf':'dog_CMT_Carcinoma',
            'CMT.mutect2.somatic.CMT-344-tumor.snpeff.vcf':'dog_CMT_Carcinoma',
            'CMT.mutect2.somatic.CMT-346-tumor.snpeff.vcf':'dog_CMT_Carcinoma',
            'CMT.mutect2.somatic.CMT-348-tumor.snpeff.vcf':'dog_CMT_Carcinoma',
            'CMT.mutect2.somatic.CMT-353-tumor.snpeff.vcf':'dog_CMT_Carcinoma',
            'CMT.mutect2.somatic.CMT-354-tumor.snpeff.vcf':'dog_CMT_Carcinoma',
            'CMT.mutect2.somatic.CMT-360-tumor.snpeff.vcf':'dog_CMT_Carcinoma',
            'CMT.mutect2.somatic.CMT-365-tumor.snpeff.vcf':'dog_CMT_Carcinoma',
            'CMT.mutect2.somatic.CMT-369-tumor.snpeff.vcf':'dog_CMT_Carcinoma',
            'CMT.mutect2.somatic.CMT-377-tumor.snpeff.vcf':'dog_CMT_Carcinoma',
            'CMT.mutect2.somatic.CMT-381-tumor.snpeff.vcf':'dog_CMT_Carcinoma',
            'CMT.mutect2.somatic.CMT-400-tumor.snpeff.vcf':'dog_CMT_Carcinoma',
            'CMT.mutect2.somatic.CMT-402-tumor.snpeff.vcf':'dog_CMT_Carcinoma',
            'CMT.mutect2.somatic.CMT-414-tumor.snpeff.vcf':'dog_CMT_Carcinoma',
            'CMT.mutect2.somatic.CMT-435-tumor.snpeff.vcf':'dog_CMT_Carcinoma',
            'CMT.mutect2.somatic.CMT-436-tumor.snpeff.vcf':'dog_CMT_Carcinoma',
            'CMT.mutect2.somatic.CMT-439-tumor.snpeff.vcf':'dog_CMT_Carcinoma',
            'CMT.mutect2.somatic.CMT-440-tumor.snpeff.vcf':'dog_CMT_Carcinoma',
            'CMT.mutect2.somatic.CMT-458-tumor.snpeff.vcf':'dog_CMT_Carcinoma',
            'CMT.mutect2.somatic.CMT-459-tumor.snpeff.vcf':'dog_CMT_Carcinoma',
            'CMT.mutect2.somatic.CMT-471-tumor.snpeff.vcf':'dog_CMT_Carcinoma',
            'CMT.mutect2.somatic.CMT-495-tumor.snpeff.vcf':'dog_CMT_Carcinoma',
            'CMT.mutect2.somatic.CMT-496-tumor.snpeff.vcf':'dog_CMT_Carcinoma',
            'CMT.mutect2.somatic.CMT-505-tumor.snpeff.vcf':'dog_CMT_Carcinoma',
            'CMT.mutect2.somatic.CMT-523-tumor.snpeff.vcf':'dog_CMT_Carcinoma',
            'CMT.mutect2.somatic.CMT-538-tumor.snpeff.vcf':'dog_CMT_Carcinoma',
            'CMT.mutect2.somatic.CMT-577-tumor.snpeff.vcf':'dog_CMT_Carcinoma',
            'CMT.mutect2.somatic.CMT-581-tumor.snpeff.vcf':'dog_CMT_Carcinoma',
            'CMT.mutect2.somatic.CMT-586-tumor.snpeff.vcf':'dog_CMT_Carcinoma',
            'CMT.mutect2.somatic.CMT-603-tumor.snpeff.vcf':'dog_CMT_Carcinoma',
            'CMT.mutect2.somatic.CMT-606-tumor.snpeff.vcf':'dog_CMT_Carcinoma',
            'CMT.mutect2.somatic.CMT-613-tumor.snpeff.vcf':'dog_CMT_Carcinoma',
            'CMT.mutect2.somatic.CMT-615-tumor.snpeff.vcf':'dog_CMT_Carcinoma',
            'CMT.mutect2.somatic.CMT-630-tumor.snpeff.vcf':'dog_CMT_Carcinoma',
            'CMT.mutect2.somatic.CMT-632-tumor.snpeff.vcf':'dog_CMT_Carcinoma',
            'CMT.mutect2.somatic.CMT-634-tumor.snpeff.vcf':'dog_CMT_Carcinoma',
            'CMT.mutect2.somatic.CMT-638-tumor.snpeff.vcf':'dog_CMT_Carcinoma',
            'CMT.mutect2.somatic.CMT-642-tumor.snpeff.vcf':'dog_CMT_Carcinoma',
            'CMT.mutect2.somatic.CMT-644-tumor.snpeff.vcf':'dog_CMT_Carcinoma',
            'CMT.mutect2.somatic.CMT-651-tumor.snpeff.vcf':'dog_CMT_Carcinoma',
            'CMT.mutect2.somatic.CMT-652-tumor.snpeff.vcf':'dog_CMT_Carcinoma',
            'CMT.mutect2.somatic.CMT-659-tumor.snpeff.vcf':'dog_CMT_Carcinoma',
            'CMT.mutect2.somatic.CMT-660-tumor.snpeff.vcf':'dog_CMT_Carcinoma',
            'CMT.mutect2.somatic.CMT-667-tumor.snpeff.vcf':'dog_CMT_Carcinoma',
            'CMT.mutect2.somatic.CMT-678-tumor.snpeff.vcf':'dog_CMT_Carcinoma',
            'CMT.mutect2.somatic.CMT-685-tumor.snpeff.vcf':'dog_CMT_Carcinoma',
            'CMT.mutect2.somatic.CMT-686-tumor.snpeff.vcf':'dog_CMT_Carcinoma',
            'CMT.mutect2.somatic.CMT-688-tumor.snpeff.vcf':'dog_CMT_Carcinoma',
            'CMT.mutect2.somatic.CMT-690-tumor.snpeff.vcf':'dog_CMT_Carcinoma',
            'CMT.mutect2.somatic.CMT-724-tumor.snpeff.vcf':'dog_CMT_Carcinoma',
            'CMT.mutect2.somatic.CMT-725-tumor.snpeff.vcf':'dog_CMT_Carcinoma',
            'CMT.mutect2.somatic.CMT-730-tumor.snpeff.vcf':'dog_CMT_Carcinoma',
            'CMT.mutect2.somatic.CMT-756-tumor.snpeff.vcf':'dog_CMT_Carcinoma',
            'CMT.mutect2.somatic.CMT-774-tumor.snpeff.vcf':'dog_CMT_Carcinoma',
            'CMT.mutect2.somatic.CMT-785-tumor.snpeff.vcf':'dog_CMT_Carcinoma',
            'CMT.mutect2.somatic.CMT-785-tumor.snpeff.vcf':'dog_CMT_Benign',
            'CMT.mutect2.somatic.CMT-182-tumor.snpeff.vcf':'dog_CMT_Benign',
            'CMT.mutect2.somatic.CMT-203-tumor.snpeff.vcf':'dog_CMT_Benign',
            'CMT.mutect2.somatic.CMT-207-tumor.snpeff.vcf':'dog_CMT_Benign',
            'CMT.mutect2.somatic.CMT-215-tumor.snpeff.vcf':'dog_CMT_Benign',
            'CMT.mutect2.somatic.CMT-232-tumor.snpeff.vcf':'dog_CMT_Benign',
            'CMT.mutect2.somatic.CMT-245-tumor.snpeff.vcf':'dog_CMT_Benign',
            'CMT.mutect2.somatic.CMT-257-tumor.snpeff.vcf':'dog_CMT_Benign',
            'CMT.mutect2.somatic.CMT-272-tumor.snpeff.vcf':'dog_CMT_Benign',
            'CMT.mutect2.somatic.CMT-276-tumor.snpeff.vcf':'dog_CMT_Benign',
            'CMT.mutect2.somatic.CMT-311-tumor.snpeff.vcf':'dog_CMT_Benign',
            'CMT.mutect2.somatic.CMT-317-tumor.snpeff.vcf':'dog_CMT_Benign',
            'CMT.mutect2.somatic.CMT-355-tumor.snpeff.vcf':'dog_CMT_Benign',
            'CMT.mutect2.somatic.CMT-374-tumor.snpeff.vcf':'dog_CMT_Benign',
            'CMT.mutect2.somatic.CMT-392-tumor.snpeff.vcf':'dog_CMT_Benign',
            'CMT.mutect2.somatic.CMT-405-tumor.snpeff.vcf':'dog_CMT_Benign',
            'CMT.mutect2.somatic.CMT-407-tumor.snpeff.vcf':'dog_CMT_Benign',
            'CMT.mutect2.somatic.CMT-417-tumor.snpeff.vcf':'dog_CMT_Benign',
            'CMT.mutect2.somatic.CMT-442-tumor.snpeff.vcf':'dog_CMT_Benign',
            'CMT.mutect2.somatic.CMT-449-tumor.snpeff.vcf':'dog_CMT_Benign',
            'CMT.mutect2.somatic.CMT-452-tumor.snpeff.vcf':'dog_CMT_Benign',
            'CMT.mutect2.somatic.CMT-456-tumor.snpeff.vcf':'dog_CMT_Benign',
            'CMT.mutect2.somatic.CMT-480-tumor.snpeff.vcf':'dog_CMT_Benign',
            'CMT.mutect2.somatic.CMT-522-tumor.snpeff.vcf':'dog_CMT_Benign',
            'CMT.mutect2.somatic.CMT-695-tumor.snpeff.vcf':'dog_CMT_Benign',
            'CMT.mutect2.somatic.CMT-711-tumor.snpeff.vcf':'dog_CMT_Benign',
            'CMT.mutect2.somatic.CMT-719-tumor.snpeff.vcf':'dog_CMT_Benign',
            'CMT.mutect2.somatic.CMT-737-tumor.snpeff.vcf':'dog_CMT_Benign',
            'CMT.mutect2.somatic.CMT-776-tumor.snpeff.vcf':'dog_CMT_Benign'
            }
'''
pediatric_columbia_2019_keep = ['PIP14-71727-T1',
                                'PIP15-11925-T1',
                                'PIP15-35162-T1',
                                'PIP15-79265-T1',
                                'PIP15-83826-T1'
                                ]

metastatic_solid_tumors_mich_2017_keep = ['MO_1024',
                                          'MO_1525',
                                          'MO_1226',
                                          'MO_1555',
                                          'TP_2013'
                                          ]

brca_mbcproject_wagle_2017_keep = ['MBC-MBCProject_0aTBfPcG-Tumor-SM-DL3W1',
                                   'MBC-MBCProject_1ps2hZH8-Tumor-SM-DL4OO',
                                    'MBC-MBCProject_1qhlIasw-Tumor-SM-AXGNH',
                                    'MBC-MBCProject_27uAugT4-Tumor-SM-DL45T',
                                    'MBC-MBCProject_2kS5fYc1-Tumor-SM-DL4O2',
                                    'MBC-MBCProject_2pt4t9cA-Tumor-SM-GQD4Q',
                                    'MBC-MBCProject_2ViLiZS8-Tumor-SM-DL3Z6',
                                    'MBC-MBCProject_2WCGcxuk-Tumor-SM-GQAVI',
                                    'MBC-MBCProject_2Wf9hPtA-Tumor-SM-GQB5X',
                                    'MBC-MBCProject_3jhES9fq-Tumor-SM-AXGIU',
                                    'MBC-MBCProject_3LI1TNs3-Tumor-SM-DL43B',
                                    'MBC-MBCProject_3MUnT7Ty-Tumor-SM-CGLLE',
                                    'MBC-MBCProject_3vhkhAcY-Tumor-SM-GQAW5',
                                    'MBC-MBCProject_43UpSwhz-Tumor-SM-DL4W8',
                                    'MBC-MBCProject_4DIpSBFV-Tumor-SM-CGM2Q',
                                    'MBC-MBCProject_4MF1FlFQ-Tumor-SM-CGM4M',
                                    'MBC-MBCProject_4OIOfAt9-Tumor-SM-CGLYL',
                                    'MBC-MBCProject_4vC6Fgio-Tumor-SM-CGLZ8',
                                    'MBC-MBCProject_57iLiJIl-Tumor-SM-CGLIV',
                                    'MBC-MBCProject_5gHasou8-Tumor-SM-GQA1C',
                                    'MBC-MBCProject_5jSPs9fl-Tumor-SM-CGL1A',
                                    'MBC-MBCProject_5vTpUpCv-Tumor-SM-GQA4H',
                                    'MBC-MBCProject_6JhKimhv-Tumor-SM-GQCJC',
                                    'MBC-MBCProject_6QhZF6ur-Tumor-SM-GQD37',
                                    'MBC-MBCProject_6vTVHzur-Tumor-SM-GQCG7',
                                    'MBC-MBCProject_6Yt3czsW-Tumor-SM-GQC9A',
                                    'MBC-MBCProject_6zfRF2fl-Tumor-SM-GQ9UO',
                                    'MBC-MBCProject_70CxiQFk-Tumor-SM-DL3TI',
                                    'MBC-MBCProject_74CYIyHw-Tumor-SM-DL3P4',
                                    'MBC-MBCProject_7AC8CoUL-Tumor-SM-GQ9WK',
                                    'MBC-MBCProject_7JtRIPI2-Tumor-SM-CGLXC',
                                    'MBC-MBCProject_7OcOivcV-Tumor-SM-DL44K',
                                    'MBC-MBCProject_7of5i3hp-Tumor-SM-CGMAA',
                                    'MBC-MBCProject_7VFnIqfw-Tumor-SM-CGL6X',
                                    'MBC-MBCProject_7wCjtKIW-Tumor-SM-GQCN4',
                                    'MBC-MBCProject_7XFmtytw-Tumor-SM-DL457',
                                    'MBC-MBCProject_8bfWIQu1-Tumor-SM-GQC81',
                                    'MBC-MBCProject_99CdCOHm-Tumor-SM-CGLF4',
                                    'MBC-MBCProject_9JI1fwuX-Tumor-SM-DL4TP',
                                    'MBC-MBCProject_9RtxtGUO-Tumor-SM-GQC4V',
                                    'MBC-MBCProject_AJF3tVc8-Tumor-SM-GQAJK',
                                    'MBC-MBCProject_Ali6SAHL-Tumor-SM-GQAH2',
                                    'MBC-MBCProject_AqtYu0IW-Tumor-SM-DL4Y4',
                                    'MBC-MBCProject_AwukckTD-Tumor-SM-GQADA',
                                    'MBC-MBCProject_AzCnt3sL-Tumor-SM-GQAY1',
                                    'MBC-MBCProject_B7f5hRfD-Tumor-SM-GQACN',
                                    'MBC-MBCProject_bAi0Slsx-Tumor-SM-CGMFB',
                                    'MBC-MBCProject_beTYUZij-Tumor-SM-GQBWY',
                                    'MBC-MBCProject_BOHWtZT8-Tumor-SM-CGM6I',
                                    'MBC-MBCProject_bpunhmim-Tumor-SM-DL4MS',
                                    'MBC-MBCProject_briXiduE-Tumor-SM-DL49L',
                                    'MBC-MBCProject_brTys8Hl-Tumor-SM-GQCYF',
                                    'MBC-MBCProject_BVCVuGhl-Tumor-SM-DL4DD',
                                    'MBC-MBCProject_BVHET8iq-Tumor-SM-CGLGD',
                                    'MBC-MBCProject_bvHXCoTY-Tumor-SM-DL3JG',
                                    'MBC-MBCProject_BvsPhwI1-Tumor-SM-DL4A8',
                                    'MBC-MBCProject_d5CbUNTb-Tumor-SM-CGM5V',
                                    'MBC-MBCProject_DDUmIZuW-Tumor-SM-DL3GX',
                                    'MBC-MBCProject_DWTkcWCW-Tumor-SM-GQB7T',
                                    'MBC-MBCProject_DxibF0hM-Tumor-SM-CGLWP',
                                    'MBC-MBCProject_DximH6Cn-Tumor-SM-DL3XW',
                                    'MBC-MBCProject_e4SaTnIW-Tumor-SM-DL4KA',
                                    'MBC-MBCProject_EBTBIGFy-Tumor-SM-DL3ZS',
                                    'MBC-MBCProject_EkHAIECZ-Tumor-SM-GQ9OE',
                                    'MBC-MBCProject_ePi1smiM-Tumor-SM-CGLOK',
                                    'MBC-MBCProject_epUYsdCE-Tumor-SM-DL3A1',
                                    'MBC-MBCProject_ErfKfJt0-Tumor-SM-AZ5GM',
                                    'MBC-MBCProject_EWubi6hd-Tumor-SM-GQ9PN',
                                    'MBC-MBCProject_eXSVsxUm-Tumor-SM-GQC2Z',
                                    'MBC-MBCProject_gASyuQSl-Tumor-SM-GQAQH',
                                    'MBC-MBCProject_gdIlF4hG-Tumor-SM-DL3ML',
                                    'MBC-MBCProject_gjhMuoha-Tumor-SM-AXGGJ',
                                    'MBC-MBCProject_gKsRsqTR-Tumor-SM-AZ5J5',
                                    'MBC-MBCProject_GKu6TliY-Tumor-SM-DL4EM',
                                    'MBC-MBCProject_glf3CWFa-Tumor-SM-GQ9N5',
                                    'MBC-MBCProject_goSDCjUM-Tumor-SM-DL4BH',
                                    'MBC-MBCProject_gotjfgf0-Tumor-SM-AXGPS',
                                    'MBC-MBCProject_GouMi0U1-Tumor-SM-DL3KP',
                                    'MBC-MBCProject_GpizsQiW-Tumor-SM-GQA8V',
                                    'MBC-MBCProject_GvHkH2Hk-Tumor-SM-AZ5H9',
                                    'MBC-MBCProject_Gxhdtafa-Tumor-SM-GQ9U2',
                                    'MBC-MBCProject_j2F2sQC0-Tumor-SM-DL3CJ',
                                    'MBC-MBCProject_JBh5TjcM-Tumor-SM-DL3EF',
                                    'MBC-MBCProject_jEhBHrS2-Tumor-SM-CGKTZ',
                                    'MBC-MBCProject_Jeu1F5cz-Tumor-SM-DL48Y',
                                    'MBC-MBCProject_JGCBTpH1-Tumor-SM-CGLR2',
                                    'MBC-MBCProject_JGcLUJiA-Tumor-SM-GQBDH',
                                    'MBC-MBCProject_JKTZhEUe-Tumor-SM-DL48C',
                                    'MBC-MBCProject_jmfDfEs8-Tumor-SM-GQ9Z3',
                                    'MBC-MBCProject_JpCASlSG-Tumor-SM-GQC7E',
                                    'MBC-MBCProject_K0UDUnFE-Tumor-SM-AZ5FD',
                                    'MBC-MBCProject_k4hlFwFg-Tumor-SM-GQD1B',
                                    'MBC-MBCProject_K7f6fdUz-Tumor-SM-AZ5MA',
                                    'MBC-MBCProject_kAFQF3t6-Tumor-SM-DL36W',
                                    'MBC-MBCProject_kduys9h5-Tumor-SM-DL47P',
                                    'MBC-MBCProject_KGs1fxH2-Tumor-SM-CGM1H',
                                    'MBC-MBCProject_KphMt4F9-Tumor-SM-DL4PB',
                                    'MBC-MBCProject_kQTqIOSP-Tumor-SM-GQAWR',
                                    'MBC-MBCProject_kRC3tlcA-Tumor-SM-DL3LY',
                                    'MBC-MBCProject_kwS1CXSR-Tumor-SM-DL4SG',
                                    'MBC-MBCProject_kzuMSZIW-Tumor-SM-DL4JN',
                                    'MBC-MBCProject_LDCbC8t9-Tumor-SM-GQCEB',
                                    'MBC-MBCProject_lGCMIGT0-Tumor-SM-GQAM3',
                                    'MBC-MBCProject_lJflt3UR-Tumor-SM-GQAEJ',
                                    'MBC-MBCProject_LnHAS5T0-Tumor-SM-GQB93',
                                    'MBC-MBCProject_LPHKFauY-Tumor-SM-GQCCF',
                                    'MBC-MBCProject_lqSlSztO-Tumor-SM-GQB9P',
                                    'MBC-MBCProject_lQtMtjFR-Tumor-SM-GQCWJ',
                                    'MBC-MBCProject_LvS2IvIy-Tumor-SM-DL3IT',
                                    'MBC-MBCProject_LVSjf8h7-Tumor-SM-DL3U5',
                                    'MBC-MBCProject_lXCmTEuv-Tumor-SM-DL4WU',
                                    'MBC-MBCProject_m1fOs2FL-Tumor-SM-GQA54',
                                    'MBC-MBCProject_m9SNc1Iq-Tumor-SM-GQAIB',
                                    'MBC-MBCProject_MlheH1iY-Tumor-SM-GQD2K',
                                    'MBC-MBCProject_mMhdcrh5-Tumor-SM-GQ9YG',
                                    'MBC-MBCProject_MmSBTJtJ-Tumor-SM-GQAGF',
                                    'MBC-MBCProject_mrhKt1Ue-Tumor-SM-DL4UY',
                                    'MBC-MBCProject_mrIQIdsx-Tumor-SM-GQATM',
                                    'MBC-MBCProject_MYuoCmI7-Tumor-SM-DL4R7',
                                    'MBC-MBCProject_N4srsKsr-Tumor-SM-CGMGK',
                                    'MBC-MBCProject_NDIbs9i6-Tumor-SM-DL4VL',
                                    'MBC-MBCProject_nEcXsyfj-Tumor-SM-GQB5B',
                                    'MBC-MBCProject_nkU3hDhB-Tumor-SM-GQCKL',
                                    'MBC-MBCProject_nwIeUoSD-Tumor-SM-GQC2D',
                                    'MBC-MBCProject_NzHBsOtg-Tumor-SM-GQBL1',
                                    'MBC-MBCProject_nZHYc4Ie-Tumor-SM-DL39E',
                                    'MBC-MBCProject_O8u9SjuW-Tumor-SM-AZ5O6',
                                    'MBC-MBCProject_oNI6SXtq-Tumor-SM-GQC8N',
                                    'MBC-MBCProject_oqupfDu7-Tumor-SM-GQCOZ',
                                    'MBC-MBCProject_OwuoTVhg-Tumor-SM-DL3OH',
                                    'MBC-MBCProject_oZT8clio-Tumor-SM-GQCVW',
                                    'MBC-MBCProject_p1CQTdIg-Tumor-SM-CGL36',
                                    'MBC-MBCProject_P2CMUEiX-Tumor-SM-DL3GB',
                                    'MBC-MBCProject_PASlTyTy-Tumor-SM-DL473',
                                    'MBC-MBCProject_pktAIxFb-Tumor-SM-GQABE',
                                    'MBC-MBCProject_PkTDsOSa-Tumor-SM-GQ9X7',
                                    'MBC-MBCProject_pMcWcrtZ-Tumor-SM-AZ5KE',
                                    'MBC-MBCProject_PXcjFYIo-Tumor-SM-DL3UR',
                                    'MBC-MBCProject_pyhbI1H5-Tumor-SM-GQCX6',
                                    'MBC-MBCProject_PzUBH7I3-Tumor-SM-CGL52',
                                    'MBC-MBCProject_QbiZiytx-Tumor-SM-CGM9N',
                                    'MBC-MBCProject_QJFdf8hQ-Tumor-SM-GQCJY',
                                    'MBC-MBCProject_qmu6TYto-Tumor-SM-DL3BW',
                                    'MBC-MBCProject_QNs8uGhM-Tumor-SM-GQCB6',
                                    'MBC-MBCProject_qQSACYuv-Tumor-SM-GQBFZ',
                                    'MBC-MBCProject_RJf1h2Sp-Tumor-SM-DL4M6',
                                    'MBC-MBCProject_rJHBiKTl-Tumor-SM-DL3K3',
                                    'MBC-MBCProject_RKf1frsr-Tumor-SM-CGLJI',
                                    'MBC-MBCProject_rMHPsbIq-Tumor-SM-DL422',
                                    'MBC-MBCProject_rmSBIguN-Tumor-SM-DL3N8',
                                    'MBC-MBCProject_rvULI3TV-Tumor-SM-DL41F',
                                    'MBC-MBCProject_VbfASefz-Tumor-SM-CGLVG',
                                    'MBC-MBCProject_VQfmS0fK-Tumor-SM-CGKZN',
                                    'MBC-MBCProject_VrsMsqTb-Tumor-SM-CGLRO',
                                    'MBC-MBCProject_W2CQswiv-Tumor-SM-GQCPM',
                                    'MBC-MBCProject_W4FBsLSx-Tumor-SM-GQA6D',
                                    'MBC-MBCProject_W5ceiLi0-Tumor-SM-GQAR4',
                                    'MBC-MBCProject_W9hmcWcz-Tumor-SM-GQ9VB',
                                    'MBC-MBCProject_wAiri7fp-Tumor-SM-AZ5DH',
                                    'MBC-MBCProject_wKuwTYhK-Tumor-SM-GQCRI',
                                    'MBC-MBCProject_wKuZuQS1-Tumor-SM-CGLW3',
                                    'MBC-MBCProject_WYHpfXf1-Tumor-SM-GQAFS',
                                    'MBC-MBCProject_wzCxuoio-Tumor-SM-GQCHG',
                                    'MBC-MBCProject_x4cLHdhw-Tumor-SM-GQ9P1',
                                    'MBC-MBCProject_xBfJfri9-Tumor-SM-CGLAP',
                                    'MBC-MBCProject_xlhkS5CG-Tumor-SM-GQ9ZP',
                                    'MBC-MBCProject_XmiOUMSm-Tumor-SM-DL3SV',
                                    'MBC-MBCProject_XVHkUVTb-Tumor-SM-GQCLU',
                                    'MBC-MBCProject_y2F3IgIE-Tumor-SM-GQBAC',
                                    'MBC-MBCProject_Y7fYC1iG-Tumor-SM-CGLFQ',
                                    'MBC-MBCProject_Y7iLIOU0-Tumor-SM-GQBY8',
                                    'MBC-MBCProject_Y8H6SghO-Tumor-SM-DL4IE',
                                    'MBC-MBCProject_ygcMFgCR-Tumor-SM-GQ9S6',
                                    'MBC-MBCProject_YNTyfDfq-Tumor-SM-GQAK7',
                                    'MBC-MBCProject_yZSnSluK-Tumor-SM-DL4T3',
                                    'MBC-MBCProject_zatdFvhp-Tumor-SM-GQAAR',
                                    'MBC-MBCProject_ZdudUNFZ-Tumor-SM-CGLPS',
                                    'MBC-MBCProject_ZeTySaU3-Tumor-SM-GQD3T',
                                    'MBC-MBCProject_ZPHRs3Uw-Tumor-SM-DL3BA',
                                    'MBC-MBCProject_zyt5TKFB-Tumor-SM-CGM59'
                                    ]

prad_mpcproject_2018_keep = ['RP-1532_PCProject_0muduPhG_T1_v2_Exome_OnPrem',
                             'RP-1532_PCProject_2JIGIWhG_T1_v3_Exome_OnPrem',
                            'RP-1532_PCProject_4ds9Sbcj_T1_v2_Exome_OnPrem',
                            'RP-1532_PCProject_55fdtjcR_BLOOD_P_v1_WGS_OnPrem',
                            'RP-1532_PCProject_5Li7ikSy_T2_v2_Exome_OnPrem',
                            'RP-1532_PCProject_79FLFYI0_BLOOD_P_v1_WGS_OnPrem',
                            'RP-1532_PCProject_7VC9Cvs0_BLOOD_P_v1_WGS_OnPrem',
                            'RP-1532_PCProject_7XsWsJin_T1_v2_Exome_OnPrem',
                            'RP-1532_PCProject_ArI5sGCG_BLOOD_P_v1_WGS_OnPrem',
                            'RP-1532_PCProject_avuEc5H8_BLOOD_P_v1_WGS_OnPrem',
                            'RP-1532_PCProject_BGS1SnHl_BLOOD_P_v1_WGS_OnPrem',
                            'RP-1532_PCProject_d5uxuKij_BLOOD_P_v1_WGS_OnPrem',
                            'RP-1532_PCProject_D7hnfXHe_BLOOD_P_v1_WGS_OnPrem',
                            'RP-1532_PCProject_d9iyS4Cx_BLOOD_P_v1_WGS_OnPrem',
                            'RP-1532_PCProject_DXhLCjtb_BLOOD_P_v1_WGS_OnPrem',
                            'RP-1532_PCProject_e8hgh8Sg_T1_v3_Exome_OnPrem',
                            'RP-1066_PCProject_EDCDTqcM_T1_v1_Exome_OnPrem',
                            'RP-1532_PCProject_epuBfyTG_T1_v2_Exome_OnPrem',
                            'RP-1532_PCProject_gZS1fRSq_BLOOD_P_v1_WGS_OnPrem',
                            'RP-1532_PCProject_JEFQcbUQ_T1_v3_Exome_OnPrem',
                            'RP-1532_PCProject_JxHQUms9_BLOOD_P_v1_WGS_OnPrem',
                            'RP-1532_PCProject_kDCJcNi0_T2_v3_Exome_OnPrem',
                            'RP-1532_PCProject_LaiMu5tY_T1_v2_Exome_OnPrem',
                            'RP-1532_PCProject_mZfgfKHn_T1_v2_Exome_OnPrem',
                            'RP-1532_PCProject_njSMF6hZ_T2_v3_Exome_OnPrem',
                            'RP-1532_PCProject_nqfBfWTL_T1_v3_Exome_OnPrem',
                            'RP-1532_PCProject_olC9Uysp_T2_v2_Exome_OnPrem',
                            'RP-1532_PCProject_oOuWh8f4_BLOOD_P_v1_WGS_OnPrem',
                            'RP-1532_PCProject_oRH7HAFg_T1_v3_Exome_OnPrem',
                            'RP-1532_PCProject_P0uOSycP_BLOOD_P_v1_WGS_OnPrem',
                            'RP-1532_PCProject_p8HMUDIj_BLOOD_P_v1_WGS_OnPrem',
                            'RP-1532_PCProject_p8hwcKFl_T1_v2_Exome_OnPrem',
                            'RP-1532_PCProject_PeS4ckiA_BLOOD_P_v1_WGS_OnPrem',
                            'RP-1532_PCProject_PvfjhmtR_T1_v2_Exome_OnPrem',
                            'RP-1532_PCProject_PyHlUlcd_T1_v2_Exome_OnPrem',
                            'RP-1532_PCProject_QYIjibtE_T1_v1_Exome_OnPrem',
                            'RP-1532_PCProject_R5hOh7u3_BLOOD_P_v1_WGS_OnPrem',
                            'RP-1532_PCProject_rKtetVUV_BLOOD_P_v1_WGS_OnPrem',
                            'RP-1532_PCProject_rliEhmUJ_BLOOD_P_v1_WGS_OnPrem',
                            'RP-1532_PCProject_Rou7FZcR_T1_v2_Exome_OnPrem',
                            'RP-1532_PCProject_V9UeINs2_T1_v3_Exome_OnPrem',
                            'RP-1532_PCProject_vLiXiBtz_T1_v3_Exome_OnPrem',
                            'RP-1532_PCProject_voixIotM_BLOOD_P_v1_WGS_OnPrem',
                            'RP-1532_PCProject_VPS8SMCO_T1_v3_Exome_OnPrem',
                            'RP-1532_PCProject_x2uxtrue_BLOOD_P_v1_WGS_OnPrem',
                            'RP-1532_PCProject_XmUahbfN_T1_v2_Exome_OnPrem',
                            'RP-1532_PCProject_xpI9ixc9_T1B_v1_Exome_OnPrem',
                            'RP-1532_PCProject_Y3hQI7Ik_BLOOD_P_v1_WGS_OnPrem',
                            'RP-1532_PCProject_YbsksWCk_T1_v2_Exome_OnPrem',
                            'RP-1532_PCProject_ymtAUBUP_BLOOD_P_v1_WGS_OnPrem',
                            'RP-1532_PCProject_yXfgCQUn_T1_v3_Exome_OnPrem',
                            'RP-1532_PCProject_z2U1UeC3_T1_v3_Exome_OnPrem',
                            'RP-1532_PCProject_zjfaFjF8_T2_v3_Exome_OnPrem'
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
                   '3_prime_UTR_truncation + exon_loss','UTR_3_DELETED',
                   '5_prime_UTR_truncation + exon_loss_variant','UTR_5_DELETED'
                   ]
LSALookup = {
    'CCB010005_0202.snp.filtered.snpeff.vcf':'Remove',
    'CCB010009_200.snp.filtered.snpeff.vcf':'dogBLSA',
    'CCB010010_200.snp.filtered.snpeff.vcf':'dogBLSA',
    'CCB010023_0201.snp.filtered.snpeff.vcf':'Remove',
    'CCB010024_0200.snp.filtered.snpeff.vcf':'Remove',
    'CCB010026_0203.snp.filtered.snpeff.vcf':'Remove',
    'CCB010028_0200.snp.filtered.snpeff.vcf':'dogTLSA',
    'CCB010035_0200.snp.filtered.snpeff.vcf':'Remove',
    'CCB010036_200.snp.filtered.snpeff.vcf':'dogBLSA',
    'CCB010040_200.snp.filtered.snpeff.vcf':'dogBLSA',
    'CCB010054_0203.snp.filtered.snpeff.vcf':'dogBLSA',
    'CCB010059_0200.snp.filtered.snpeff.vcf':'Remove',
    'CCB010060_0201.snp.filtered.snpeff.vcf':'Remove',
    'CCB010069_0201.snp.filtered.snpeff.vcf':'Remove',
    'CCB010089_0203.snp.filtered.snpeff.vcf':'Remove',
    'CCB010099_0200.snp.filtered.snpeff.vcf':'dogBLSA',
    'CCB010102_0200.snp.filtered.snpeff.vcf':'Remove',
    'CCB010153_201.snp.filtered.snpeff.vcf':'dogBLSA',
    'CCB010174_0200.snp.filtered.snpeff.vcf':'Remove',
    'CCB010182_201.snp.filtered.snpeff.vcf':'dogBLSA',
    'CCB010214_200.snp.filtered.snpeff.vcf':'dogBLSA',
    'CCB010217_0200.snp.filtered.snpeff.vcf':'Remove',
    'CCB010233_201.snp.filtered.snpeff.vcf':'dogTLSA',
    'CCB010243_200.snp.filtered.snpeff.vcf':'dogTLSA',
    'CCB010253_201.snp.filtered.snpeff.vcf':'dogTLSA',
    'CCB010289_201.snp.filtered.snpeff.vcf':'dogTLSA',
    'CCB010299_0200.snp.filtered.snpeff.vcf':'Remove',
    'CCB010307_201.snp.filtered.snpeff.vcf':'dogBLSA',
    'CCB010327_0201.snp.filtered.snpeff.vcf':'Remove',
    'CCB010330_0201.snp.filtered.snpeff.vcf':'Remove',
    'CCB010332_0201.snp.filtered.snpeff.vcf':'Remove',
    'CCB010346_200.snp.filtered.snpeff.vcf':'dogBLSA',
    'CCB020009_200.snp.filtered.snpeff.vcf':'dogBLSA',
    'CCB020012_0207.snp.filtered.snpeff.vcf':'Remove',
    'CCB020022_0204.snp.filtered.snpeff.vcf':'Remove',
    'CCB020030_0200.snp.filtered.snpeff.vcf':'dogTLSA',
    'CCB020062_0200.snp.filtered.snpeff.vcf':'Remove',
    'CCB020071_0202.snp.filtered.snpeff.vcf':'Remove',
    'CCB020081_0200.snp.filtered.snpeff.vcf':'Remove',
    'CCB020104_200.snp.filtered.snpeff.vcf':'dogBLSA',
    'CCB020156_0203.snp.filtered.snpeff.vcf':'Remove',
    'CCB020249_200.snp.filtered.snpeff.vcf':'dogBLSA',
    'CCB030003_200.snp.filtered.snpeff.vcf':'dogTLSA',
    'CCB030016_201.snp.filtered.snpeff.vcf':'dogTLSA',
    'CCB030021_0200.snp.filtered.snpeff.vcf':'dogBLSA',
    'CCB030026_201.snp.filtered.snpeff.vcf':'dogTLSA',
    'CCB030059_201.snp.filtered.snpeff.vcf':'dogTLSA',
    'CCB030061_0201.snp.filtered.snpeff.vcf':'Remove',
    'CCB030078_0201.snp.filtered.snpeff.vcf':'Remove',
    'CCB030108_0200.snp.filtered.snpeff.vcf':'Remove',
    'CCB030117_201.snp.filtered.snpeff.vcf':'dogBLSA',
    'CCB030130_200.snp.filtered.snpeff.vcf':'dogBLSA',
    'CCB030145_201.snp.filtered.snpeff.vcf':'dogBLSA',
    'CCB030156_200.snp.filtered.snpeff.vcf':'dogTLSA',
    'CCB030166_0200.snp.filtered.snpeff.vcf':'dogTLSA',
    'CCB030193_0200.snp.filtered.snpeff.vcf':'dogBLSA',
    'CCB030212_0201.snp.filtered.snpeff.vcf':'dogTLSA',
    'CCB030224_0201.snp.filtered.snpeff.vcf':'Remove',
    'CCB030251_200.snp.filtered.snpeff.vcf':'dogBLSA',
    'CCB030258_200.snp.filtered.snpeff.vcf':'dogTLSA',
    'CCB030259_200.snp.filtered.snpeff.vcf':'dogTLSA',
    'CCB030265_201.snp.filtered.snpeff.vcf':'dogTLSA',
    'CCB030268_201.snp.filtered.snpeff.vcf':'dogBLSA',
    'CCB030271_201.snp.filtered.snpeff.vcf':'dogTLSA',
    'CCB030292_0201.snp.filtered.snpeff.vcf':'dogTLSA',
    'CCB030293_0201.snp.filtered.snpeff.vcf':'Remove',
    'CCB030300_201.snp.filtered.snpeff.vcf':'dogTLSA',
    'CCB030317_200.snp.filtered.snpeff.vcf':'dogTLSA',
    'CCB030346_201.snp.filtered.snpeff.vcf':'dogBLSA',
    'CCB030357_200.snp.filtered.snpeff.vcf':'Remove',
    'CCB030362_200.snp.filtered.snpeff.vcf':'dogBLSA',
    'CCB040007_0200.snp.filtered.snpeff.vcf':'Remove',
    'CCB040033_0200.snp.filtered.snpeff.vcf':'Remove',
    'CCB040038_201.snp.filtered.snpeff.vcf':'dogBLSA',
    'CCB040105_201.snp.filtered.snpeff.vcf':'Remove',
    'CCB040143_200.snp.filtered.snpeff.vcf':'dogBLSA',
    'CCB040180_0201.snp.filtered.snpeff.vcf':'dogBLSA',
    'CCB040184_0200.snp.filtered.snpeff.vcf':'Remove',
    'CCB040207_0200.snp.filtered.snpeff.vcf':'Remove',
    'CCB040222_200.snp.filtered.snpeff.vcf':'dogBLSA',
    'CCB040232_0201.snp.filtered.snpeff.vcf':'Remove',
    'CCB040265_0200.snp.filtered.snpeff.vcf':'dogTLSA',
    'CCB040287_0200.snp.filtered.snpeff.vcf':'Remove',
    'CCB040292_201.snp.filtered.snpeff.vcf':'dogTLSA',
    'CCB040313_200.snp.filtered.snpeff.vcf':'Remove',
    'CCB040425_201.snp.filtered.snpeff.vcf':'dogTLSA',
    'CCB040464_0201.snp.filtered.snpeff.vcf':'Remove',
    'CCB040497_201.snp.filtered.snpeff.vcf':'dogBLSA',
    'CCB050015_0200.snp.filtered.snpeff.vcf':'Remove',
    'CCB050025_0200.snp.filtered.snpeff.vcf':'Remove',
    'CCB050063_0200.snp.filtered.snpeff.vcf':'Remove',
    'CCB050102_0201.snp.filtered.snpeff.vcf':'Remove',
    'CCB050193_200.snp.filtered.snpeff.vcf':'dogBLSA',
    'CCB050194_0200.snp.filtered.snpeff.vcf':'Remove',
    'CCB050195_0200.snp.filtered.snpeff.vcf':'Remove',
    'CCB050204_0201.snp.filtered.snpeff.vcf':'Remove',
    'CCB050206_200.snp.filtered.snpeff.vcf':'dogBLSA',
    'CCB050315_0200.snp.filtered.snpeff.vcf':'Remove',
    'CCB050320_0200.snp.filtered.snpeff.vcf':'Remove',
    'CCB060008_200.snp.filtered.snpeff.vcf':'dogTLSA',
    'CCB060021_0200.snp.filtered.snpeff.vcf':'Remove',
    'CCB060027_0201.snp.filtered.snpeff.vcf':'dogTLSA',
    'CCB060032_0201.snp.filtered.snpeff.vcf':'dogTLSA',
    'CCB060033_0201.snp.filtered.snpeff.vcf':'Remove',
    'CCB060034_0201.snp.filtered.snpeff.vcf':'dogTLSA',
    'CCB060047_201.snp.filtered.snpeff.vcf':'dogBLSA',
    'CCB060052_200.snp.filtered.snpeff.vcf':'dogBLSA',
    'CCB060066_0201.snp.filtered.snpeff.vcf':'Remove',
    'CCB060086_0201.snp.filtered.snpeff.vcf':'Remove',
    'CCB060087_0201.snp.filtered.snpeff.vcf':'Remove',
    'CCB060088_0200.snp.filtered.snpeff.vcf':'Remove',
    'CCB060100_0200.snp.filtered.snpeff.vcf':'Remove',
    'CCB060135_0201.snp.filtered.snpeff.vcf':'dogBLSA',
    'CCB060151.snp.filtered.snpeff.vcf':'dogTLSA',
    'CCB060153_200.snp.filtered.snpeff.vcf':'dogBLSA',
    'CCB070011_0200.snp.filtered.snpeff.vcf':'Remove',
    'CCB070012_200.snp.filtered.snpeff.vcf':'dogBLSA',
    'CCB070060_0201.snp.filtered.snpeff.vcf':'Remove',
    'CCB070080_0201.snp.filtered.snpeff.vcf':'dogTLSA',
    'CCB070081_201.snp.filtered.snpeff.vcf':'dogBLSA',
    'CCB070086_0201.snp.filtered.snpeff.vcf':'dogBLSA',
    'CCB070115_200.snp.filtered.snpeff.vcf':'dogBLSA',
    'CCB070124_0200.snp.filtered.snpeff.vcf':'dogBLSA',
    'CCB070131_0200.snp.filtered.snpeff.vcf':'dogBLSA',
    'CCB070139_0201.snp.filtered.snpeff.vcf':'dogTLSA',
    'CCB070162_0200.snp.filtered.snpeff.vcf':'dogTLSA',
    'CCB070177_0201.snp.filtered.snpeff.vcf':'Remove',
    'CCB070230_201.snp.filtered.snpeff.vcf':'dogBLSA',
    'CCB070243_0200.snp.filtered.snpeff.vcf':'dogBLSA',
    'CCB070252_201.snp.filtered.snpeff.vcf':'dogBLSA',
    'CCB070253_0200.snp.filtered.snpeff.vcf':'Remove',
    'CCB070280_200.snp.filtered.snpeff.vcf':'dogBLSA',
    'CCB070289_0200.snp.filtered.snpeff.vcf':'dogBLSA',
    'CCB070292_0201.snp.filtered.snpeff.vcf':'dogTLSA',
    'CCB080009_0200.snp.filtered.snpeff.vcf':'Remove',
    'CCB080016_0200.snp.filtered.snpeff.vcf':'Remove',
    'CCB080023_0201.snp.filtered.snpeff.vcf':'Remove',
    'CasLea_Normal.snp.filtered.snpeff.vcf':'dogTLSA',
    'MB140_Normal.snp.filtered.snpeff.vcf':'dogBLSA',
    'MB2A11_Normal.snp.filtered.snpeff.vcf':'dogBLSA',
    'MB2C13_Normal.snp.filtered.snpeff.vcf':'dogBLSA',
    'MurVoe_Normal.snp.filtered.snpeff.vcf':'dogTLSA',
    'NEWMET_Normal.snp.filtered.snpeff.vcf':'dogTLSA',
    'NygRoy_Normal.snp.filtered.snpeff.vcf':'dogTLSA',
    'RoxRho_Normal.snp.filtered.snpeff.vcf':'dogTLSA',
    'SHAKIT_Normal.snp.filtered.snpeff.vcf':'dogTLSA',
    'SM-2977C.snp.filtered.snpeff.vcf':'Remove',
    'SM-2977E.snp.filtered.snpeff.vcf':'dogBLSA',
    'SM-2977G.snp.filtered.snpeff.vcf':'dogBLSA',
    'SM-2977I.snp.filtered.snpeff.vcf':'dogBLSA',
    'SM-2977K.snp.filtered.snpeff.vcf':'dogBLSA',
    'SM-2977S.snp.filtered.snpeff.vcf':'dogBLSA',
    'SM-2977W.snp.filtered.snpeff.vcf':'dogBLSA',
    'SM-29781.snp.filtered.snpeff.vcf':'dogBLSA',
    'SM-29783.snp.filtered.snpeff.vcf':'dogBLSA',
    'SM-29785.snp.filtered.snpeff.vcf':'dogBLSA',
    'SM-2978R.snp.filtered.snpeff.vcf':'dogBLSA',
    'SM-2978V.snp.filtered.snpeff.vcf':'Remove',
    'SM-2978X.snp.filtered.snpeff.vcf':'dogBLSA',
    'SM-3JNGX.snp.filtered.snpeff.vcf':'dogBLSA',
    'SM-3LE5G.snp.filtered.snpeff.vcf':'dogBLSA',
    'SM-3LE5H.snp.filtered.snpeff.vcf':'dogBLSA',
    'SM-3LE5I.snp.filtered.snpeff.vcf':'dogBLSA',
    'SM-3LE5L.snp.filtered.snpeff.vcf':'dogBLSA',
    'SM-3LE5N.snp.filtered.snpeff.vcf':'Remove',
    'SopCas_Normal.snp.filtered.snpeff.vcf':'dogTLSA',
    'SopDix_Normal.snp.filtered.snpeff.vcf':'dogTLSA',
    'SteHel_Normal.snp.filtered.snpeff.vcf':'dogTLSA',
    }                  
BOCA_UK_lookup = [] 
ESCA_lookup = []
SARC_lookup = []

with open('/seq/vgb/cancer_r01/icgc/release_28/BOCA-UK/BOCA-UK_OSA.txt', 'r') as inFile:
    for line in inFile:
            fields = line.strip().split('\t')
            sample_id = fields[0].strip() #should perhaps be 0?
            if sample_id != 'icgc_sample_id':
                #print('lookup file name: ' + sample_id)
                BOCA_UK_lookup.append(sample_id)
                
with open('/seq/vgb/cancer_r01/tcga_data/annotations/ESCA/ESCA_adenosarcoma.txt', 'r') as inFileOne, open('/seq/vgb/cancer_r01/tcga_data/annotations/ESCA/ESCA_scc.txt', 'r') as inFileTwo:
    for line in inFileOne:
        sample_id = line.strip()
        ESCA_lookup.append(sample_id)
    for line in inFileTwo:
        sample_id = line.strip()
        ESCA_lookup.append(sample_id)
        
with open('/seq/vgb/cancer_r01/tcga_data/annotations/SARC/SARC_leiomyosarcomas.txt', 'r') as inFile:
    for line in inFile:
        sample_id = line.strip()
        SARC_lookup.append(sample_id)
                    
variantIDLookup = {}
with open('/seq/vgb/swofford/R01/convertingToSKLInput/testVCFs/singleBase/hg19CancerHg38.bed', 'r') as inFile:
    for line in inFile:
        fields = line.strip().split('\t')
        hg38 = '@'.join([fields[0], fields[2], 'hg38'])
        hg19 = fields[3]
        variantIDLookup[hg19] = hg38

with open('/seq/vgb/swofford/R01/convertingToSKLInput/testVCFs/singleBase/canFam3CancerHg38.bed', 'r') as inFile:
    for line in inFile:
        fields = line.strip().split('\t')
        hg38 = '@'.join([fields[0], fields[2], 'hg38'])
        canFam3 = fields[3]
        variantIDLookup[canFam3] = hg38

with open('/seq/vgb/swofford/R01/convertingToSKLInput/testVCFs/singleBase/hg38Cancer.bed', 'r') as inFile:
    for line in inFile:
        fields = line.strip().split('\t')
        hg38 = '@'.join([fields[0], fields[2], 'hg38'])
        variantIDLookup[hg38] = hg38               
# Go through individual files and collect variant info for each sample
previouslySeen = []

'''
variantToCodon = {}
with open('/seq/vgb/swofford/R01/convertingToSKLInput/testVCFs/singleBase/allCancerHg38LookupLabeled.bed', 'r') as inFile:
    for line in inFile:
        fields = line.strip().split('\t')
        try:
            variantToCodon[fields[3]].add('@'.join(fields[0:3]) + '@' + fields[7])
        except KeyError:
            variantToCodon[fields[3]] = set()
            variantToCodon[fields[3]].add('@'.join(fields[0:3]) + '@' + fields[7])
'''            
            
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
                    
                    
                    if inCDS('dogFast', chrom, pos):
                        newVariant = variant(Gene_Name = Gene_Name, Transcript_ID = Transcript_ID, 
                                             Genotype = Genotype, Effect = Effect, chrom = chrom, pos = pos, ref = ref, assembly = 'canFam3') 
                        try:
                            sampleVariants[name].append(newVariant)
                        except KeyError:
                            sampleVariants[name] = [newVariant]

                for name in sampleVariants:
                    if len(sampleVariants[name]) < maxVar and len(sampleVariants[name]) >= minVar and name not in noFlyList and label != 'Remove':
                        samples.append(sample(variants = sampleVariants[name], name = name, species = 'dog', label = label, source = fileName.strip()))
                    elif name in noFlyList:
                        with open('noFly.out', 'a') as outFile:
                            outFile.write(name + '\n')
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
                                         Genotype = Genotype, Effect = Effect, chrom = chrom, pos = pos, ref = ref, assembly = 'hg38')
                     
                    try:
                        sampleVariants[name].append(newVariant)
                    except KeyError:
                        sampleVariants[name] = [newVariant]

                for name in sampleVariants:
                    if len(sampleVariants[name]) < maxVar:
                        samples.append(sample(sampleVariants[name], name, 'Human', label))
        '''            
                
        #Ensembl_so_term
        if 'data_mutations' in fileName:
            
            if 'dlbc_broad' in fileName: #special parsing for txt file, while regular file missing
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
                            names = line.strip().split('\t')
                            break
                        
                    phenoReader = csv.DictReader(pheno_in, delimiter = '\t', fieldnames = names)
                    for line in phenoReader:
                        sampleID = line['SAMPLE_ID']
                        study = line['ONCOTREE_CODE']
                        sampleLabel[sampleID] = studyClass[study]
                        
                    for line in maf_in:
                        if line[:1] != '#':
                            names = line.strip().split('\t')
                            break
                    
                    mafReader = csv.DictReader(maf_in, delimiter = '\t', fieldnames=names)
    
                    
                    csv.field_size_limit(int(sys.maxsize/10))
                    for line in mafReader:
                        name = line['Tumor_Sample_Barcode']
                        Exon_Rank = '' #missing
                        Gene_Name = line['Hugo_Symbol'] # good
                        Amino_Acid_Change = '' #missing
                        Transcript_ID = line['Transcript_ID'] #missing
                        Genotype = line['Tumor_Seq_Allele2'] #good
                        Effect = line['Consequence'] #missing - replaced
                        chrom = 'chr' + str(line['Chromosome']) #good
                        try:
                            pos = int(line['Start_Position']) #slightly different - replaced
                            #print('pos is working: ' + str(pos))
                        except ValueError:
                            pos = 'NA'
                        ref = line['Reference_Allele'] #good
                        #print('ref is working: ' + ref)
                        
                        
                        if inCDS('hg19Fast',chrom,pos):
                            #print('ever in cds?')
                            newVariant = variant(Amino_Acid_Change = Amino_Acid_Change, 
                                                 Gene_Name = Gene_Name, Exon_Rank = Exon_Rank, Transcript_ID = Transcript_ID, 
                                                 Genotype = Genotype, Effect = Effect, chrom = chrom, pos = pos, ref = ref, assembly = 'hg19')
                            try:
                                sampleVariants[name].append(newVariant)
                            except:
                                sampleVariants[name] = [newVariant]
    
                    for name in sampleVariants:                       
                        try:
                            if len(sampleVariants[name]) < maxVar and len(sampleVariants[name]) >= minVar and sampleLabel[name] != 'Remove' and name not in previouslySeen:
                            #print(name + ' has ' + str(len(sampleVariants[name])) + ' variants')
                                samples.append(sample(variants = sampleVariants[name], name = name, species = 'Human', label = sampleLabel[name], source = fileName.strip()))
                                previouslySeen.append(name)
                        except KeyError:
                            print('Sample name not found in phenotype file: ' + name)
            
            elif 'dlbcl_dfci_2018' in fileName: 
                print('here!')                           
                basePath = ''
                fullPath = fileName.split('/')
                for i in range(len(fullPath) - 1):
                    basePath = basePath + fullPath[i] + '/'
                phenotypeFile = basePath + 'data_clinical_sample.txt'
                print(phenotypeFile)
                
    
                with open(fileName, 'rt') as maf_in, open(phenotypeFile, 'rt') as pheno_in:
                    sampleLabel = {}
                    sampleVariants = {}
                    
                    # skip header, if there is one, and get column names
                    for line in pheno_in:
                        if line[:1] != '#':
                            line = line.strip()
                            names = line.strip().split('\t')
                            break
                        
                    phenoReader = csv.DictReader(pheno_in, delimiter = '\t', fieldnames = names)
                    for line in phenoReader:
                        #pediatric pan-cancer sample names end in '-R' when sample is from a relapse, rather than primary
                        #We're going to ignore the relapse samples
                        
                        sampleID = line['SAMPLE_ID'].strip()
                        study = line['CANCER_TYPE'].strip()
                        sampleLabel[sampleID] = studyClass[study]
                        
                    for line in maf_in:
                        if line[:1] != '#':
                            names = line.strip().split('\t')
                            break
                    
                    mafReader = csv.DictReader(maf_in, delimiter = '\t', fieldnames=names)
    
                    
                    
                    for line in mafReader:
                        name = line['Tumor_Sample_Barcode'].strip()
                        Exon_Rank = line['Protein_position'] #missing
                        Gene_Name = line['Hugo_Symbol'] # good
                        Amino_Acid_Change = line['HGVSp'] #missing
                        Transcript_ID = line['Transcript_ID'] #missing
                        Genotype = line['Tumor_Seq_Allele2'] #good
                        Effect = line['Consequence'] #missing
                        chrom = 'chr' + str(line['Chromosome']) #good
                        try:
                            pos = int(line['Start_Position']) #good
                        except ValueError:
                            pos = 'NA'
                        ref = line['Reference_Allele'] #good
                        
                        
                        if inCDS('hg19Fast', chrom, pos):
                            newVariant = variant(Amino_Acid_Change = Amino_Acid_Change, 
                                                 Gene_Name = Gene_Name, Exon_Rank = Exon_Rank, Transcript_ID = Transcript_ID, 
                                                 Genotype = Genotype, Effect = Effect, chrom = chrom, pos = pos, ref = ref, assembly = 'hg19')
                        
                            try:
                                sampleVariants[name].append(newVariant)
                            except:
                                sampleVariants[name] = [newVariant]
    
                    for name in sampleVariants:
                        try:
                            if sampleLabel[name] != 'Remove' and len(sampleVariants[name]) < maxVar and len(sampleVariants[name]) >= minVar:
                            #print(name)
                                samples.append(sample(variants = sampleVariants[name], name = name, species = 'Human', label = sampleLabel[name], source = fileName.strip()))
                                previouslySeen.append(name)
                        except KeyError:
                            print('Sample name not found in phenotype file: ' + name)
            
            
            
            
            
            
            
            
            
            elif 'angs_painter_2020' in fileName: #special parsing for extended file, while regular file missing
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
                            names = line.strip().split('\t')
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
                            names = line.strip().split('\t')
                            break
                    
                    mafReader = csv.DictReader(maf_in, delimiter = '\t', fieldnames=names)
    
                    
                    csv.field_size_limit(int(sys.maxsize/10))
                    for line in mafReader:
                        name = line['Tumor_Sample_Barcode']
                        Exon_Rank = '' #missing
                        Gene_Name = line['Hugo_Symbol'] # good
                        Amino_Acid_Change = '' #missing
                        Transcript_ID = '' #missing
                        Genotype = line['Tumor_Seq_Allele2'] #good
                        Effect = line['Ensembl_so_term'] #missing - replaced
                        chrom = 'chr' + str(line['Chromosome']) #good
                        try:
                            pos = int(line['Start_position']) #slightly different - replaced
                            #print('pos is working: ' + str(pos))
                        except ValueError:
                            pos = 'NA'
                        ref = line['Reference_Allele'] #good
                        #print('ref is working: ' + ref)
                        
                        
                        if inCDS('hg19Fast',chrom,pos):
                            #print('ever in cds?')
                            newVariant = variant(Amino_Acid_Change = Amino_Acid_Change, 
                                                 Gene_Name = Gene_Name, Exon_Rank = Exon_Rank, Transcript_ID = Transcript_ID, 
                                                 Genotype = Genotype, Effect = Effect, chrom = chrom, pos = pos, ref = ref, assembly = 'hg19')
                            try:
                                sampleVariants[name].append(newVariant)
                            except:
                                sampleVariants[name] = [newVariant]
    
                    for name in sampleVariants:                       
                        try:
                            if len(sampleVariants[name]) < maxVar and len(sampleVariants[name]) >= minVar and sampleLabel[name] != 'Remove':
                            #print(name + ' has ' + str(len(sampleVariants[name])) + ' variants')
                                samples.append(sample(variants = sampleVariants[name], name = name, species = 'Human', label = sampleLabel[name], source = fileName.strip()))
                        except KeyError:
                            print('Sample name not found in phenotype file: ' + name)
            
            elif 'pediatric_columbia_2019' in fileName:
                basePath = ''
                fullPath = fileName.split('/')
                for i in range(len(fullPath) - 1):
                    basePath = basePath + fullPath[i] + '/'
                phenotypeFile = basePath + 'data_clinical_sample.txt'
                print(phenotypeFile)
                
    
                with open(fileName, 'rt') as maf_in, open(phenotypeFile, 'rt') as pheno_in:
                    sampleLabel = {}
                    sampleVariants = {}
                    
                    # skip header, if there is one, and get column names
                    for line in pheno_in:
                        if line[:1] != '#':
                            line = line.strip()
                            names = line.strip().split('\t')
                            break
                        
                    phenoReader = csv.DictReader(pheno_in, delimiter = '\t', fieldnames = names)
                    for line in phenoReader:
                        #pediatric pan-cancer sample names end in '-R' when sample is from a relapse, rather than primary
                        #We're going to ignore the relapse samples
                        
                        if line['SAMPLE_ID'].strip() in pediatric_columbia_2019_keep:
                            sampleID = line['SAMPLE_ID'].strip()
                            study = line['CANCER_TYPE'].strip()
                            sampleLabel[sampleID] = studyClass[study]
                        else: print(line['SAMPLE_ID'].strip() + ' not in keep')
                            
                        
                    for line in maf_in:
                        if line[:1] != '#':
                            names = line.strip().split('\t')
                            break
                    
                    mafReader = csv.DictReader(maf_in, delimiter = '\t', fieldnames=names)
    
                    
                    
                    for line in mafReader:
                        name = line['Tumor_Sample_Barcode']
                        Exon_Rank = line['Protein_position'] #missing
                        Gene_Name = line['Hugo_Symbol'] # good
                        Amino_Acid_Change = line['HGVSp'] #missing
                        Transcript_ID = line['Transcript_ID'] #missing
                        Genotype = line['Tumor_Seq_Allele2'] #good
                        Effect = line['Consequence'] #missing
                        chrom = 'chr' + str(line['Chromosome']) #good
                        try:
                            pos = int(line['Start_Position']) #good
                        except ValueError:
                            pos = 'NA'
                        ref = line['Reference_Allele'] #good
                        
                        
                        if inCDS('hg19Fast', chrom, pos):
                            newVariant = variant(Amino_Acid_Change = Amino_Acid_Change, 
                                                 Gene_Name = Gene_Name, Exon_Rank = Exon_Rank, Transcript_ID = Transcript_ID, 
                                                 Genotype = Genotype, Effect = Effect, chrom = chrom, pos = pos, ref = ref, assembly = 'hg19')
                            try:
                                sampleVariants[name].append(newVariant)
                            except:
                                sampleVariants[name] = [newVariant]
    
                    for name in sampleVariants:                      
                        try:
                            if len(sampleVariants[name]) < maxVar and len(sampleVariants[name]) >= minVar and sampleLabel[name] != 'Remove':
                            #print(name)
                                samples.append(sample(variants = sampleVariants[name], name = name, species = 'Human', label = sampleLabel[name], source = fileName.strip()))
                        except KeyError:
                            print('Sample name not found in phenotype file: ' + name)
            
            elif 'brca_mbcproject_wagle_2017' in fileName:
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
                            names = line.strip().split('\t')
                            break
                        
                    phenoReader = csv.DictReader(pheno_in, delimiter = '\t', fieldnames = names)
                    for line in phenoReader:
                        #pediatric pan-cancer sample names end in '-R' when sample is from a relapse, rather than primary
                        #We're going to ignore the relapse samples
                        
                        if line['SAMPLE_ID'] in brca_mbcproject_wagle_2017_keep:
                            sampleID = line['SAMPLE_ID']
                            study = line['CANCER_TYPE']
                            try:
                                sampleLabel[sampleID] = studyClass[study]
                            except KeyError:
                                sampleLabel = 'Remove'
                        
                    for line in maf_in:
                        if line[:1] != '#':
                            names = line.strip().split('\t')
                            break
                    
                    mafReader = csv.DictReader(maf_in, delimiter = '\t', fieldnames=names)
    
                    
                    
                    for line in mafReader:
                        name = line['Tumor_Sample_Barcode']
                        Exon_Rank = line['Protein_position'] #missing
                        Gene_Name = line['Hugo_Symbol'] # good
                        Amino_Acid_Change = line['HGVSp'] #missing
                        Transcript_ID = line['Transcript_ID'] #missing
                        Genotype = line['Tumor_Seq_Allele2'] #good
                        Effect = line['Consequence'] #missing
                        chrom = 'chr' + str(line['Chromosome']) #good
                        try:
                            pos = int(line['Start_Position']) #good
                        except ValueError:
                            pos = 'NA'
                        ref = line['Reference_Allele'] #good
                        
                        
                        
                        if inCDS('hg19Fast',chrom,pos):
                            newVariant = variant(Amino_Acid_Change = Amino_Acid_Change, 
                                                 Gene_Name = Gene_Name, Exon_Rank = Exon_Rank, Transcript_ID = Transcript_ID, 
                                                 Genotype = Genotype, Effect = Effect, chrom = chrom, pos = pos, ref = ref, assembly = 'hg19')
                            try:
                                sampleVariants[name].append(newVariant)
                            except:
                                sampleVariants[name] = [newVariant]
    
                    for name in sampleVariants:
                        try:
                            #print(name)
                            if sampleLabel[name] == 'Remove' or len(sampleVariants[name]) > maxVar or len(sampleVariants[name]) < minVar:
                                continue
                            else:
                                samples.append(sample(variants = sampleVariants[name], name = name, species = 'Human', label = sampleLabel[name], source = fileName.strip()))
                        except KeyError:
                            print('Sample name not found in phenotype file: ' + name)
                            
            elif 'prad_mpcproject_2018' in fileName:
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
                            names = line.strip().split('\t')
                            break
                        
                    phenoReader = csv.DictReader(pheno_in, delimiter = '\t', fieldnames = names)
                    for line in phenoReader:
                        #pediatric pan-cancer sample names end in '-R' when sample is from a relapse, rather than primary
                        #We're going to ignore the relapse samples
                        
                        if line['SAMPLE_ID'] in prad_mpcproject_2018_keep:
                            sampleID = line['SAMPLE_ID']
                            study = line['CANCER_TYPE']
                            sampleLabel[sampleID] = studyClass[study]
                            
                        
                    for line in maf_in:
                        if line[:1] != '#':
                            names = line.strip().split('\t')
                            break
                    
                    mafReader = csv.DictReader(maf_in, delimiter = '\t', fieldnames=names)
    
                    
                    
                    for line in mafReader:
                        name = line['Tumor_Sample_Barcode']
                        Exon_Rank = line['Protein_position'] #missing
                        Gene_Name = line['Hugo_Symbol'] # good
                        Amino_Acid_Change = line['HGVSp'] #missing
                        Transcript_ID = line['Transcript_ID'] #missing
                        Genotype = line['Tumor_Seq_Allele2'] #good
                        Effect = line['Consequence'] #missing
                        chrom = 'chr' + str(line['Chromosome']) #good
                        try:
                            pos = int(line['Start_Position']) #good
                        except ValueError:
                            pos = 'NA'
                        ref = line['Reference_Allele'] #good
                        
                        
                        if inCDS('hg19Fast', chrom, pos):
                        
                            newVariant = variant(Amino_Acid_Change = Amino_Acid_Change, 
                                                 Gene_Name = Gene_Name, Exon_Rank = Exon_Rank, Transcript_ID = Transcript_ID, 
                                                 Genotype = Genotype, Effect = Effect, chrom = chrom, pos = pos, ref = ref, assembly = 'hg19')
                    
                            try:
                                sampleVariants[name].append(newVariant)
                            except:
                                sampleVariants[name] = [newVariant]
    
                    for name in sampleVariants:                      
                        try:
                            if sampleLabel[name] != 'Remove' and len(sampleVariants[name]) < maxVar and len(sampleVariants[name]) >= minVar:
                            #print(name)
                                samples.append(sample(variants = sampleVariants[name], name = name, species = 'Human', label = sampleLabel[name], source = fileName.strip()))
                        except KeyError:
                            print('Sample name not found in phenotype file: ' + name)
            
            elif 'metastatic_solid_tumors_mich_2017' in fileName:
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
                            names = line.strip().split('\t')
                            break
                        
                    phenoReader = csv.DictReader(pheno_in, delimiter = '\t', fieldnames = names)
                    for line in phenoReader:
                        #pediatric pan-cancer sample names end in '-R' when sample is from a relapse, rather than primary
                        #We're going to ignore the relapse samples
                        
                        if line['SAMPLE_ID'] in metastatic_solid_tumors_mich_2017_keep:
                            sampleID = line['SAMPLE_ID']
                            study = line['CANCER_TYPE']
                            sampleLabel[sampleID] = studyClass[study]
                            
                        
                    for line in maf_in:
                        if line[:1] != '#':
                            names = line.strip().split('\t')
                            break
                    
                    mafReader = csv.DictReader(maf_in, delimiter = '\t', fieldnames=names)
    
                    
                    
                    for line in mafReader:
                        name = line['Tumor_Sample_Barcode']
                        Exon_Rank = line['Protein_position'] #missing
                        Gene_Name = line['Hugo_Symbol'] # good
                        Amino_Acid_Change = line['HGVSp'] #missing
                        Transcript_ID = line['Transcript_ID'] #missing
                        Genotype = line['Tumor_Seq_Allele2'] #good
                        Effect = line['Consequence'] #missing
                        chrom = 'chr' + str(line['Chromosome']) #good
                        try:
                            pos = int(line['Start_Position']) #good
                        except ValueError:
                            pos = 'NA'
                        ref = line['Reference_Allele'] #good
                        
                        
                        if inCDS('hg19Fast', chrom, pos):
                        
                            newVariant = variant(Amino_Acid_Change = Amino_Acid_Change, 
                                                 Gene_Name = Gene_Name, Exon_Rank = Exon_Rank, Transcript_ID = Transcript_ID, 
                                                 Genotype = Genotype, Effect = Effect, chrom = chrom, pos = pos, ref = ref, assembly = 'hg19')
                        
                            try:
                                sampleVariants[name].append(newVariant)
                            except:
                                sampleVariants[name] = [newVariant]
    
                    for name in sampleVariants:
                        try:
                            if sampleLabel[name] != 'Remove' and len(sampleVariants[name]) < maxVar and len(sampleVariants[name]) >= minVar:
                            #print(name)
                                samples.append(sample(variants = sampleVariants[name], name = name, species = 'Human', label = sampleLabel[name], source = fileName.strip()))
                        except KeyError:
                            print('Sample name not found in phenotype file: ' + name)
                
            
            else:                            
                basePath = ''
                fullPath = fileName.split('/')
                for i in range(len(fullPath) - 1):
                    basePath = basePath + fullPath[i] + '/'
                phenotypeFile = basePath + 'data_clinical_sample.txt'
                print(phenotypeFile)
                
    
                with open(fileName, 'rt') as maf_in, open(phenotypeFile, 'rt') as pheno_in:
                    sampleLabel = {}
                    sampleVariants = {}
                    
                    # skip header, if there is one, and get column names
                    for line in pheno_in:
                        if line[:1] != '#':
                            line = line.strip()
                            names = line.strip().split('\t')
                            break
                        
                    phenoReader = csv.DictReader(pheno_in, delimiter = '\t', fieldnames = names)
                    for line in phenoReader:
                        #pediatric pan-cancer sample names end in '-R' when sample is from a relapse, rather than primary
                        #We're going to ignore the relapse samples
                        
                        if studyClass[line['CANCER_TYPE']] != 'AS':
                            if '-R' not in line['SAMPLE_ID']:
                                sampleID = line['SAMPLE_ID'].strip()
                                study = line['CANCER_TYPE'].strip()
                                sampleLabel[sampleID] = studyClass[study]
                        else:
                            if line['SAMPLE_ID'] in ASKeep:
                                sampleID = line['SAMPLE_ID']
                                study = line['CANCER_TYPE']
                                sampleLabel[sampleID] = studyClass[study]
                        
                    for line in maf_in:
                        if line[:1] != '#':
                            names = line.strip().split('\t')
                            break
                    
                    mafReader = csv.DictReader(maf_in, delimiter = '\t', fieldnames=names)
    
                    
                    
                    for line in mafReader:
                        name = line['Tumor_Sample_Barcode'].strip()
                        Exon_Rank = line['Protein_position'] #missing
                        Gene_Name = line['Hugo_Symbol'] # good
                        Amino_Acid_Change = line['HGVSp'] #missing
                        Transcript_ID = line['Transcript_ID'] #missing
                        Genotype = line['Tumor_Seq_Allele2'] #good
                        Effect = line['Consequence'] #missing
                        chrom = 'chr' + str(line['Chromosome']) #good
                        try:
                            pos = int(line['Start_Position']) #good
                        except ValueError:
                            pos = 'NA'
                        ref = line['Reference_Allele'] #good
                        
                        
                        if inCDS('hg19Fast', chrom, pos):
                            newVariant = variant(Amino_Acid_Change = Amino_Acid_Change, 
                                                 Gene_Name = Gene_Name, Exon_Rank = Exon_Rank, Transcript_ID = Transcript_ID, 
                                                 Genotype = Genotype, Effect = Effect, chrom = chrom, pos = pos, ref = ref, assembly = 'hg19')
                        
                            try:
                                sampleVariants[name].append(newVariant)
                            except:
                                sampleVariants[name] = [newVariant]
    
                    for name in sampleVariants:
                        try:
                            if sampleLabel[name] != 'Remove' and len(sampleVariants[name]) < maxVar and len(sampleVariants[name]) >= minVar:
                            #print(name)
                                samples.append(sample(variants = sampleVariants[name], name = name, species = 'Human', label = sampleLabel[name], source = fileName.strip()))
                        except KeyError:
                            print('Sample name not found in phenotype file: ' + name)
                        
                        
        elif 'TCGA' in fileName:
            #assumes that TCGA filenames all follow the same naming convention
            #species = 'human'
            sampleVariants = {}
            try:
                label = studyClass[fileName.split('.')[1]]
            except KeyError:
                print('missing study class: ' + studyClass[fileName.split('.')[1]])
                label = 'Remove'
            if label == 'Remove':
                continue
            
            print(label)
            
            with gzip.open(fileName, 'rt') as maf_in:
                for line in maf_in:
                    if line[:1] != '#':
                        names = line.strip().split('\t')
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
                    
                    if inCDS('humFast', chrom, pos):
                        newVariant = variant(Amino_Acid_Change = Amino_Acid_Change, 
                                             Gene_Name = Gene_Name, Exon_Rank = Exon_Rank, Transcript_ID = Transcript_ID, 
                                             Genotype = Genotype, Effect = Effect, chrom = chrom, pos = pos, ref = ref, assembly = 'hg38')
                        #print(newVariant.chrom + '\t' + str(newVariant.pos) + '\t' + newVariant.ref + '\t' + str(newVariant.breaksGene))
                        try:
                            sampleVariants[name].append(newVariant)
                        except KeyError:
                            sampleVariants[name] = [newVariant]

                for name in sampleVariants:
                    if len(sampleVariants[name]) < maxVar and len(sampleVariants[name]) >= minVar and name not in noFlyList:
                        if fileName.split('.')[1] == 'SARC':
                            if name in SARC_lookup:
                                samples.append(sample(variants = sampleVariants[name], name = name, species = 'Human', label = label, source = fileName.strip()))
                            else:
                                continue
                        elif fileName.split('.')[1] == 'ESCA':
                            if name in ESCA_lookup:
                                samples.append(sample(variants = sampleVariants[name], name = name, species = 'Human', label = label, source = fileName.strip()))
                            else:
                                continue
                        else:
                            samples.append(sample(variants = sampleVariants[name], name = name, species = 'Human', label = label, source = fileName.strip()))
                        #print(name + '\t' + str(sample(sampleVariants[name], name, 'Human', label).brokenGenes))
                    elif name in noFlyList:
                        with open('noFly.out', 'a') as outFile:
                            outFile.write(name + '\n')
                        
        elif 'nogistic' in fileName:
            #assumes that TCGA filenames all follow the same naming convention
            #species = 'human'
            sampleVariants = {}
            
            label = 'dogGlioma'
            
            print(label)
            
            with open(fileName, 'rt') as maf_in:
                for line in maf_in:
                    if line[:1] != '#':
                        names = line.strip().split('\t')
                        break
                    
                mafReader = csv.DictReader(maf_in, delimiter = '\t', fieldnames=names)
                variantList = []
                oldName = ''
                
                
                for line in mafReader:
                    name = line['Tumor_Sample_Barcode']
                    if name in canineGliomaKeep:
                        if oldName == '':
                            oldName = name
                                     
                        Amino_Acid_Change = line['HGVSp_Short']
                        Exon_Rank = ''
                        Gene_Name = line['Hugo_Symbol']
                        #Amino_Acid_Change = line['HGVSp']
                        Transcript_ID = line['ENSG_ID'] # this will be gene id rather than transcript
                        Genotype = line['Tumor_Seq_Allele2']
                        Effect = line['Variant_Classification_Alt']
                        if line['Chromosome'] == 'MT':
                            chrom = 'chrM'
                        else:
                            chrom = 'chr' + line['Chromosome']
                        pos = int(line['Start_Position'])
                        ref = line['Reference_Allele']
                        
                        
                        if inCDS('dogFast', chrom, pos):
                        
                            newVariant = variant(Amino_Acid_Change = Amino_Acid_Change, 
                                                 Gene_Name = Gene_Name, Exon_Rank = Exon_Rank, Transcript_ID = Transcript_ID, 
                                                 Genotype = Genotype, Effect = Effect, chrom = chrom, pos = pos, ref = ref, assembly = 'canFam3')

                            try:
                                sampleVariants[name].append(newVariant)
                            except KeyError:
                                sampleVariants[name] = [newVariant]

                for name in sampleVariants:
                    if len(sampleVariants[name]) < maxVar and len(sampleVariants[name]) >= minVar and name not in noFlyList and label != 'Remove':
                        samples.append(sample(variants = sampleVariants[name], name = name, species = 'dog', label = label, source = fileName.strip()))
                    elif name in noFlyList:
                        with open('noFly.out', 'a') as outFile:
                            outFile.write(name + '\n')
        
        elif 'vcf' in fileName:
            #print(fileName)
            name = os.path.basename(fileName)
            if 'LSA' in fileName:
                try:
                    label = LSALookup[name]
                except KeyError:
                    label = 'Remove'
            elif 'MELANOMA' in fileName:
                label = 'dogMELANOMA'
            elif 'HSA' in fileName:
                label = 'dogHSA'
            elif 'OSA' in fileName:
                label = 'dogOSA'
            elif 'CMT' in fileName:
                if fileName.split('/')[-1] in cmtLookup:
                    label = 'dogCMT'
                else:
                    label = 'skip'
            elif 'BCL' in fileName:
                label = 'dogBLSA'
            if label != 'skip':
                
                #with pysam.Fastafile(args.humanFasta) as humFast, pysam.Fastafile(args.dogFasta) as dogFast, pysam.Fastafile(args.hg19Fasta) as hg19Fast:
                with open(fileName, 'r') as vcf_in:
                    #print('here')
                    name = os.path.basename(fileName)
                    variants = []
                    for line in vcf_in:
                        Effect = ''
                        if line[0] == '#':
                            continue
                        fields = line.strip().split()
                        chrom = fields[0]
                        if 'chr' not in chrom:
                            chrom = 'chr' + chrom
                        ref = fields[3]
                        pos = int(fields[1])
                        Genotype = fields[4]
                        meta = fields[7].split(';')
                        snpEff = ''
                        for i in meta:
                            if i[:3] == 'EFF':
                                snpEff = i[4:].split(',')
                        Effect = snpEff[0].split('(')[0]
                        transcript = snpEff[0].split('|')[8]
                        try:
                            dogGene = dogdb[transcript]['gene_id'][0]
                        except:
                            dogGene = snpEff[0].split('|')[5]
                        if inCDS('dogFast', chrom, pos):
                            newVariant = variant(Transcript_ID = dogGene, Genotype = Genotype, Effect = Effect, chrom = chrom, pos = pos, ref = ref, assembly = 'canFam3')
                        #print('here3')
                            variants.append(newVariant)
                    if (len(variants) < maxVar and len(variants) >= minVar and label != 'Remove'):
                                    #print('here4')
                                    samples.append(sample(variants = variants, name = name, species = 'dog', label = label, source = fileName.strip()))
                        
                '''        
                with VariantFile(fileName) as vcf_in:
                    variants = []
                #split each record into its various transcripts (each record can list
                #effects for multiple transcripts
        
                    for rec in vcf_in.fetch():
                    
                        ref = rec.ref
                        chrom = rec.contig
                        if 'chr' not in chrom:
                            chrom = 'chr' + chrom
                        pos = int(rec.pos)
                        
                        
                        #snpEFF = re.split('[,\'()]', str(rec.info['EFF']))
                        
                        
                        #print('effect: ' + str(rec.info['EFF']))
                        transcripts = re.split('[,\'()]', str(rec.info['EFF']))
                        #print(transcripts)
                        #print(str(transcripts) + '\n\n\n')
                        
                    #split out fields from each transcript. Push to variants list
                        for transcript in transcripts:
                            #Effect = ''
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
                                #Gene_Name = fields[5]
                                Transcript_ID = fields[8]
                                try:
                                    Gene_Name = dogdb[Transcript_ID]['gene_id'][0]
                                except:
                                    Gene_Name = ''
                                #print('gene name ' + Gene_Name)
                                #print('Effect ' + Effect)
                                Transcript_Biotype = fields[6]
                                Gene_Coding = fields[7]
                                Exon_Rank = fields[9]
                                Genotype = fields[10]
                                if getHumanGeneID(Gene_Name) != 'NA':
                                    if inCDS('dogFast', chrom, pos):
                                        newVariant = variant(Effect_Impact, Functional_Class, 
                                                         Codon_change, Amino_Acid_Change, 
                                                         Amino_Acid_length, Gene_Name, 
                                                         Transcript_Biotype, Gene_Coding, 
                                                         Transcript_ID, Exon_Rank, 
                                                         Genotype, Effect, chrom, 
                                                         pos, ref, assembly = 'canFam3')
                                        variants.append(newVariant)
                                #print(Genotype)
                                #if len(fields[5]) > 0:
                                #    genes.add(fields[5])
                            if (len(variants) > maxVar):
                                break
                if (len(variants) < maxVar and len(variants) >= minVar and label != 'Remove'):
                    samples.append(sample(variants,name,'dog',label))
                '''
    
        elif 'tsv' in fileName:
            sampleVariants = {}
            #print(fileName)
            with gzip.open(fileName, 'rt') as ssm:
                try:
                    label = studyClass[fileName.split('.')[2]]
                except KeyError:
                    print('missing study class: ' + studyClass[fileName.split('.')[2]])
                    label = 'Remove'
                #label = studyClass[fileName.split('.')[2]]
                if label == 'Remove':
                    continue
                #label = label.split('-')[0]
                print(label)
                

                for line in ssm:
                    if line[:1] != '#':
                        names = line.strip().split('\t')
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
                    mutationID = row['icgc_mutation_id'].strip()
                    #mutationList.add(mutationID)
                    name = row['icgc_sample_id'].strip()
                    Effect = row['consequence_type'].strip()
                    Gene_Name = row['gene_affected'].strip()
                    Transcript_ID = row['transcript_affected'].strip()
                    Amino_Acid_Change = row['aa_mutation'].strip()
                    Codon_change = row['cds_mutation'].strip()
                    chrom = 'chr' + row['chromosome'].strip()
                    pos = int(row['chromosome_start'].strip())
                    ref = row['reference_genome_allele'].strip()
                    Genotype = row['mutated_to_allele'].strip()
                    #mutation type is something I don't think we have for dog, but seems important (ex. 'single base substitution') should look into getting for dog
                    #also need "Gene_Coding, Exon_Rank, functional class
                    mutationType = row['mutation_type'].strip()
                    
                    try:
                        test = sampleVariants[name]
                    except KeyError:
                        sampleVariants[name] = []
                    
                    if inCDS('hg19Fast', chrom, pos):
                        newVariant = variant(Codon_change = Codon_change, 
                                             Amino_Acid_Change = Amino_Acid_Change, 
                                             Gene_Name = Gene_Name, 
                                             Transcript_ID = Transcript_ID, 
                                             Genotype = Genotype, Effect = Effect, 
                                             chrom = chrom, pos = pos, ref = ref, 
                                             assembly = 'hg19')
                        try:
                            sampleVariants[name].append(newVariant)
                        except KeyError:
                            sampleVariants[name] = [newVariant]
                        
                for name in sampleVariants:
                    keep = True
                    if fileName.split('.')[2] == 'BOCA-UK':
                        print('mutation file name: ' + name)
                        if name not in BOCA_UK_lookup:
                            keep = False
                    #if fileName.split('.')[2] == 'ESCA':
                    #    keep = ESCA_lookup[name]
                    #if fileName.split('.')[2] == 'SARC':
                    #    keep = SARC_lookup[name]
                        
                    if len(sampleVariants[name]) < maxVar and len(sampleVariants[name]) >= minVar and name not in noFlyList and label != 'Remove':
                        samples.append(sample(variants = sampleVariants[name], name = name, species = 'Human', label = label, source = fileName.strip()))
                    elif name in noFlyList:
                        with open('noFly.out', 'a') as outFile:
                            outFile.write(name + '\n')
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

    
if args.triNucCon:
    with open(args.outFile + '_triCon.txt', 'w') as contextOutFile:
        
        contextOutFile.write('ID\tdisease')
        
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
        '''
        for context in contextKeys:
            featureList.append(context)
        '''
        
        for context in contextKeys:
            featureList.append(context)
            contextOutFile.write('\t' + context)
            
        contextOutFile.write('\n')
        
        for sample in samples:
            contextOutFile.write(sample.name + '\t' + sample.label)
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
            
            if not os.path.exists('triNucleotideContexts'):
                os.mkdir('triNucleotideContexts')
            assembly = sample.assembly
            with open('triNucleotideContexts/' + sample.name + '_' + assembly + '_tri.txt', 'w') as outFile:
                outFile.write('Sample\tRef\tAlt\tTrinuc\n')
            
                for variant in  sample.variants:
                    #if variant.variantID not in variantSet:
                    try:
                        if len(variant.refTrinucleotide) != 0 and variant.ref != variant.Genotype:
                            outFile.write(sample.name + '\t' + variant.ref + '\t' + variant.Genotype + '\t' + variant.refTrinucleotide + '\n')
                            triContextDict[variant.trinucleotide] += 1
                        #variantSet.add(variant.variantID)
                            #snpCount += 1
                    except KeyError:
                        if variant.trinucleotide != 'NA':
                            print('This ain\'t got no fitting in gud!')
                            print(variant.trinucleotide)
                            print(sample.name)
                            print(sample.label)
            #sample.VariantCount = len(variantSet)
            #print(sample.variants)
            mutCount = 0
            for context in contextKeys:
                mutCount += triContextDict[context]
                #print(mutCount)
            for context in contextKeys:
                try:
                    #print(str(triContextDict[context]/mutCount))
                    #removing the averaging so as to report the integer context counts for Diane               
                    contextOutFile.write('\t' + str(triContextDict[context]))
                    #rounding added to eliminate the possibility that the algorithm is figuring out how many
                    #samples in a class by the level of precision. Still possible with 2 digits, but less bad, probably
                    sample.toWrite.append(round(triContextDict[context]/mutCount,2))
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

if args.nonDamage:
    for geneID in nonDamageGenes:
        featureList.append(geneID + '_nonDamaging') 
        
    print('there are ' + str(len(samples)) + ' samples')
    for sample in samples:
            nonDamageInSample = []
            for geneID in nonDamageGenes:
                print('here?')
                if isGeneNonDamage(sample, geneID):
                    nonDamageInSample.append(1)
                    sample.toWrite.append(1)
                else:
                    nonDamageInSample.append(0)
                    sample.toWrite.append(0)


if args.brokenPathways:
    broked = brokenGenes
    
    brokenPathways = set()
    for gene in broked:
        try:
            for pathway in reference[gene].pathways:
                brokenPathways.add(pathway)
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

if args.driverMutCount:
    featureList.append('driverMutCount')
    for sample in samples:
        sample.toWrite.append(sample.driverMutCount)
        
if args.source:
    featureList.append('source')
    for sample in samples:
        sample.toWrite.append(sample.source)
        
        
if args.nonDamageGeneCount:
    for geneID in nonDamageGenes:
        featureList.append(geneID + '_nonDamaging_count')
        
    for sample in samples:
        for geneID in nonDamageGenes:             
            if isGeneNonDamage(sample, geneID):
                sample.toWrite.append(sample.nonDamageGeneCount[geneID])
            else:
                sample.toWrite.append(0)
            
if args.pathwayCount:
    brokenPathways = set()
    for gene in brokenGenes:
        try:
            for pathway in reference[gene].pathways:
                brokenPathways.add(pathway)
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
    learnOut.write('label\tID\tspecies')
    for feature in featureList:
        learnOut.write('\t' + feature)
    learnOut.write('\n')
    for sample in samples:
        learnOut.write(sample.label + '\t' + sample.name + '\t' + sample.species)
        for feature in sample.toWrite:
            learnOut.write('\t' + str(feature))
        learnOut.write('\n')


'''            
outVariantSet = list(outVariantSet)
outCodonSet = list(outCodonSet)
with open(args.outFile + 'codonsOnly.out', 'w') as codonOut:
    codonOut.write('label\tID\tspecies')
    for codon in outCodonSet:
        codonOut.write('\t' + codon)
    codonOut.write('\n')
    for sample in samples:
        codonOut.write(sample.label + '\t' + sample.name + '\t' + sample.species)
        for codon in outCodonSet:
            try:
                codonOut.write('\t' + str(sample.outCodons[codon]))
            except KeyError:
                codonOut.write('\t0')
        codonOut.write('\n')
        
with open(args.outFile + 'VariantsOnly.out', 'w') as varOut:
    varOut.write('label\tID\tspecies')
    for var in outVariantSet:
        varOut.write('\t' + var)
    varOut.write('\n')
    for sample in samples:
        varOut.write(sample.label + '\t' + sample.name + '\t' + sample.species)
        for var in outVariantSet:
            try:
                varOut.write('\t' + '/'.join(sample.outVariants[var]))
            except KeyError:
                varOut.write('\t0')
        varOut.write('\n')
    
    

with pysam.Fastafile(args.humanFasta) as humFast, pysam.Fastafile(args.dogFasta) as dogFast, pysam.Fastafile(args.hg19Fasta) as hg19Fast:
    with open(fileName, 'r') as vcf_in:
        name = os.path.basename(fileName)
        label = 'dogLSA'
        variants = []
        CDSCount = 0
        outCount = 0
        for line in vcf_in:
            Effect = ''
            if line[0] == '#':
                continue
            fields = line.strip().split()
            chrom = fields[0]
            if 'chr' not in chrom:
                chrom = 'chr' + chrom
            ref = fields[3]
            pos = int(fields[1])
            if not inCDS('dogFast', chrom, pos):
                outCount += 1
                continue
            else:
                CDSCount += 1
            Genotype = fields[4]
            meta = fields[7].split(';')
            snpEff = ''
            for i in meta:
                if i[:3] == 'EFF':
                    snpEff = i[4:].split(',')
            Effect = snpEff[0].split('(')[0]
            transcript = snpEff[0].split('|')[8]
            try:
                dogGene = dogdb[transcript]['gene_id'][0]
            except:
                dogGene = snpEff[0].split('|')[5]
            newVariant = variant(Transcript_ID = dogGene, Genotype = Genotype, Effect = Effect, chrom = chrom, pos = pos, ref = ref, assembly = 'canFam3')
            variants.append(newVariant)
        if (len(variants) < maxVar and len(variants) >= minVar and label != 'Remove'):
            samples.append(sample(variants,name,'dog',label))
        else:
            print(len(variants))

Ok, so we'll use the above. Need to check if variant in CDS, if so add, 
then check if effect is damaging and use the transcript id to look up the gene.
can have multipl effects per variant, so technically could have multiple genes with
damaging effects, though I don't think we handle that anywhere else, so *shrug*?
         
with VariantFile(fileName) as vcf_in:
                    variants = []
                #split each record into its various transcripts (each record can list
                #effects for multiple transcripts
        
                    for rec in vcf_in.fetch():
                    
                        ref = rec.ref
                        chrom = rec.contig
                        if 'chr' not in chrom:
                            chrom = 'chr' + chrom
                        pos = int(rec.pos)
                        
                        
                        #snpEFF = re.split('[,\'()]', str(rec.info['EFF']))
                        
                        
                        #print('effect: ' + str(rec.info['EFF']))
                        transcripts = re.split('[,\'()]', str(rec.info['EFF']))
                        #print(transcripts)
                        #print(str(transcripts) + '\n\n\n')
                        
                    #split out fields from each transcript. Push to variants list
                        for transcript in transcripts:
                            #Effect = ''
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
                                #Gene_Name = fields[5]
                                Transcript_ID = fields[8]
                                try:
                                    Gene_Name = dogdb[Transcript_ID]['gene_id'][0]
                                except:
                                    Gene_Name = ''
                                #print('gene name ' + Gene_Name)
                                #print('Effect ' + Effect)
                                Transcript_Biotype = fields[6]
                                Gene_Coding = fields[7]
                                Exon_Rank = fields[9]
                                Genotype = fields[10]
                                if getHumanGeneID(Gene_Name) != 'NA':
                                    if inCDS('dogFast', chrom, pos):
                                        newVariant = variant(Effect_Impact, Functional_Class, 
                                                         Codon_change, Amino_Acid_Change, 
                                                         Amino_Acid_length, Gene_Name, 
                                                         Transcript_Biotype, Gene_Coding, 
                                                         Transcript_ID, Exon_Rank, 
                                                         Genotype, Effect, chrom, 
                                                         pos, ref, assembly = 'canFam3')
                                        variants.append(newVariant)
                                #print(Genotype)
                                #if len(fields[5]) > 0:
                                #    genes.add(fields[5])
                            if (len(variants) > maxVar):
                                break
                if (len(variants) < maxVar and len(variants) >= minVar and label != 'Remove'):
                    samples.append(sample(variants,name,'dog',label))
                

inFile = 'codonTest120922codonsOnly.out'                
samples ={}
codonCounts = {} 
codonList = []
counts = []
header = ''         

with open(inFile, 'r') as codons:
    line = next(codons)
    codonLine = line.strip().split('\t')[3:]    
    for codon in codonLine:
        counts.append(0)
        codonCounts[codon] = 0
        codonList.append(codon)
    for line in codons:
        codonLine = line.strip().split('\t')[3:]
        for i in range(len(codonList)):
            codonCounts[codonList[i]] += int(codonLine[i])
        
        codonList.append(codon)
    for line in codons:
        codonLine = line.strip().split('\t')[3:]
    for i in range(len(codons)):
        counts[i] += int(codonLine[i])
    gc.collect()
        
while True:
    line = next(codons)
    codonLine = line.strip().split('\t')[3:]
    for i in range(len(codonLine)):
        counts[i] += int(codonLine[i])
with open(inFile, 'r') as inHeader:
    header = next(codons).strip().split('\t')

import gc 
keep = []
count = 0
with open('codonKeep.out', 'w') as outFile:
    for i in range(len(codonList)):
        codons = open(inFile, 'r')
        line = next(codons)
        for line in codons:
            codonLine = line.strip().split('\t')[3:]
            count += int(codonLine[i])
        if count > 5:
            outFile.write(codonList[i] + '\n')
            #keep.append(codonList[i])
        count = 0
        codons.close()

current = 0
lineLength = 1264455
with open('outFix.txt','w') as outFile:
    outFile.write('\t'.join(line[:lineLength]) + '\t' + line[lineLength][0] + '\n')
    current = lineLength
    while True:
        try:
            #outFile.write(line[current][1:] + '\t' + 
            outFile.write(line[current][1:] + '\t' + '\t'.join(line[(current + 1):(current + lineLength - 1)]) + '\t' + line[current + lineLength][0] + '\n')
            current = current + lineLength
        except IndexError:
            pass
line2 = line[current][0:] + line[current:current + lineLength]
line3 = line[2528910][0:] + line[]
with open('outFix.txt' ,'w') as outFile:


with open('codonTest120922codonsOnly.out', 'r') as inFile, open('codonTest5orMoreSamples', 'w') as outFile:
    for line in inFile:
        labs = np.array(line.strip().split('\t')[:3])
        subset = np.array(line.strip().split('\t')[3:])
        outLine = np.append(labs, subset[keepCols])
        outFile.write('\t'.join(outLine) + '\n')


fromPass = set()
docm = set()
file1 = 'fromPassengerHotspotMutationsInCancer_qLNP0.1.toHg38.VEPOutput.codonIntervals.bed'
file2 = 'docm_variants_hg19.toHg38.VEPoutPut.codonIntervals.bed'

with open(file1, 'r') as inFile:
    for line in inFile:
        fields = line.strip().split('\t')
        fromPass.add(fields[3])
with open(file2, 'r') as inFile:
    for line in inFile:
        fields = line.strip().split('\t')
        docm.add(fields[3][3:])
        
        
inFile = '/seq/vgb/swofford/R01/convertingToSKLInput/testVCFs/codonTest120922codonsOnly.out'
docm_variants_loc = []
fromPass_loc = []
hotspot_loc = []

with open(inFile, 'r') as codonFile:
    count = 0
    codons = next(codonFile).strip().split('\t')[3:]
    for i in codons:
        if i in fromPass:
            fromPass_loc.append(count)
            hotspot_loc.append(count)
        if i in docm:
            docm_variants_loc.append(count)
            if count not in hotspot_loc:
                hotspot_loc.append(count)
        count += 1


with open(inFile, 'r') as codonFile:
    count = 0
    codons = next(codonFile).strip().split('\t')[3:]
    for i in codons:
        if i in fromPass:
            fromPass_loc.append(count)
            hotspot_loc.append(count)
        if i in docm:
            docm_variants_loc.append(count)
            if count not in hotspot_loc:
                hotspot_loc.append(count)
        count += 1


with open('/seq/vgb/swofford/R01/convertingToSKLInput/testVCFs/codonTest120922codonsOnly.out', 'r') as inFile, open('fromPassAndDocmCodons.txt', 'w') as allOut, open('fromPassCodons.txt', 'w') as fromOut, open('docmCodons.txt', 'w') as docmOut:
    for line in inFile:
        labs = np.array(line.strip().split('\t')[:3])
        subset = np.array(line.strip().split('\t')[3:])
        allOutLine = np.append(labs, subset[hotspot_loc])
        allOut.write('\t'.join(allOutLine) + '\n')
        fromLine = np.append(labs, subset[fromPass_loc])
        fromOut.write('\t'.join(fromLine) + '\n')
        docmLine = np.append(labs, subset[docm_variants_loc])
        docmOut.write('\t'.join(docmLine) + '\n')
        
with open('/seq/vgb/swofford/R01/convertingToSKLInput/testVCFs/codonTest120922codonsOnly.out', 'r') as inFile, open('codonsOnly_5SamplesPlus_fromPass_docm.txt', 'w') as outFile:
    for line in inFile:
        labs = np.array(line.strip().split('\t')[:3])
        subset = np.array(line.strip().split('\t')[3:])
        outLine = np.append(labs, subset[keep_loc])
        outFile.write('\t'.join(outLine) + '\n')
'''    
        
        