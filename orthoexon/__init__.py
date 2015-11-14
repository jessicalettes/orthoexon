# -*- coding: utf-8 -*-

__author__ = 'Jessica Lettes'
__email__ = 'jlettes@ucsd.edu'
__version__ = '0.1.0'


#imports
#for gffutils database
import gffutils

#for pandas data frame
import pandas as pd

#import Seq for biopython translation
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

#to write files
import csv

#for FASTA file
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


#import gffutils databases and Ensembl compare data as a table
#species1DB = gffutils.FeatureDB('/Users/rhythmicstar/projects/exon_evolution'
#                                '//gencode.v19.annotation.gtf.db',
#                                keep_order=True)
#species1DB = gffutils.FeatureDB('/Users/rhythmicstar/projects/exon_evolution'
#                                '//gencode.v19.annotation.chrX.gtf.db',
#                                keep_order=True)
species1DB = gffutils.FeatureDB('/Users/rhythmicstar/projects/exon_evolution//'
                             'gencode.v19.annotation.humanrbfox2andfmr1.gtf.db'
                             , keep_order=True)
#species2DB = gffutils.FeatureDB('/Users/rhythmicstar/projects/exon_evolution'
#                                '//gencode.vM5.annotation.gtf.db',
#                                keep_order=True)
#species2DB = gffutils.FeatureDB('/Users/rhythmicstar/projects/exon_evolution'
#                                '//gencode.vM5.annotation.chrX.gtf.db',
#                                keep_order=True)
species2DB = gffutils.FeatureDB('/Users/rhythmicstar/projects/exon_evolution//'
                             'gencode.v19.annotation.mouserbfox2andfmr1.gtf.db'
                             , keep_order=True)
compara = pd.read_table('/Users/rhythmicstar/projects/exon_evolution/'
                        'EnsemblHumanMouse.txt')

#Fasta files with sequence
species1Fasta ='/Users/rhythmicstar/projects/exon_evolution//' \
               'GRCh37.p13.genome.fa'
species2Fasta ='/Users/rhythmicstar/projects/exon_evolution//' \
               'GRCm38.p3.genome.fa'




#to change gencode gene ids into ensembl gene ids
def separate(geneId):
    sep = '.'
    splitGeneId = geneId.partition(sep)
    splittingGeneId = splitGeneId[0]
    sep2 = "'"
    splitGeneId2 = splittingGeneId.partition(sep2)
    newGeneId = splitGeneId2[2]
    return newGeneId


# to translate sequence with correct strand and frame
def translate(exon, fasta):
    exonFrame = int(exon.frame)
    exonSeq = exon.sequence(fasta, use_strand=False)
    exonSeq = Seq(exonSeq, alphabet = generic_dna)
    if exon.strand == '-':
        exonSeq = exonSeq.reverse_complement()
    exonProtein = exonSeq[exonFrame:].translate(to_stop=True)
    #print("{}:{}-{}:{} {}\t{}\t{}".format(exon.chrom, exon.start, exon.stop,
    #                                      exon.frame, exon.strand, exonSeq,
    #                                      exonProtein))
    return exonProtein


# to make an array of all exons to make FASTA file
def transferarray(proteindfsize, finalproteindf):
    transferArray = []
    for x in range (0, int((proteindfsize/4))):
        transferArray.append(x)
        transferArray[x] = SeqRecord(Seq(finalproteindf.iat[x,3]),
                                     finalproteindf.iat[x,1], description='')
    return transferArray





def orthoexon():
    #def orthoexon(species1DB, species2DB, compara, species1Fasta,
    # species2Fasta):
    #Drop compara duplicates
    dropDuplicates = compara.drop_duplicates(['Ensembl Gene ID',
                                          'Ensembl Gene ID.1'])
    newCompara = dropDuplicates.reset_index()
    comparaGeneIdIndex = newCompara.set_index('Ensembl Gene ID')
    numcompara = int((newCompara.size/5))


    #dataframe for proteins
    proteindf = pd.DataFrame(columns=['Proteins', 'Gene Id', 'Exon Id'])

    #for each human gene in gffutils database get gene id
    for index, species1gene in enumerate(species1DB.features_of_type('gene')):
        species1GFFUtilsGeneId = str(species1gene['gene_id'])
        species1geneid = separate(species1GFFUtilsGeneId)
        #if(index % 10 == 0):
        print(index)

        #if gene ID equals one from ensembl, get mouse gene ID & exons at point
        try:
            species1EnsGeneId = comparaGeneIdIndex.loc[species1geneid]
        except KeyError:
            continue

        for x in range (0, numcompara):
            if (species1geneid == newCompara.iat[x, 1]):
                species2EnsGeneId = newCompara.iat[x, 3]
        #get human exons
        for species1Exon in species1DB.children(species1gene,
                                                featuretype = 'CDS',
                                                order_by = 'start'):
            #print(species1Exon.chrom)
            species1ExonProtein = translate(species1Exon, species1Fasta)

            #create table with protein seq, geneId and exon Id
            appendproteindf = pd.DataFrame({'Proteins' :
                                                [str(species1ExonProtein)],
                                            'Gene Id' : [str(species1geneid)],
                                            'Exon Id' : [str(species1Exon
                                                             ['exon_id']
                                                         [0])]})
            proteindf = proteindf.append(appendproteindf, ignore_index=True)



        #CHECK THIS!
        #for each mouse gene ID from database
        for species2gene in species2DB.features_of_type('gene'):
            species2GFFUtilsGeneId = str(species2gene['gene_id']) #gffutils id
            species2geneid = separate(species2GFFUtilsGeneId)
            if (species2EnsGeneId == species2geneid):
                for species2Exon in species2DB.children(species2gene,
                                                    featuretype = 'CDS',
                                                    order_by = 'start'):
                    species2ExonProtein = translate(species2Exon,
                                                    species2Fasta)

                    #create dataframe with protein seq, geneId and exon Id
                    appendproteindf=pd.DataFrame({'Proteins'
                                                  : [str(species2ExonProtein)],
                                                'Gene Id' :
                                                     [str(species2geneid)],
                                                'Exon Id' : [str(species2Exon
                                                                    ['exon_id']
                                                                    [0])]})

                    proteindf = proteindf.append(appendproteindf,
                                                 ignore_index=True)


        #drop duplicates for protein seq
        newproteindf = proteindf.drop_duplicates('Exon Id')
        finalproteindf = newproteindf.reset_index()
        arraytransfer = transferarray(finalproteindf.size, finalproteindf)

        #Write to FASTA
        output_handle = open('Human Gene {} - Mouse Gene {}.fasta'.format
                             (str(species1geneid),str(species2EnsGeneId)), "w")
        SeqIO.write(arraytransfer, output_handle, "fasta")
        output_handle.close()

        #Reset the protein dataframe
        proteindf = pd.DataFrame(columns=['Proteins', 'Gene Id', 'Exon Id'])