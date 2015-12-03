# -*- coding: utf-8 -*-

__author__ = 'Jessica Lettes'
__email__ = 'jlettes@ucsd.edu'
__version__ = '0.1.0'

import sys

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
species1DB = gffutils.FeatureDB('/Users/rhythmicstar/projects/exon_evolution//'
                             'gencode.v19.annotation.humanrbfox2andfmr1.gtf.db'
                             , keep_order=True)
#species2DB = gffutils.FeatureDB('/Users/rhythmicstar/projects/exon_evolution'
#                                '//gencode.vM5.annotation.gtf.db',
#                                keep_order=True)
species2DB = gffutils.FeatureDB('/Users/rhythmicstar/projects/exon_evolution//'
                             'gencode.vM5.annotation.mouserbfox2andfmr1.gtf.db'
                             , keep_order=True)
compara = pd.read_table('/Users/rhythmicstar/projects/exon_evolution/'
                        'EnsemblHumanMouse.txt')

#Fasta files with sequence
#species1Fasta ='/Users/rhythmicstar/projects/exon_evolution//' \
#               'GRCh37.p13.genome.fa'
species1Fasta = '/Users/rhythmicstar/projects/exon_evolution//GRCh37.p13.' \
                'genome.x22.fa'

#species2Fasta ='/Users/rhythmicstar/projects/exon_evolution//' \
#               'GRCm38.p3.genome.fa'
species2Fasta = '/Users/rhythmicstar/projects/exon_evolution//GRCm38.p3.' \
                'genome.x15.fa'



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
    print("{}:{}-{}:{} {}\t{}\t{}".format(exon.chrom, exon.start, exon.stop,
                                          exon.frame, exon.strand, exonSeq,
                                          exonProtein))
    return exonProtein


# to get sequence with correct strand and frame
def getsequence(exon, fasta):
    exonFrame = int(exon.frame)
    exonSeq = exon.sequence(fasta, use_strand=False)
    exonSeq = Seq(exonSeq, alphabet = generic_dna)
    if exon.strand == '-':
        exonSeq = exonSeq.reverse_complement()
    exonSeq = exonSeq[exonFrame:]
    print("{}:{}-{}:{} {}\t{}".format(exon.chrom, exon.start, exon.stop,
                                      exon.frame, exon.strand, exonSeq))
    return exonSeq


# to make an array of all exons to make FASTA file
def make_sequence_array(finalsequencedf):
    sequence_array = []
    for index, row in finalsequencedf.iterrows():
            sequence_array.append(SeqRecord(Seq(row['Sequences']), 
                                            id=row['Exon ID'], description=''))
    return sequence_array


def orthoexon(species1name, species2name, species1DB, species2DB, compara,
              species1Fasta, species2Fasta, translateFlag):
    #Drop compara duplicates
    dropDuplicates = compara.drop_duplicates(['Ensembl Gene ID',
                                              'Ensembl Gene ID.1'])
    newCompara = dropDuplicates.reset_index()
    comparaGeneIdIndex = newCompara.set_index('Ensembl Gene ID')
    numcompara = int((newCompara.size/5))

    gene_dfs = []

    #dataframe for proteins
    data = []

    #for each human gene in gffutils database get gene id
    for index, species1gene in enumerate(species1DB.features_of_type('gene')):
        species1GFFUtilsGeneId = str(species1gene['gene_id'])
        species1geneid = separate(species1GFFUtilsGeneId)
        #if(index % 10 == 0):
        print(index)

        #if gene ID equals one from ensembl, get mouse gene ID & exons at point
        try:
            species1ensgeneid = comparaGeneIdIndex.loc[species1geneid]
        except KeyError:
            continue

        for x in range (0, numcompara):
            if (species1geneid == newCompara.iat[x, 1]):
                species2EnsGeneId = newCompara.iat[x, 3]
        #get human exons
        for species1Exon in species1DB.children(species1gene,
                                                featuretype = 'CDS',
                                                order_by = 'start'):
            if (translateFlag == True):
                species1ExonProtein = translate(species1Exon, species1Fasta)
                row = [species1name, str(species1ExonProtein),
                       str(species1geneid), str(species1Exon['exon_id'][0])]
                data.append(row)
            else:
                species1ExonSeq = getsequence(species1Exon, species1Fasta)
                row = [species1name, str(species1ExonSeq), str(species1geneid),
                       str(species1Exon['exon_id'][0])]
                data.append(row)

        #CHECK THIS!
        #for each mouse gene ID from database
        for species2gene in species2DB.features_of_type('gene'):
            species2GFFUtilsGeneId = str(species2gene['gene_id']) #gffutils id
            species2geneid = separate(species2GFFUtilsGeneId)
            if (species2EnsGeneId == species2geneid):
                for species2Exon in species2DB.children(species2gene,
                                                    featuretype = 'CDS',
                                                    order_by = 'start'):

                    if (translateFlag == True):
                        species2ExonProtein = translate(
                            species2Exon, species2Fasta)
                        row = [species1name, str(species2ExonProtein),
                               str(species2geneid),
                               str(species2Exon['exon_id'][0])]
                        data.append(row)

                    else:
                        species2ExonSeq = getsequence(species2Exon,
                                                      species2Fasta)
                        row = [species2name, str(species2ExonSeq), str(
                            species2geneid), str(species2Exon['exon_id'][0])]
                        data.append(row)

        sequencedf = pd.DataFrame(data, columns=['Species', 'Sequences',
                                                'Gene ID', 'Exon ID'])

    #drop duplicates for protein seq
    sequencedf_noduplicates = sequencedf.drop_duplicates('Sequences')

    sequencearray = make_sequence_array(sequencedf_noduplicates)

    #Write to FASTA
    if (translateFlag == True):
        output_filename = '{}_{}_Proteins.fasta'.format(species1name,
                                                        species2name)
    else:
        output_filename = '{}_{}_Sequences.fasta'.format(species1name,
                                                        species2name)
    sys.stdout.write('Writing sequences to "{}" ...\n'.format(output_filename))
    output_handle = open(output_filename, "w")
    SeqIO.write(sequencearray, output_handle, "fasta")
    output_handle.close()
    sys.stdout.write('\tDone.\n')

    gene_dfs.append(sequencedf)
    all_gene_dfs = pd.concat(gene_dfs, ignore_index=True)
    all_gene_dfs = all_gene_dfs.drop_duplicates()
    return all_gene_dfs

    #savesequencedf.to_csv("AllExons.csv", columns= ['Reset Index', 'Proteins',
    #'Gene Id', 'Exon Id'])
    # savedroppedsequencedf.to_csv("AlldroppedExons.csv", columns=
    # ['Reset Index', 'Proteins', 'Gene Id', 'Exon Id'])