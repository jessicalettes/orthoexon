__author__ = 'rhythmicstar'

from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

def separate(geneId):
    sep = '.'
    splitGeneId = geneId.partition(sep)
    splittingGeneId = splitGeneId[0]
    splittingGeneId = splittingGeneId.strip("'")
    splittingGeneId = splittingGeneId.strip('"')
    return splittingGeneId


def splitstart(location):
    sep = ':'
    splitstart = location.partition(sep)
    splittingstart = splitstart[2]
    sep2 = "-"
    splitstart2 = splittingstart.partition(sep2)
    startLoc = splitstart2[0]
    return startLoc


def splitend(location):
    sep = '-'
    splitend = location.partition(sep)
    splittingend = splitend[2]
    sep2 = ":"
    splitend2 = splittingend.partition(sep2)
    endLoc = splitend2[0]
    return endLoc

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


# to change gencode gene ids into ensembl gene ids
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


def make_sequence_array(finalsequencedf):
    sequence_array = []
    for index, row in finalsequencedf.iterrows():
            sequence_array.append(SeqRecord(Seq(row['Sequences']),
                                            id=row['Exon ID'], description=''))
    return sequence_array