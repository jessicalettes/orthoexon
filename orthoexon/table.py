
__author__ = 'rhythmicstar'

import copy

import gffutils
import pandas as pd
import numpy as np

from .util import separate, separate, splitstart, splitend

species1DB = gffutils.FeatureDB('/Users/rhythmicstar/projects/exon_evolution//'
                                'gencode.v19.annotation.humanrbfox2and'
                                'fmr1andsnap25.gtf.db', keep_order=True)
species2DB = gffutils.FeatureDB('/Users/rhythmicstar/projects/exon_evolution//'
                                'gencode.vM5.annotation.mouserbfox2andfmr1and'
                                'snap25.gtf.db', keep_order=True)


class OrthologyTable(object):

    exon1 = 'exon1'
    exon2 = 'exon2'

    def __init__(self, nucleotide_table, protein_table, species1_name,
                 species1_version, species2_name, species2_version):
        self.nucleotide = self.read_blast_table(nucleotide_table, 'nucleotide')
        self.protein = self.read_blast_table(protein_table, 'protein')

    def read_blast_table(self, filename, sequence='protein'):
        """Parse blast alignment table, with only non-identical comparisons

        Parameters
        ----------
        filename : str
            Path of the blast output
        sequence : "protein" | "nucleotide"
            Specifies whether the table being read is from nucleotide or
            protein comparisons, which modifies the column names of the output
            pandas DataFrame

        Output
        ------
        table : pandas.DataFrame
            Parsed table of the exon alignments
        """

        abbrev = 'prot' if sequence == 'protein' else "nuc"

        columns = [self.exon1, self.exon2, '{}_length_of_overlap'.format(abbrev),
                   '{}_PID'.format(abbrev)]
        table = pd.read_table(filename, names=columns)

        # Get only rows where exon1 and exon2 don't have the same ID
        different = table[self.exon1] != table[self.exon2]
        table = table.loc[different]
        return table

    def remove_duplicate_comparisons(self, table):
        """Check for the same comparison in the opposite order

        E.g. (exonA, exonB) is the same as (exonB, exonA), so we want to
        remove that

        Parameters
        ----------
        table : pandas.DataFrame
            Modified blast orthology table

        Output
        ------
        table: pandas.DataFrame
            Parsed table with same comparison in opposite order removed
        """
        seen = set([])
        rows_to_use = []

        for i, row in table.iterrows():
            exon1 = row[self.exon1]
            exon2 = row[self.exon2]
            pair = tuple(sorted([exon1, exon2]))

            if pair not in seen:
                seen.update({pair})
                rows_to_use.append(i)

        table = table.loc[rows_to_use]
        table = table.reset_index(drop=True)
        return table

    def blast_table_setup(self, ortho_table):
        ref_blast_table = copy.deepcopy(ortho_table)
        rows = int(ref_blast_table.size/7)

        for row in range (0, rows):
            ref_blast_table.ix[(row + rows), 0] = ref_blast_table.ix[row, 1]

        cols = [1,2,3]
        ref_blast_table.drop(ref_blast_table.columns[cols],axis=1,inplace=True)
        ref_blast_table = ref_blast_table.drop_duplicates('Exon')
        ref_blast_table = ref_blast_table.reset_index(drop=True)
        return ref_blast_table

    def add_ortho_table_columns(self, nucleotide_table, protein_table,
                                rows, protrow):
        """add Prot_PID, Prot_Length_of_Overlap to ortho_table

        Parameters
        ----------
        nucleotide_table : pandas.DataFrame
            Modified blastn orthology table
        protein_table : pandas.DataFrame
            Modified blastp orthology table

        Output
        ------
        ortho_table: pandas.DataFrame
            Orthology table with rows from both protein and nucleotide blast
        """
        ortholog_columns = pd.DataFrame(columns=['Prot_PID'])
        nucleotide_table.insert(3, 'Prot_Length_of_Overlap', 'Below Threshold of ...')
        nucleotide_table = pd.concat([nucleotide_table, ortholog_columns], axis=1)
        nucleotide_table = nucleotide_table.replace(np.nan,'Below Threshold of ...', regex=True)

        for rowN in range (0, rows):
            for rowP in range (0, protrow):
                exon1N = nucleotide_table.iat[rowN, 0]
                exon1P = protein_table.iat[rowP, 0]
                exon2N = nucleotide_table.iat[rowN, 1]
                exon2P = protein_table.iat[rowP, 1]
                if ((exon1N == exon1P) & (exon2N == exon2P)):
                    nucleotide_table.iat[rowN, 5] = protein_table.iat[rowP, 3]
                    nucleotide_table.iat[rowN, 3] = protein_table.iat[rowP, 2]
                    protein_table.ix[rowP, 0] = np.nan
                    break

        protein_table.dropna(axis=0,inplace=True)
        protein_table = protein_table.reset_index(drop=True)

        protein_table.insert(2, 'Nuc_Length_of_Overlap', 'Below Threshold of ...')
        protein_table.insert(4, 'Nuc_PID', 'Below Threshold of ...')

        ortho_table = pd.concat([nucleotide_table, protein_table])
        ortho_table = ortho_table.reset_index(drop=True)
        return ortho_table


    def fill_in_blast_table(self, exon, speciesname, speciesversion):
        GFFUtilsExonId = str(exon['exon_id']) #gffutils id
        GFFUtilsExonId = separate(GFFUtilsExonId)
        if (GFFUtilsExonId == exon_ID2):
            exon_geneid = str(exon.attributes['gene_id'])
            exon_geneid = separate(exon_geneid)

            exon_transcriptid = str(exon.attributes['transcript_id'])
            exon_transcriptid = separate(exon_transcriptid)

            new_blast_row = [exon_ID, ("{}{}{}{}").format(
                speciesname,'(', speciesversion, ')'), (
                "{}:{}-{}:{}:{}".format(exon.chrom, exon.start,
                                        exon.stop, exon.strand,
                                        exon.frame)), exon_geneid,
                             exon_transcriptid, exon.stop -
                             exon.start + 1]
        return new_blast_row



    # figure out how to loop and get location
    #  def location(self):
    #     ("{}:{}-{}:{}:{}".format(exon.chrom, exon.start, exon.stop,
    #                              exon.strand, exon.frame)

    # fix this!!! get rid of blast and call location
    #  def orthology(self, ortho_table, blast_table, start1 = '0', start2 = '1',
    #               end1 = '0', end2 = '1'):
    #      """add orthology classification to ortho_table
    #
    #     Parameters
    #     ----------
    #     ortho_table : pandas.DataFrame
    #         Modified blast orthology table
    #     blast_table :
    #         Table with
    #
    #     Output
    #     ------
    #     table: pandas.DataFrame
    #         Orthology table with added columns for protein pid and overlap
    #     """
    #
    #     if (((ortho_table.iat[row, 0])[0:7]) != ((ortho_table.iat[row, 1])[0:7])):
    #         ortho_table.ix[row, 'Relationship'] = 'Orthologous'
    #     else:
    #         for inrow in range (0, blastRows):
    #             if (ortho_table.iat[row, 0] == blast_table.iat[inrow, 0]):
    #                 start1 = splitstart(blast_table.iat[inrow, 2])
    #                 end1 = splitend(blast_table.iat[inrow, 2])
    #             if (ortho_table.iat[row, 1] == blast_table.iat[inrow, 0]):
    #                 start2 = splitstart(blast_table.iat[inrow, 2])
    #                 end2 = splitend(blast_table.iat[inrow, 2])
    #             if (((start1 >= start2) and (end2 >= end1)) or
    #                 ((start2 >= start1) and (end1 >= end2)) or
    #                 ((start2 >= start1) and (end2 >= end1) and
    #                  (start2 <= end1)) or
    #                 ((start1 >= start2) and (end1 >= end2) and
    #                  (start1 <= end2))):
    #                 ortho_table.iat[row, 6] = 'Overlapping Genomic Loci'
    #                 if (blast_table.iat[inrow, 2] == blast_table.iat[inrow, 2]):
    #                     ortho_table.iat[row, 6] = 'Identical Genomic Loci'
    #             else:
    #                 ortho_table.iat[row, 6] = 'Paralogous'
    #     return ortho_table




    def save_to_csv(self, table, table_type):
        """Save table to csv format

        Parameters
        ----------
        table_type : type of table being saved
            Either ortho or blast
        """
        blast_columns = ['Exon', 'Species(Version)',
                         'Chrom:Start-Stop:Strand:Offset', 'Gene', 'Transcript',
                         'Exon_Length']
        ortho_columns = ['Exon', 'Exon2', 'Nuc_Length_of_Overlap',
                         'Prot_Length_of_Overlap', 'Nuc_PID', 'Prot_PID',
                         'Relationship']
        if table_type == 'blast':
            filename = "BLAST_Table.csv"
            column = blast_columns
        else:
            filename = "Ortho_Table.csv"
            column = ortho_columns
        table.to_csv(filename, columns= column, index=False)






def create_ortho_table(species1name, species2name, species1version,
                       species2version, blastnucfilename, blastprotfilename):
    read_blast_table(self, blastnucfilename, sequence='nucleotide')
    read_blast_table(self, blastprotfilename, sequence='protein')

    remove_duplicate_comparisons(self, Ortho_Table)
    remove_duplicate_comparisons(self, Prot_Table)

    #add the other columns to the Ortho_Table
    Ortholog_Columns = pd.DataFrame(columns=['Prot_PID'])
    #Overlap_Column = pd.DataFrame(columns=['Prot_Length_of_Overlap'])
    Ortho_Table.insert(3, 'Prot_Length_of_Overlap', 'Below Threshold of ...')

    #append new dataframe to Ortho_Table
    Ortho_Table = pd.concat([Ortho_Table, Ortholog_Columns], axis=1)
    Ortho_Table = Ortho_Table.replace(np.nan,'Below Threshold of ...', regex=True)

    BLAST_Table = []

    #to drop a row that has the same two exons
    rows = int(Ortho_Table.size/6)
    #loop through each row
    for row in range (0, rows):
        if(Ortho_Table.iat[row, 0] == Ortho_Table.iat[row, 1]):
            Ortho_Table.ix[row, 0] = np.nan
    #remove each duplicate
    Ortho_Table.dropna(axis=0,inplace=True)
    Ortho_Table = Ortho_Table.reset_index(drop=True)

    protrow = int(Prot_Table.size/4)
    #loop through each row
    for row in range (0, protrow):
        if(Prot_Table.iat[row, 0] == Prot_Table.iat[row, 1]):
            Prot_Table.ix[row, 0] = np.nan

    #remove each duplicate
    Prot_Table.dropna(axis=0,inplace=True)

    Prot_Table = Prot_Table.reset_index(drop=True)

    #to remove rows in nuc table that have the same two exons in the opposite order
    seen = set([])
    rows_to_use = []

    for i, row in Ortho_Table.iterrows():
        exon1 = row['Exon']
        exon2 = row['Exon2']
        pair = tuple(sorted([exon1, exon2]))

        if pair not in seen:
            seen.update({pair})
            rows_to_use.append(i)

    ortho_table_no_duplicate_comparisons = Ortho_Table.loc[rows_to_use]
    ortho_table_no_duplicate_comparisons = ortho_table_no_duplicate_comparisons.reset_index(drop=True)
    Ortho_Table = ortho_table_no_duplicate_comparisons

    #to remove rows in protein table that have the same two exons in the opposite order
    seen = set([])
    rows_to_use = []

    for i, row in Prot_Table.iterrows():
        exon1 = row['Exon']
        exon2 = row['Exon2']
        pair = tuple(sorted([exon1, exon2]))

        if pair not in seen:
            seen.update({pair})
            rows_to_use.append(i)

    prot_table_no_duplicate_comparisons = Prot_Table.loc[rows_to_use]
    prot_table_no_duplicate_comparisons = prot_table_no_duplicate_comparisons.reset_index(drop=True)
    Prot_Table = prot_table_no_duplicate_comparisons

    rows = int(Ortho_Table.size/6)
    protrow = int(Prot_Table.size/4)

    for rowN in range (0, rows):
        for rowP in range (0,protrow):
            exon1N = Ortho_Table.iat[rowN, 0]
            exon1P = Prot_Table.iat[rowP, 0]
            exon2N = Ortho_Table.iat[rowN, 1]
            exon2P = Prot_Table.iat[rowP, 1]
            if ((exon1N == exon1P) & (exon2N == exon2P)):
                #put protein pid in nuc table
                Ortho_Table.iat[rowN, 5] = Prot_Table.iat[rowP, 3]
                Ortho_Table.iat[rowN, 3] = Prot_Table.iat[rowP, 2]
                #put na into protein table
                Prot_Table.ix[rowP, 0] = np.nan
                break

    #remove each duplicate
    Prot_Table.dropna(axis=0,inplace=True)
    Prot_Table = Prot_Table.reset_index(drop=True)

    # #in Prot_Table, add col with 'Nuc_PID' between 'Length_of_Overlap' and 'Prot_PID'
    Prot_Table.insert(2, 'Nuc_Length_of_Overlap', 'Below Threshold of ...')
    Prot_Table.insert(4, 'Nuc_PID', 'Below Threshold of ...')

    #add protein table at end
    Ortho_Table = pd.concat([Ortho_Table, Prot_Table])
    Ortho_Table = Ortho_Table.reset_index(drop=True)

    #create new dataframe
    Relationship_Column = pd.DataFrame(columns=['Relationship'])

    #append new dataframe to Ortho_Table
    Ortho_Table = pd.concat([Ortho_Table, Relationship_Column], axis=1)

    #create BLAST_Table dataframe with exon ids
    Ref_BLAST_Table = copy.deepcopy(Ortho_Table)
    rows = int(Ref_BLAST_Table.size/7)

    for row in range (0, rows):
        Ref_BLAST_Table.ix[(row + rows), 0] = Ref_BLAST_Table.ix[row, 1]

    cols = [1,2,3,4,5,6]
    Ref_BLAST_Table.drop(Ref_BLAST_Table.columns[cols],axis=1,inplace=True)
    Ref_BLAST_Table = Ref_BLAST_Table.drop_duplicates('Exon')
    Ref_BLAST_Table = Ref_BLAST_Table.reset_index(drop=True)

    #in form to do
    #get data to fill in table
    rows = int(Ref_BLAST_Table.size)

    for exonCode in species1DB.features_of_type('CDS'):
        Species1Code = str(exonCode['exon_id'])
        Species1Code = separate(Species1Code)
        Species1Code = Species1Code[0:7]
        break

    #loop through each exon and get its gene and length
    for row in range (0, rows):
        exon_ID = Ref_BLAST_Table.ix[row, 'Exon']
        exon_ID2 = separate(exon_ID)

        #determine which file to look in
        if (exon_ID2[0:7] == Species1Code):
            for exon in species1DB.features_of_type('CDS'):
                GFFUtilsExonId = str(exon['exon_id']) #gffutils id
                GFFUtilsExonId = separate(GFFUtilsExonId)
                if (GFFUtilsExonId == exon_ID2):
                    exon_geneid = str(exon.attributes['gene_id'])
                    exon_geneid = separate(exon_geneid)

                    exon_transcriptid = str(exon.attributes['transcript_id'])
                    exon_transcriptid = separate(exon_transcriptid)

                    New_BLAST_Row = [exon_ID, ("{}{}{}{}").format(species1name,'(', species1version, ')'),
                                     ("{}:{}-{}:{}:{}".format(exon.chrom, exon.start, exon.stop,
                                                              exon.strand, exon.frame)), exon_geneid,
                                     exon_transcriptid, exon.stop - exon.start + 1]

                    BLAST_Table.append(New_BLAST_Row)

        else:
            #if exon is mouse
            for exon in species2DB.features_of_type('CDS'):
                GFFUtilsExonId = str(exon['exon_id']) #gffutils id
                GFFUtilsExonId = separate(GFFUtilsExonId)
                if (GFFUtilsExonId == exon_ID2):
                    exon_geneid = str(exon.attributes['gene_id'])
                    exon_geneid = separate(exon_geneid)

                    exon_transcriptid = str(exon.attributes['transcript_id'])
                    exon_transcriptid = separate(exon_transcriptid)

                    New_BLAST_Row = [exon_ID, ("{}{}{}{}").format(species2name,'(', species2version, ')'),
                                     ("{}:{}-{}:{}:{}".format(exon.chrom, exon.start, exon.stop,
                                                              exon.strand, exon.frame)), exon_geneid,
                                     exon_transcriptid, exon.stop - exon.start + 1]

                    BLAST_Table.append(New_BLAST_Row)

    BLAST_Table = pd.DataFrame(BLAST_Table)
    BLAST_Table.columns = ['Exon', 'Species(Version)', 'Chrom:Start-Stop:Strand:Offset', 'Gene', 'Transcript', 'Exon_Length']

    #fill in row for paralogous or orthologous
    rows = int(Ortho_Table.size/7)
    blastRows = int(BLAST_Table.size/6)
    for row in range (0, rows):
    #paralogous
        if (((Ortho_Table.iat[row, 0])[0:7]) != ((Ortho_Table.iat[row, 1])[0:7])):
            Ortho_Table.ix[row, 'Relationship'] = 'Orthologous'

    #orthologous
        else:
            exonOne = Ortho_Table.iat[row, 0]
            exonTwo = Ortho_Table.iat[row, 1]
            locationOne = 0
            locationTwo = 1
            start1 = '0'
            start2 = '1'
            end1 = '0'
            end2 = '1'
            for innerrow in range (0, blastRows):
                if (exonOne == BLAST_Table.iat[innerrow, 0]):
                    locationOne = BLAST_Table.iat[innerrow, 2]
                    start1 = splitstart(locationOne)
                    end1 = splitend(locationOne)
                if (exonTwo == BLAST_Table.iat[innerrow, 0]):
                    locationTwo = BLAST_Table.iat[innerrow, 2]
                    start2 = splitstart(locationTwo)
                    end2 = splitend(locationTwo)
                if (((start1 >= start2) and (end2 >= end1)) or
                    ((start2 >= start1) and (end1 >= end2)) or
                    ((start2 >= start1) and (end2 >= end1) and
                     (start2 <= end1)) or
                    ((start1 >= start2) and (end1 >= end2) and
                     (start1 <= end2))):
                    Ortho_Table.iat[row, 6] = 'Overlapping Genomic Loci'
                    if (locationOne == locationTwo):
                        Ortho_Table.iat[row, 6] = 'Identical Genomic Loci'
                else:
                    Ortho_Table.iat[row, 6] = 'Paralogous'

    #to save BLAST_Table
    BLAST_Table.to_csv("BLAST_Table.csv", columns= ['Exon', 'Species(Version)', 'Chrom:Start-Stop:Strand:Offset',
                                                    'Gene', 'Transcript', 'Exon_Length'], index=False)

    #to save Ortho_Table
    Ortho_Table.to_csv("Ortho_Table.csv", columns= ['Exon', 'Exon2', 'Nuc_Length_of_Overlap', 'Prot_Length_of_Overlap',
                                                    'Nuc_PID', 'Prot_PID', 'Relationship'], index=False)






    def create_blast_table(species1name, species2name, species1version, species2version):
        Ortho_Table = pd.read_table('/Users/rhythmicstar/blast/db//nuctable.html',
                                    names=['Exon', 'Exon2', 'Nuc_Length_of_Overlap', 'Nuc_PID'])

        Prot_Table = pd.read_table('/Users/rhythmicstar/blast/db//protable.html',
                                    names=['Exon', 'Exon2', 'Prot_Length_of_Overlap', 'Prot_PID'])

        #add the other columns to the Ortho_Table
        Ortholog_Columns = pd.DataFrame(columns=['Prot_PID'])
        #Overlap_Column = pd.DataFrame(columns=['Prot_Length_of_Overlap'])
        Ortho_Table.insert(3, 'Prot_Length_of_Overlap', 'Below Threshold of ...')

        #append new dataframe to Ortho_Table
        Ortho_Table = pd.concat([Ortho_Table, Ortholog_Columns], axis=1)
        Ortho_Table = Ortho_Table.replace(np.nan,'Below Threshold of ...', regex=True)

        BLAST_Table = []

        #to drop a row that has the same two exons
        rows = int(Ortho_Table.size/6)
        #loop through each row
        for row in range (0, rows):
            if(Ortho_Table.iat[row, 0] == Ortho_Table.iat[row, 1]):
                Ortho_Table.ix[row, 0] = np.nan
        #remove each duplicate
        Ortho_Table.dropna(axis=0,inplace=True)
        Ortho_Table = Ortho_Table.reset_index(drop=True)

        protrow = int(Prot_Table.size/4)
        #loop through each row
        for row in range (0, protrow):
            if(Prot_Table.iat[row, 0] == Prot_Table.iat[row, 1]):
                Prot_Table.ix[row, 0] = np.nan

        #remove each duplicate
        Prot_Table.dropna(axis=0,inplace=True)

        Prot_Table = Prot_Table.reset_index(drop=True)

        #to remove rows in nuc table that have the same two exons in the opposite order
        seen = set([])
        rows_to_use = []

        for i, row in Ortho_Table.iterrows():
            exon1 = row['Exon']
            exon2 = row['Exon2']
            pair = tuple(sorted([exon1, exon2]))

            if pair not in seen:
                seen.update({pair})
                rows_to_use.append(i)

        ortho_table_no_duplicate_comparisons = Ortho_Table.loc[rows_to_use]
        ortho_table_no_duplicate_comparisons = ortho_table_no_duplicate_comparisons.reset_index(drop=True)
        Ortho_Table = ortho_table_no_duplicate_comparisons

        #to remove rows in protein table that have the same two exons in the opposite order
        seen = set([])
        rows_to_use = []

        for i, row in Prot_Table.iterrows():
            exon1 = row['Exon']
            exon2 = row['Exon2']
            pair = tuple(sorted([exon1, exon2]))

            if pair not in seen:
                seen.update({pair})
                rows_to_use.append(i)

        prot_table_no_duplicate_comparisons = Prot_Table.loc[rows_to_use]
        prot_table_no_duplicate_comparisons = prot_table_no_duplicate_comparisons.reset_index(drop=True)
        Prot_Table = prot_table_no_duplicate_comparisons

        rows = int(Ortho_Table.size/6)
        protrow = int(Prot_Table.size/4)

        for rowN in range (0, rows):
            for rowP in range (0,protrow):
                exon1N = Ortho_Table.iat[rowN, 0]
                exon1P = Prot_Table.iat[rowP, 0]
                exon2N = Ortho_Table.iat[rowN, 1]
                exon2P = Prot_Table.iat[rowP, 1]
                if ((exon1N == exon1P) & (exon2N == exon2P)):
                    #put protein pid in nuc table
                    Ortho_Table.iat[rowN, 5] = Prot_Table.iat[rowP, 3]
                    Ortho_Table.iat[rowN, 3] = Prot_Table.iat[rowP, 2]
                    #put na into protein table
                    Prot_Table.ix[rowP, 0] = np.nan
                    break

        #remove each duplicate
        Prot_Table.dropna(axis=0,inplace=True)
        Prot_Table = Prot_Table.reset_index(drop=True)

        # #in Prot_Table, add col with 'Nuc_PID' between 'Length_of_Overlap' and 'Prot_PID'
        Prot_Table.insert(2, 'Nuc_Length_of_Overlap', 'Below Threshold of ...')
        Prot_Table.insert(4, 'Nuc_PID', 'Below Threshold of ...')

        #add protein table at end
        Ortho_Table = pd.concat([Ortho_Table, Prot_Table])
        Ortho_Table = Ortho_Table.reset_index(drop=True)

        #create new dataframe
        Relationship_Column = pd.DataFrame(columns=['Relationship'])

        #append new dataframe to Ortho_Table
        Ortho_Table = pd.concat([Ortho_Table, Relationship_Column], axis=1)

        #create BLAST_Table dataframe with exon ids
        Ref_BLAST_Table = copy.deepcopy(Ortho_Table)
        rows = int(Ref_BLAST_Table.size/7)

        for row in range (0, rows):
            Ref_BLAST_Table.ix[(row + rows), 0] = Ref_BLAST_Table.ix[row, 1]

        cols = [1,2,3,4,5,6]
        Ref_BLAST_Table.drop(Ref_BLAST_Table.columns[cols],axis=1,inplace=True)
        Ref_BLAST_Table = Ref_BLAST_Table.drop_duplicates('Exon')
        Ref_BLAST_Table = Ref_BLAST_Table.reset_index(drop=True)

        #in form to do
        #get data to fill in table
        rows = int(Ref_BLAST_Table.size)

        for exonCode in species1DB.features_of_type('CDS'):
            Species1Code = str(exonCode['exon_id'])
            Species1Code = separate(Species1Code)
            Species1Code = Species1Code[0:7]
            break

        #loop through each exon and get its gene and length
        for row in range (0, rows):
            exon_ID = Ref_BLAST_Table.ix[row, 'Exon']
            exon_ID2 = separate(exon_ID)

            #determine which file to look in
            if (exon_ID2[0:7] == Species1Code):
                for exon in species1DB.features_of_type('CDS'):
                    GFFUtilsExonId = str(exon['exon_id']) #gffutils id
                    GFFUtilsExonId = separate(GFFUtilsExonId)
                    if (GFFUtilsExonId == exon_ID2):
                        exon_geneid = str(exon.attributes['gene_id'])
                        exon_geneid = separate(exon_geneid)

                        exon_transcriptid = str(exon.attributes['transcript_id'])
                        exon_transcriptid = separate(exon_transcriptid)

                        New_BLAST_Row = [exon_ID, ("{}{}{}{}").format(species1name,'(', species1version, ')'),
                                         ("{}:{}-{}:{}:{}".format(exon.chrom, exon.start, exon.stop,
                                                                  exon.strand, exon.frame)), exon_geneid,
                                         exon_transcriptid, exon.stop - exon.start + 1]

                        BLAST_Table.append(New_BLAST_Row)

            else:
                #if exon is mouse
                for exon in species2DB.features_of_type('CDS'):
                    GFFUtilsExonId = str(exon['exon_id']) #gffutils id
                    GFFUtilsExonId = separate(GFFUtilsExonId)
                    if (GFFUtilsExonId == exon_ID2):
                        exon_geneid = str(exon.attributes['gene_id'])
                        exon_geneid = separate(exon_geneid)

                        exon_transcriptid = str(exon.attributes['transcript_id'])
                        exon_transcriptid = separate(exon_transcriptid)

                        New_BLAST_Row = [exon_ID, ("{}{}{}{}").format(species2name,'(', species2version, ')'),
                                         ("{}:{}-{}:{}:{}".format(exon.chrom, exon.start, exon.stop,
                                                                  exon.strand, exon.frame)), exon_geneid,
                                         exon_transcriptid, exon.stop - exon.start + 1]

                        BLAST_Table.append(New_BLAST_Row)

        BLAST_Table = pd.DataFrame(BLAST_Table)
        BLAST_Table.columns = ['Exon', 'Species(Version)', 'Chrom:Start-Stop:Strand:Offset', 'Gene', 'Transcript', 'Exon_Length']

        #fill in row for paralogous or orthologous
        rows = int(Ortho_Table.size/7)
        blastRows = int(BLAST_Table.size/6)
        for row in range (0, rows):
        #paralogous
            if (((Ortho_Table.iat[row, 0])[0:7]) != ((Ortho_Table.iat[row, 1])[0:7])):
                Ortho_Table.ix[row, 'Relationship'] = 'Orthologous'

        #orthologous
            else:
                exonOne = Ortho_Table.iat[row, 0]
                exonTwo = Ortho_Table.iat[row, 1]
                locationOne = 0
                locationTwo = 1
                start1 = '0'
                start2 = '1'
                end1 = '0'
                end2 = '1'
                for innerrow in range (0, blastRows):
                    if (exonOne == BLAST_Table.iat[innerrow, 0]):
                        locationOne = BLAST_Table.iat[innerrow, 2]
                        start1 = splitstart(locationOne)
                        end1 = splitend(locationOne)
                    if (exonTwo == BLAST_Table.iat[innerrow, 0]):
                        locationTwo = BLAST_Table.iat[innerrow, 2]
                        start2 = splitstart(locationTwo)
                        end2 = splitend(locationTwo)
                    if (((start1 >= start2) and (end2 >= end1)) or
                        ((start2 >= start1) and (end1 >= end2)) or
                        ((start2 >= start1) and (end2 >= end1) and
                         (start2 <= end1)) or
                        ((start1 >= start2) and (end1 >= end2) and
                         (start1 <= end2))):
                        Ortho_Table.iat[row, 6] = 'Overlapping Genomic Loci'
                        if (locationOne == locationTwo):
                            Ortho_Table.iat[row, 6] = 'Identical Genomic Loci'
                    else:
                        Ortho_Table.iat[row, 6] = 'Paralogous'

        #to save BLAST_Table
        BLAST_Table.to_csv("BLAST_Table.csv", columns= ['Exon', 'Species(Version)', 'Chrom:Start-Stop:Strand:Offset',
                                                        'Gene', 'Transcript', 'Exon_Length'], index=False)

        #to save Ortho_Table
        Ortho_Table.to_csv("Ortho_Table.csv", columns= ['Exon', 'Exon2', 'Nuc_Length_of_Overlap', 'Prot_Length_of_Overlap',
                                                        'Nuc_PID', 'Prot_PID', 'Relationship'], index=False)