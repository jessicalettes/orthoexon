#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
test_orthoexon
----------------------------------

Tests for `orthoexon` module.
"""
import os
import pytest


@pytest.fixture
def exon_id_with_quotes():
    return "'ENSE00001229068.1'"

@pytest.fixture
def exon_id():
    return "ENSE00001229068.1"

def test_separate_with_quotes(exon_id_with_quotes):
    from orthoexon.util import separate

    test = separate(exon_id_with_quotes)
    true = "ENSE00001229068"

    assert test == true



def test_separate(exon_id):
    from orthoexon.util import separate

    test = separate(exon_id)
    true = "ENSE00001229068"

    assert test == true


@pytest.fixture
def location():
    return "chr20:10256140-10256211:+:0"

def test_splitstart(location):
    from orthoexon.util import splitstart

    test = splitstart(location)
    true = '10256140'

    assert test == true


def test_splitend(location):
    from orthoexon.util import splitend

    test = splitend(location)
    true = '10256211'

    assert test == true


@pytest.fixture
def human_gtf_filename(table_folder):
    return os.path.join(table_folder, 'humanrbfox2andfmr1andsnap25.gtf')

@pytest.fixture
def human_gtf_database(table_folder):
    return os.path.join(table_folder, 'humanrbfox2andfmr1andsnap25.gtf.db')

@pytest.fixture
def human_fasta(table_folder):
    return os.path.join(table_folder, 'GRCm38.p3.genome.fa')

def test_translate(exon_id, human_fasta, human_gtf_database):
    from orthoexon.util import translate
    from orthoexon.util import separate
    for index, species1gene in enumerate(human_gtf_database.features_of_type('gene')):
        species1gffutilsgeneid = str(species1gene['gene_id'])
        species1geneid = separate(species1gffutilsgeneid)

        for exon in human_gtf_database.children(species1geneid,
                                                featuretype='CDS',
                                                order_by='start'):
            if exon_id == exon:
                test = translate(exon, human_fasta)
            break
        break
    true = 'MAEDADMRNELEEMQRRADQLADE'
    assert test == true



# def test_getsequence(exon, human_gtf_database):
#     from orthoexon.util import getsequence
#
#     test = getsequence(exon, human_gtf_database)
#     true = 'ATGGCCGAAGACGCAGACATGCGCAATGAGCTGGAGGAGATGCAGCGAAGGGCTGACCAGTT' \
#            'GGCTGATGAG'
#
#     assert test == true

# def test_make_sequence_array(finalsequencedf):
#     from orthoexon.util import make_sequence_array
#
#     test = make_sequence_array(finalsequencedf)
#     true = ......
#
#     assert test == true