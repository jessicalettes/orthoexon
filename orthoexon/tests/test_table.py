__author__ = 'rhythmicstar'

import os

import pandas as pd
import pandas.util.testing as pdt
import pytest

class TestOrthologyTable(object):
    
    @pytest.fixture
    def nucleotide_filename(self, table_folder):
        return os.path.join(table_folder, 'nuctable.html')

    @pytest.fixture
    def protein_filename(self, table_folder):
        return os.path.join(table_folder, 'protable.html')

    def test___init(self, nucleotide_filename, protein_filename):
        from orthoexon.table import OrthologyTable

        true_nucleotide = pd.read_csv(
            nucleotide_filename.replace('.html', '.csv'))

        t = OrthologyTable(nucleotide_filename, protein_filename, None, None,
                           None, None)

        pdt.assert_frame_equal(t.nucleotide, true_nucleotide)