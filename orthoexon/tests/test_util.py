#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
test_orthoexon
----------------------------------

Tests for `orthoexon` module.
"""

import pytest

@pytest.fixture
def gene_id_with_quotes():
    return "'ENSMUSG00000064842.1'"

@pytest.fixture
def gene_id():
    return "ENSMUSG00000064842.1"

def test_separate_with_quotes(gene_id_with_quotes):
    from orthoexon.util import separate

    test = separate(gene_id_with_quotes)
    true = "ENSMUSG00000064842"

    assert test == true



def test_separate(gene_id):
    from orthoexon.util import separate

    test = separate(gene_id)
    true = "ENSMUSG00000064842"

    assert test == true


@pytest.fixture
def location():
    return "chr20:10256140-10256211:+:0"

def test_splitstart(location):
    from orthoexon.util import splitstart

    test = splitstart(location)
    true = '10256140'

    assert test == true

