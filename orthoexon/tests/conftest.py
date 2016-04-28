__author__ = 'rhythmicstar'

import os

import pytest

@pytest.fixture
def data_folder():
    return os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data')


@pytest.fixture(params=[True, False])
def maybe(request):
    return request.param