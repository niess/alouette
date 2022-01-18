import pytest

import numpy
import alouette
from alouette.tauola import ipcht, wtmax


def test_tauola():
    '''Test the TAUOLA API'''

    # Check initial weights
    z15 = numpy.zeros(15)
    assert wtmax.dadmaa == 0
    assert wtmax.dadmel == 0
    assert wtmax.dadmks == 0
    assert wtmax.dadmmu == 0
    assert wtmax.dadmro == 0
    assert numpy.all(wtmax.dadnew == z15)

    # Initialise and check new weights
    alouette.initialise()

    assert wtmax.dadmaa > 0
    assert wtmax.dadmel > 0
    assert wtmax.dadmks > 0
    assert wtmax.dadmmu > 0
    assert wtmax.dadmro > 0
    assert numpy.all(wtmax.dadnew > z15)

    # Check mode flag
    assert ipcht.iver == 1
    ipcht.iver = 0
    assert ipcht.iver == 0
    ipcht.iver = 1
    assert ipcht.iver == 1

    with pytest.raises(TypeError) as exc_info:
        ipcht.iver = 'toto'
    assert 'bad iver' in exc_info.exconly()

    with pytest.raises(ValueError) as exc_info:
        ipcht.iver = 2
    assert 'bad iver' in exc_info.exconly()
