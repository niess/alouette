import pytest

import numpy
import alouette
from alouette.tauola import ipcht, parmas, wtmax


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

    # Check mass getters
    assert parmas.amtau > 0
    assert parmas.amnuta > 0
    assert parmas.amel > 0
    assert parmas.amnue == 0
    assert parmas.ammu > 0
    assert parmas.amnumu == 0
    assert parmas.ampi > 0
    assert parmas.ampiz > 0
    assert parmas.amro > 0
    assert parmas.gamro > 0
    assert parmas.ama1 > 0
    assert parmas.gama1 > 0
    assert parmas.amk > 0
    assert parmas.amkz > 0
    assert parmas.amkst > 0
    assert parmas.gamkst > 0

    with pytest.raises(AttributeError):
        parmas.amtau = 0
