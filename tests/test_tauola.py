import numpy
import alouette
from alouette.tauola import wtmax


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
    alouette.core.initialise()

    assert wtmax.dadmaa > 0
    assert wtmax.dadmel > 0
    assert wtmax.dadmks > 0
    assert wtmax.dadmmu > 0
    assert wtmax.dadmro > 0
    assert numpy.all(wtmax.dadnew > z15)
