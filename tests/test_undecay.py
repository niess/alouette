import numpy
import pytest
import alouette
from .utils import almost_equal


def test_undecay():
    '''Test the undecay wrapper'''

    # Disable radiative corrections
    alouette.initialise(xk0dec=0)

    # Check the default settings
    assert alouette.undecay.mother == 0
    assert alouette.undecay.bias == 1

    # Check the pid and mother arguments
    momentum = numpy.array((0, 0, 1))
    products = alouette.undecay(momentum=momentum)
    assert products.pid[0] == 15

    products = alouette.undecay(pid=16, momentum=momentum)
    assert products.pid[0] == 15

    products = alouette.undecay(pid=-16, momentum=momentum)
    assert products.pid[0] == -15
    alouette.undecay.mother = 15
    products = alouette.undecay(mode=1, pid=11, momentum=momentum)
    assert products.pid[0] == 15
    assert products.pid[1] == 16
    assert products.pid[2] == -12

    alouette.undecay.mother = -15
    products = alouette.undecay(mode=1, pid=12, momentum=momentum)
    assert products.pid[0] == -15
    assert products.pid[1] == -16
    assert products.pid[2] == -11

    alouette.undecay.mother = 15
    with pytest.raises(ValueError) as exc_info:
        alouette.undecay(mode=1, pid=12, momentum=momentum)
    assert (
        'inconsistent values for mother (15), daughter (12) and mode (1)'
        in exc_info.exconly()
    )

    # Check the null momentum case
    alouette.undecay.mother = 0
    with pytest.raises(ValueError) as exc_info:
        alouette.undecay()
    assert 'bad daughter momentum (0)' in exc_info.exconly()

    # Check the energy-momentum-conservation
    mnu = 0.01
    energy = numpy.sqrt(sum(momentum ** 2) + mnu ** 2)
    products = alouette.undecay(mode=1, momentum=momentum)
    P = numpy.sum(products.P[1:], axis=0)
    P[:3] += momentum
    P[3] += energy
    for i in range(4):
        almost_equal(P[i], products.P[0, i])

    # Check the polarisation argument
    def polarisation_cb(pid, momentum):
        return -momentum / numpy.linalg.norm(momentum)

    alouette.undecay(momentum=momentum, polarisation=polarisation_cb)
    alouette.undecay(momentum=(0, 0, 1), polarisation=polarisation_cb)

    # Check the mother parameter
    alouette.undecay.mother = 15
    assert alouette.undecay.mother == 15
    alouette.undecay.mother = -15
    assert alouette.undecay.mother == -15
    with pytest.raises(TypeError):
        alouette.undecay.mother = 'tau'
    with pytest.raises(ValueError):
        alouette.undecay.mother = 11
    alouette.undecay.mother = 0
    assert alouette.undecay.mother == 0

    # Check the bias parameter
    alouette.undecay.bias = 0
    assert alouette.undecay.bias == 0
    alouette.undecay.bias = -1
    assert alouette.undecay.bias == -1
    alouette.undecay(momentum=momentum)
    with pytest.raises(TypeError):
        alouette.undecay.bias = '1'
    with pytest.raises(ValueError):
        alouette.undecay.bias = 2
    alouette.undecay.bias = 1
    assert alouette.undecay.bias == 1

    # Check the scheme parameter
    assert alouette.undecay.scheme == 'cartesian'
    alouette.undecay.scheme = 'spherical'
    assert alouette.undecay.scheme == 'spherical'
    alouette.undecay.scheme = 'energy'
    assert alouette.undecay.scheme == 'energy'
    alouette.undecay.scheme = 'cartesian'
    assert alouette.undecay.scheme == 'cartesian'
    with pytest.raises(TypeError):
        alouette.undecay.scheme = 1
    with pytest.raises(ValueError):
        alouette.undecay.scheme = 'nothing'

