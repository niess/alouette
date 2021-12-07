import alouette
import numpy
import pytest


def almost_equal(a, b, epsilon=None):
    '''Compare floats up to epsilon'''
    epsilon = epsilon or 10 * numpy.finfo(numpy.float32).eps
    assert pytest.approx(a, abs=epsilon) == pytest.approx(b, abs=epsilon)


def test_undecay():
    '''Test the undecay function'''

    # Disable radiative corrections
    alouette.core.initialise(xk0dec=0)

    # Check the pid and mother arguments
    momentum = numpy.array((0, 0, 1))
    products = alouette.undecay(momentum=momentum)
    assert products.pid[0] == 15

    products = alouette.undecay(pid=16, momentum=momentum)
    assert products.pid[0] == 15

    products = alouette.undecay(pid=-16, momentum=momentum)
    assert products.pid[0] == -15

    products = alouette.undecay(mode=1, pid=11, mother=15, momentum=momentum)
    assert products.pid[0] == 15
    assert products.pid[1] == 16
    assert products.pid[2] == -12

    products = alouette.undecay(mode=1, pid=12, mother=-15, momentum=momentum)
    assert products.pid[0] == -15
    assert products.pid[1] == -16
    assert products.pid[2] == -11

    with pytest.raises(ValueError) as exc_info:
        alouette.undecay(mode=1, pid=12, mother=15, momentum=momentum)
    assert (
        'inconsistent values for mother (15), daugther (12) and mode (1)'
        in exc_info.exconly()
    )

    # Check the null momentum case
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

    # Check the bias argument
    alouette.undecay(momentum=momentum, bias=0)
    alouette.undecay(momentum=momentum, bias=-1)
    with pytest.raises(TypeError):
        alouette.undecay(momentum=momentum, bias='1')
