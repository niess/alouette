import alouette
import numpy
import pytest


def almost_equal(a, b, epsilon=None):
    '''Compare floats up to epsilon
    '''
    epsilon = epsilon or 10 * numpy.finfo(numpy.float32).eps
    assert(pytest.approx(a, abs=epsilon) == pytest.approx(b, abs=epsilon))


def test_decay():
    '''Test the decay function
    '''

    # Disable radiative corrections
    alouette.core.initialise(xk0dec=0)

    # Check the pid argument
    product = alouette.decay()
    assert(product.pid[0] == 16)

    product = alouette.decay(pid=15)
    assert(product.pid[0] == 16)

    product = alouette.decay(pid=-15)
    assert(product.pid[0] == -16)

    # Check the energy-momentum-conservation
    mtau = 1.777
    product = alouette.decay(mode=1)
    P = numpy.sum(product.P, axis=0)
    for i in range(3):
        almost_equal(P[i], 0)
    almost_equal(P[3], mtau, epsilon=1E-03)

    momentum = numpy.array((0, 0, 1))
    energy = numpy.sqrt(sum(momentum**2) + mtau**2)
    product = alouette.decay(mode=1, momentum=momentum)
    P = numpy.sum(product.P, axis=0)
    for i in range(3):
        almost_equal(P[i], momentum[i])
    almost_equal(P[3], energy, epsilon=1E-03)

    # Check the polarisation argument
    polarisation = -momentum / numpy.linalg.norm(momentum)
    alouette.decay(momentum=momentum, polarisation=polarisation)
    alouette.decay(momentum=(0, 0, 1), polarisation=(0, 0, -1))
