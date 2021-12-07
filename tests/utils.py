'''Utilities for unit tests'''
import numpy
import pytest


def almost_equal(a, b, epsilon=None):
    '''Compare floats up to epsilon'''
    epsilon = epsilon or 10 * numpy.finfo(numpy.float32).eps
    assert pytest.approx(a, abs=epsilon) == pytest.approx(b, abs=epsilon)
