import alouette
import numpy
import pytest


def test_products():
    '''Test the decay products wrapper'''

    # Check the accessors
    alouette.core.initialise(xk0dec=0)  # Disable radiative corrections
    product = alouette.decay(mode=1)

    assert isinstance(product.pid, numpy.ndarray)
    assert product.pid.shape == (3,)
    assert product.pid.dtype == numpy.int32

    assert isinstance(product.P, numpy.ndarray)
    assert product.P.shape == (3, 4)
    assert product.P.dtype == numpy.float64

    assert isinstance(product.polarimeter, numpy.ndarray)
    assert product.polarimeter.shape == (3,)
    assert product.polarimeter.dtype == numpy.float64

    assert isinstance(product.size, int)
    assert product.size == 3

    assert isinstance(product.weight, float)
    assert product.weight == 1

    # Check the caching
    assert product.pid is product.pid
    assert product.P is product.P
    assert product.polarimeter is product.polarimeter

    # Check the representation (roughly)
    rep = str(product)
    assert rep.startswith('Products')
    assert 'pid=[ 16 ' in rep

    # Check the comparison operator
    assert product == product
    assert product != alouette.decay(mode=1)
    assert product != 1

    # Check that data are read-only
    with pytest.raises(AttributeError) as exc_info:
        product.weight = 0
    assert "can't set attribute" in exc_info.exconly()

    with pytest.raises(ValueError) as exc_info:
        product.pid[0] = 0
    assert "assignment destination is read-only" in exc_info.exconly()

    with pytest.raises(ValueError) as exc_info:
        product.P[0][0] = 0
    assert "assignment destination is read-only" in exc_info.exconly()

    with pytest.raises(ValueError) as exc_info:
        product.polarimeter[0] = 0
    assert "assignment destination is read-only" in exc_info.exconly()
