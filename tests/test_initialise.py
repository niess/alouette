import pytest
import alouette


def test_initialise():
    '''Test explicit initialisation'''

    # Check that the PRNG is preserved
    alouette.random.set(12345)
    alouette.random()
    u1 = alouette.random()
    alouette.random.set(12345)
    alouette.random()

    alouette.core.initialise()
    assert alouette.random.seed == 12345
    assert alouette.random() == u1

    # Check re-initialisation error
    with pytest.raises(RuntimeError) as exc_info:
        alouette.core.initialise(xk0dec=0)
    assert 'TAUOLA already initialised' in exc_info.exconly()
