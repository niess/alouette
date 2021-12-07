import alouette


def test_random():
    '''Test the random API'''

    # Check the default PRNG seed
    seed = alouette.random.seed
    assert isinstance(seed, int)
    assert seed != 0  # Unlikely, yet it could happen

    # Get the first event for latter comparison
    products0 = alouette.decay()

    # Check the seed setter
    alouette.random.set(0)
    assert alouette.random.seed == 0
    assert alouette.decay() != products0

    alouette.random.set(None)
    assert alouette.random.seed != 0
    assert alouette.decay() != products0

    # Check the reset of the random stream
    alouette.random.set(seed)
    assert alouette.decay() == products0

    # Check the PRNG getter
    u = alouette.random()
    assert isinstance(u, float)
    assert 0 <= u <= 1
    assert u != alouette.random()
