from ._core import ffi, lib

__all__ = ('Random', 'random')


class Random:
    '''Proxy to ALOUETTE's Pseudo Random Numbers Generator (PRNG)
    '''

    @property
    def seed(self):
        return lib.alouette_random_seed()

    def set(self, seed=None):
        '''Reset the state of the PRNG
        '''
        if seed is None:
            seed = ffi.NULL
        else:
            seed = ffi.new('unsigned long *', seed)
        lib.alouette_random_set(seed)

    def __call__(self):
        return float(lib.alouette_random())


random = Random()
