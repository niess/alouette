from ._core import ffi, lib

__all__ = ('random',)


class random:
    '''Alouette's Pseudo Random Numbers Generator (PRNG)'''

    @property
    def seed(self):
        '''The PRNG current seed
        '''
        return lib.alouette_random_seed()

    @staticmethod
    def set(seed=None):
        '''Reset the state of the PRNG'''
        if seed is None:
            seed = ffi.NULL
        else:
            seed = ffi.new('unsigned long *', seed)
        lib.alouette_random_set(seed)

    def __call__(self):
        return float(lib.alouette_random())


# Override as class singleton
random = random()
