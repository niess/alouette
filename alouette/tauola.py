import numpy

from ._core import ffi, lib

__all__ = ('wtmax',)


class wtmax:
    '''TAUOLA maximum weights, computed during initialisation'''

    @property
    def dadmaa(self):
        '''Maximum weight for the DADMAA procedure
        '''
        return float(lib.tauola_weight_dadmaa.wtmax)

    @property
    def dadmel(self):
        '''Maximum weight for the DADMEL procedure
        '''
        return float(lib.tauola_weight_dadmel.wtmax)

    @property
    def dadmmu(self):
        '''Maximum weight for the DADMMU procedure
        '''
        return float(lib.tauola_weight_dadmmu.wtmax)

    @property
    def dadmks(self):
        '''Maximum weight for the DADMKS procedure
        '''
        return float(lib.tauola_weight_dadmks.wtmax)

    @property
    def dadmro(self):
        '''Maximum weight for the DADMRO procedure
        '''
        return float(lib.tauola_weight_dadmro.wtmax)

    @property
    def dadnew(self):
        '''Maximum weight for the DADNEW procedure
        '''
        return numpy.frombuffer(
            ffi.buffer(lib.tauola_weight_dadnew.wtmax)[:], dtype='f4'
        )

# Override as class singleton
wtmax = wtmax()
