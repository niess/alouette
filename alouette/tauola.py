from numbers import Number
import numpy

from ._core import ffi, lib

__all__ = ('ipcht', 'wtmax',)


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


class ipcht:
    '''TAUOLA new currents flags'''

    def __init__(self):
        self.wtmax = wtmax

    @property
    def iver(self):
        '''Parametrisation version for decays to 2 and 3 pions.

           This flag allows to switch forth and back between (legacy) CLEO
           parameterisation and (new currents) Belle + RChL for decays to 2 and
           3 pions.  By default, TAUOLOA is set with iver=1, i.e. Belle + RChL.
        '''
        return int(lib.tauola_ipcht.iver)

    @iver.setter
    def iver(self, v):
        if isinstance(v, Number):
            if v in (0, 1):
                lib.tauola_ipcht.iver = v
            else:
                raise ValueError('bad iver')
        else:
            raise TypeError('bad iver')

# Override as class singleton
ipcht = ipcht()
