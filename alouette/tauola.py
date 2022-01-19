from numbers import Number
import numpy

from ._core import ffi, lib

__all__ = ('ipcht', 'parmas', 'wtmax',)


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


class parmas:
    '''TAUOLA mass parameters'''

    @property
    def amtau(self):
        '''Tau mass, in GeV'''
        return float(lib.tauola_parmas.amtau)

    @property
    def amnuta(self):
        '''Tau neutrino mass, in GeV'''
        return float(lib.tauola_parmas.amnuta)

    @property
    def amel(self):
        '''Electron mass, in GeV'''
        return float(lib.tauola_parmas.amel)

    @property
    def amnue(self):
        '''Electron neutrino mass, in GeV'''
        return float(lib.tauola_parmas.amnue)

    @property
    def ammu(self):
        '''Muon mass, in GeV'''
        return float(lib.tauola_parmas.ammu)

    @property
    def amnumu(self):
        '''Muon neutrino mass, in GeV'''
        return float(lib.tauola_parmas.amnumu)

    @property
    def ampiz(self):
        '''Neutral pion mass, in GeV'''
        return float(lib.tauola_parmas.ampiz)

    @property
    def ampi(self):
        '''Charged pion mass, in GeV'''
        return float(lib.tauola_parmas.ampi)

    @property
    def amro(self):
        '''Rho mass, in GeV'''
        return float(lib.tauola_parmas.amro)

    @property
    def gamro(self):
        '''Rho decay width, in GeV'''
        return float(lib.tauola_parmas.gamro)

    @property
    def ama1(self):
        '''A1 mass, in GeV'''
        return float(lib.tauola_parmas.ama1)

    @property
    def gama1(self):
        '''A1 decay width, in GeV'''
        return float(lib.tauola_parmas.gama1)

    @property
    def amk(self):
        '''Charged kaon mass, in GeV'''
        return float(lib.tauola_parmas.amk)

    @property
    def amkz(self):
        '''Neutral kaon mass, in GeV'''
        return float(lib.tauola_parmas.amkz)

    @property
    def amkst(self):
        '''K*(892) mass, in GeV'''
        return float(lib.tauola_parmas.amkst)

    @property
    def gamkst(self):
        '''K*(892) decay width, in GeV'''
        return float(lib.tauola_parmas.gamkst)

# Override as class singleton
parmas = parmas()
