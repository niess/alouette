from ._core import ffi, lib
import numpy

__all__ = ('decay', 'initialise', 'Products', 'undecay')


class Products:
    '''Proxy for decay products

    This class is a proxy for accessing the products of a decay. It exposes
    the C struct data as read-only numpy arrays.
    '''

    def __init__(self, c_struct):
        self._c = c_struct

    def __repr__(self):
        return f'{type(self).__name__}(size={self.size}, pid={self.pid}, '     \
            f'polarimetric={self.polarimetric}, P={self.P})'

    @property
    def P(self):
        '''Four momenta of the N decay products as an N x 4 numpy array
        '''
        try:
            return self._P
        except AttributeError:
            n = self._c.size
            self._P = numpy.frombuffer(
                ffi.buffer(self._c.P)[:])[:4 * n].reshape((n, 4))
            return self._P

    @property
    def pid(self):
        '''Particle IDs of the N decay products, using PDG numbering scheme
        '''
        try:
            return self._pid
        except AttributeError:
            self._pid = numpy.frombuffer(ffi.buffer(self._c.pid)[:],
                dtype='i4')[:self._c.size]
            return self._pid

    @property
    def polarimetric(self):
        '''Polarimetric vector of the decay
        '''
        try:
            return self._polarimetric
        except AttributeError:
            self._polarimetric = numpy.frombuffer(
                ffi.buffer(self._c.polarimetric)[:])
            return self._polarimetric

    @property
    def size(self):
        '''Number, N, of decay products
        '''
        return int(self._c.size)

    @property
    def weight(self):
        '''Monte Carlo weight of the decay
        '''
        return float(self._c.weight)


# Mapping between ALOUETTE exceptions codes and Python classes
_index_to_exception = (
        Exception,          # ALOUETTE_RETURN_SUCCESS
        ValueError,         # ALOUETTE_RETURN_VALUE_ERROR
        FloatingPointError, # ALOUETTE_RETURN_FLOATING_ERROR
        RuntimeError        # ALOUETTE_RETURN_TAUOLA_ERROR
)


def _call(f, *args):
    '''Library call with error check
    '''
    rc = f(*args)
    if rc != 0:
        exception = _index_to_exception[int(rc)]
        message = ffi.string(lib.alouette_message()).decode()
        raise exception(message)


def initialise(xk0dec=None):
    '''Initialise TAUOLA using custom settings
    '''

    if xk0dec is None:
        xk0dec = ffi.NULL
    else:
        xk0dec = ffi.new('double [1]', (xk0dec,))

    _call(lib.alouette_initialise, xk0dec)


def decay(mode=None, pid=None, momentum=None, polarisation=None):
    '''Decay a tau particle using tauola
    '''

    if mode is None:
        mode = 0

    if pid is None:
        pid = 15

    if momentum is None:
        momentum = ffi.new('double [3]', (0,0,0))
    else:
        momentum = ffi.new('double [3]', momentum)

    if polarisation is None:
        polarisation = ffi.NULL
    else:
        polarisation = ffi.new('double [3]', polarisation)

    products = ffi.new('struct alouette_products *')
    _call(lib.alouette_decay, mode, pid, momentum, polarisation, products)

    return Products(products)


'''User supplied polarisation callback, for backward decays
'''
_polarisation = None


@ffi.def_extern()
def _polarisation_callback(pid, momentum, polarisation):
    '''Callback wrapper, for setting the polarisation in backward decays
    '''
    m = numpy.frombuffer(ffi.buffer(momentum, 24)[:])
    polarisation[0:3] = _polarisation(int(pid), m)


def undecay(mode=None, pid=None, momentum=None, polarisation=None, bias=None):
    '''Backward Monte Carlo decay to a tau particle using tauola
    '''

    if mode is None:
        mode = 0

    if pid is None:
        pid = 16

    if momentum is None:
        momentum = ffi.new('double [3]', (0,0,0))
    else:
        momentum = ffi.new('double [3]', momentum)

    if polarisation is None:
        polar_cb = ffi.NULL
    else:
        global _polarisation
        _polarisation = polarisation
        polar_cb = lib._polarisation_callback

    if bias is None:
        bias = 0

    products = ffi.new('struct alouette_products *')
    _call(lib.alouette_undecay, mode, pid, momentum, polar_cb, bias, products)

    return Products(products)
