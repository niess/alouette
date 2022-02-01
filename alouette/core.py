import numpy
from numbers import Number

from ._core import ffi, lib

__all__ = ('decay', 'initialise', 'Products', 'undecay')


class Products:
    '''Proxy for decay products

    This class is a proxy for accessing the products of a decay. It exposes
    the C struct data as read-only numpy arrays.
    '''

    def __init__(self, c_struct):
        self._c = c_struct
        self._P = None
        self._pid = None
        self._polarimeter = None

    def __repr__(self):
        return (
            f'{type(self).__name__}(size={self.size}, pid={self.pid}, '
            f'polarimeter={self.polarimeter}, P={self.P})'
        )

    def __eq__(self, other):
        try:
            return ffi.buffer(self._c) == ffi.buffer(other._c)
        except AttributeError:
            return False

    @property
    def P(self):
        '''Four momenta of the N decay products as an N x 4 numpy array'''
        if self._P is None:
            n = self._c.size
            self._P = numpy.frombuffer(ffi.buffer(self._c.P)[:])[
                : 4 * n
            ].reshape((n, 4))
        return self._P

    @property
    def pid(self):
        '''Particle IDs of the N decay products, using PDG numbering scheme'''
        if self._pid is None:
            self._pid = numpy.frombuffer(
                ffi.buffer(self._c.pid)[:], dtype='i4'
            )[: self._c.size]
        return self._pid

    @property
    def polarimeter(self):
        '''Polarimeter vector of the decay'''
        if self._polarimeter is None:
            self._polarimeter = numpy.frombuffer(
                ffi.buffer(self._c.polarimeter)[:]
            )
        return self._polarimeter

    @property
    def size(self):
        '''Number, N, of decay products'''
        return int(self._c.size)

    @property
    def weight(self):
        '''Monte Carlo weight of the decay'''
        return float(self._c.weight)


# Mapping between Alouette exceptions codes and Python classes
_index_to_exception = (
    Exception,  # ALOUETTE_RETURN_SUCCESS
    ValueError,  # ALOUETTE_RETURN_VALUE_ERROR
    RuntimeError,  # ALOUETTE_RETURN_TAUOLA_ERROR
)


def _call(f, *args):
    '''Library call with error check'''
    rc = f(*args)
    if rc != 0:
        exception = _index_to_exception[int(rc)]
        message = ffi.string(lib.alouette_message()).decode()
        raise exception(message)


def initialise(seed=None, xk0dec=None):
    '''Initialise TAUOLA using custom settings'''

    if seed is None:
        seed = ffi.NULL
    else:
        seed = ffi.new('unsigned long [1]', (seed,))

    if xk0dec is None:
        xk0dec = ffi.NULL
    else:
        xk0dec = ffi.new('double [1]', (xk0dec,))

    _call(lib.alouette_initialise, seed, xk0dec)


def decay(mode=None, pid=None, momentum=None, polarisation=None):
    '''Forward Monte Carlo tau decay'''

    if mode is None:
        mode = 0

    if pid is None:
        pid = 15

    if momentum is None:
        momentum = ffi.new('double [3]', (0, 0, 0))
    else:
        if not isinstance(momentum, (list, tuple)):
            momentum = tuple(momentum)
        momentum = ffi.new('double [3]', momentum)

    if polarisation is None:
        polarisation = ffi.NULL
    else:
        if not isinstance(polarisation, (list, tuple)):
            polarisation = tuple(polarisation)
        polarisation = ffi.new('double [3]', polarisation)

    products = ffi.new('struct alouette_products *')
    _call(lib.alouette_decay, mode, pid, momentum, polarisation, products)

    return Products(products)


# User supplied polarisation callback, for backward decays
_POLARISATION = None


@ffi.def_extern()
def _polarisation_callback(pid, momentum, polarisation):
    '''Callback wrapper, for setting the polarisation in backward decays'''
    m = numpy.frombuffer(ffi.buffer(momentum, 24)[:])
    polarisation[0:3] = _POLARISATION(int(pid), m)


class undecay:
    '''Backward Monte Carlo tau decay'''

    @property
    def mother(self):
        '''Mother particle(s) for backward decays'''
        return int(lib.alouette_undecay_mother)

    @mother.setter
    def mother(self, value):
        if value not in (-15, 0, 15):
            exc = ValueError if isinstance(value, Number) else TypeError
            raise exc(
                'bad mother pid. Must be one of 15 (tau-), -15 (tau+) or '
                '0 (tau- or tau+)'
            )
        else:
            lib.alouette_undecay_mother = value

    @property
    def bias(self):
        '''Tuning parameter for the spin bias in backward decays'''
        return float(lib.alouette_undecay_bias)

    @bias.setter
    def bias(self, value):
        if (not isinstance(value, Number)) or (value < -1) or (value > 1):
            exc = ValueError if isinstance(value, Number) else TypeError
            raise exc('bad bias value. Must be within [-1, 1].')
        else:
            lib.alouette_undecay_bias = value

    @property
    def scheme(self):
        '''Scheme for the backward Monte Carlo weight'''
        if lib.alouette_undecay_scheme == lib.ALOUETTE_UNDECAY_CARTESIAN:
            return 'cartesian'
        elif lib.alouette_undecay_scheme == lib.ALOUETTE_UNDECAY_SPHERICAL:
            return 'spherical'
        else:
            return 'energy'

    @scheme.setter
    def scheme(self, value):
        if value == 'cartesian':
            lib.alouette_undecay_scheme = lib.ALOUETTE_UNDECAY_CARTESIAN
        elif value == 'spherical':
            lib.alouette_undecay_scheme = lib.ALOUETTE_UNDECAY_SPHERICAL
        elif value == 'energy':
            lib.alouette_undecay_scheme = lib.ALOUETTE_UNDECAY_ENERGY
        else:
            exc = ValueError if isinstance(value, (str, bytes)) else TypeError
            raise exc('bad scheme value.')

    def __call__(self, mode=None, pid=None, momentum=None, polarisation=None):
        '''Perform a backward Monte Carlo decay'''
        if mode is None:
            mode = 0

        if pid is None:
            pid = 16

        if momentum is None:
            momentum = ffi.new('double [3]', (0, 0, 0))
        else:
            if not isinstance(momentum, (list, tuple)):
                momentum = tuple(momentum)
            momentum = ffi.new('double [3]', momentum)

        if polarisation is None:
            polar_cb = ffi.NULL
        else:
            global _POLARISATION
            _POLARISATION = polarisation
            polar_cb = lib._polarisation_callback

        products = ffi.new('struct alouette_products *')
        _call(
            lib.alouette_undecay,
            mode,
            pid,
            momentum,
            polar_cb,
            products,
        )

        return Products(products)

# Override as class singleton
undecay = undecay()
