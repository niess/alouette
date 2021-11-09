from ._core import ffi, lib
import numpy

__all__ = ('decay', 'initialise', 'Products')


class Products:
    '''Container for decay products
    '''

    def __init__(self, c):
        self._c = c

    def __str__(self):
        return f'{type(self).__name__}(size={self.size}, pid={self.pid}, '     \
            f'polarimetric={self.polarimetric}, P={self.P})'

    @property
    def P(self):
        try:
            return self._P
        except AttributeError:
            n = self._c.size
            self._P = numpy.frombuffer(
                ffi.buffer(self._c.P)[:])[:4 * n].reshape((n, 4))
            return self._P

    @property
    def pid(self):
        try:
            return self._pid
        except AttributeError:
            self._pid = numpy.frombuffer(ffi.buffer(self._c.pid)[:],
                dtype='i4')[:self._c.size]
            return self._pid

    @property
    def polarimetric(self):
        try:
            return self._polarimetric
        except AttributeError:
            self._polarimetric = numpy.frombuffer(
                ffi.buffer(self._c.polarimetric)[:])
            return self._polarimetric

    @property
    def size(self):
        return int(self._c.size)

    @property
    def weight(self):
        return float(self._c.weight)


def initialise(xk0dec=None):
    '''Initialise TAUOLA with custom settings
    '''

    if xk0dec is None:
        xk0dec = ffi.NULL
    else:
        xk0dec = ffi.new('double [1]', (xk0dec,))

    lib.alouette_initialise(xk0dec)


def decay(mode=None, pid=None, momentum=None, polarisation=None):
    '''Decay a tau particle with TAUOLA
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
    lib.alouette_decay(mode, pid, momentum, polarisation, products)

    return Products(products)
