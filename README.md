# Alouette
(Y-**ET** another encapsula**TE**d TA-**UOLA**, backwards)

## Description

Alouette is a [TAUOLA][TAUOLA] wrapper allowing to simulate single tau decays.
It can operate in forward or in backward Monte Carlo mode. The Alouette `C`
library is built over `tauola-fortran` from the [TAUOLA universal
interface][tauolapp] (version 1.1.8, for LHC).  Python3 bindings are also
provided as the [`alouette`][alouette_py] package.

The Alouette library exposes the following `C` functions:
```c
/* Perform a forward Monte Carlo tau decay. */
enum alouette_return alouette_decay(
    int mode,
    int pid, const double momentum[3], const double * polarisation,
    struct alouette_products * products);

/* Perform a backward Monte Carlo tau decay. */
enum alouette_return alouette_undecay(
    int mode,
    int pid, const double momentum[3], alouette_polarisation_cb * polarisation,
    double bias,
    struct alouette_products * products);
```

For detailed usage, please refer to the `C` [header][alouette_h] file,
to the [`alouette`][alouette_py] Python package, or to the [examples][examples]
shipped with the source.

## Portability

The Alouette wrapper is meant to be used on UNIX systems. It relies on POSIX
specifics, e.g. it assumes `/dev/urandom` to be available.

## Thread safety

TAUOLA and the Alouette wrapper do **not** support **multi-threading**.

## License

The Alouette wrapper is  under the **GNU LGPLv3** license. See the provided
`LICENSE` and `COPYING.LESSER` files.


[alouette_h]: https://github.com/niess/alouette/blob/master/include/alouette.h
[alouette_py]: https://pypi.org/project/alouette
[examples]: https://github.com/niess/alouette/tree/master/examples
[TAUOLA]: https://www.sciencedirect.com/science/article/pii/001046559190038M
[tauolapp]: http://tauolapp.web.cern.ch/tauolapp/
