# ALOUETTE
(Y-**ET** another encapsula**TE**d TA-**UOLA**, backwards)

## Description

ALOUETTE is a C wrapper over [TAUOLA][TAUOLA] allowing to perform both forward
and backward Monte-Carlo tau decays. It is built over tauola-fortran
distribution from the [TAUOLA universal
interface](http://tauolapp.web.cern.ch/tauolapp/) (1.1.8, LHC).  The library API
provides the following C functions:

```c
/* Perform a forward Monte-Carlo tau decay. */
enum alouette_return alouette_decay(
    int mode,
    int pid, const double momentum[3], const double * polarisation,
    struct alouette_products * products);

/* Perform a backward Monte-Carlo tau decay. */
enum alouette_return alouette_undecay(
    int mode,
    int pid, const double momentum[3], alouette_polarisation_cb * polarisation,
    double bias,
    struct alouette_products * products);
```

Check the [header](include/alouette.h) file and the provided example files for
more detailed usage.

## Portability

This wrapper is meant to be used on UNIX systems. It relies on POSIX specifics,
e.g. it assumes `/dev/urandom` to be available.

## Thread safety

TAUOLA and the ALOUETTE wrapper do **not** support **multi-threading**.

## License

The ALOUETTE wrapper is  under the **GNU LGPLv3** license. See the provided
`LICENSE` and `COPYING.LESSER` files.


[TAUOLA]: https://www.sciencedirect.com/science/article/pii/001046559190038M
