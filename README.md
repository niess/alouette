# ALOUETTE
(Y-**ET** another encapsula**TE**d TA-**UOLA**, backwards)

## Description

ALOUETTE is a backward Monte-Carlo library for tau decays. It is built over
tauola-fortran (v2.9, LHC) from the [TAUOLA universal
interface](http://tauolapp.web.cern.ch/tauolapp/).  The library API provides
the following C functions:

```c
/* Perform a forward Monte-Carlo tau decay. */
enum alouette_return alouette_decay(
    int pid, const double momentum[3], const double * polarisation);

/* Perform a backward Monte-Carlo tau decay. */
enum alouette_return alouette_undecay(int pid, const double momentum[3],
    alouette_polarisation_cb * polarisation, double bias, double * weight);

/* Iterator over the decay products. */
enum alouette_return alouette_product(int * pid, double momentum[3]);

/* Getter for the polarimteric vector of the last decay. */
enum alouette_return alouette_polarimetric(double polarimetric[3]);
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
