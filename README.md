# ALOUETTE
(**E**nvironment for **T**aus and anti-**T**aus **E**ncapsulation with
ta-**UOLA**, backwards)

## Description

ALOUETTE is a backward Monte-Carlo library for tau decays. It is built over
tauola-fortran (v2.9) from the [C++ TAUOLA interface](http://tauolapp.web.cern.ch/tauolapp/).
The library's API provides the following C functions :

```c
/* Initialise the wrapper and TAUOLA. */
enum alouette_return alouette_initialise(int mute, int * seed);

/* Clean the wrapper. */
enum alouette_return tauola_finalise(void);

/* Return a string describing an `alouette_return` code. */
const char * alouette_strerror(enum alouette_return rc);

/* Perform a forward Monte-Carlo tau decay. */
enum alouette_return alouette_decay(int pid, const double momentum[3], const double * polarisation);

/* Perform a backward Monte-Carlo tau decay. */
enum alouette_return alouette_undecay(int pid, const double momentum[3], alouette_polarisation_cb * polarisation, double * weight);

/* Iterator over the decay products. */
enum alouette_return alouette_product(int * pid, double momentum[3]);

/* Getter for the polarimteric vector of the last decay. */
enum alouette_return alouette_polarimetric(double polarimetric[3]);
```

Check the [header](include/alouette.h) file and the provided example files for
detailed usage.

## Portability

This wrapper is meant to be used on a linux box. It relies on POSIX specifics,
e.g. it assumes `/dev/urandom` to be available.

## Thread safety

TAUOLA nor this wrapper do **not** support **multi-threading**.

## License

The ALOUETTE wrapper is  under the **GNU LGPLv3** license. See the provided
`LICENSE` and `COPYING.LESSER` files.
