# TAUOLA-C

## Description

TAUOLA-C is a minimalist C wrapper for the [C++ TAUOLA interface](http://tauolapp.web.cern.ch/tauolapp/). It requires TAUOLA++ to
be installed with HEPEVT enabled, e.g. `./configure --without-hepmc`. The
wrapper provides the following four interface functions to C :

```c
/* Initialise TAUOLA++ through the wrapper. */
void tauola_initialise(int mute, int * seed);
/* Clean the wrapper. */
void tauola_finalise(void);
/* Perform a Monte-Carlo tau decay with TAUOLA. */
void tauola_decay(int pid, double momentum[3], double * polarisation);
/* Iterator over the decay products. */
int tauola_product(int * pid, double momentum[3]);
```

Check the *header* file and the provided *example* file for detailed usage.

## Portability
This wrapper is meant to be used on a linux box. It relies on POSIX specifics,
e.g. `dup2` and it assumes `/dev/urandom` to be available.

## Thread safety
TAUOLA nor this wrapper do **not** support **multi-threading**.

## License
The TAUOLA-C wrapper is  under the **GPL-3.0** license. See the provided
`LICENSE` file.
