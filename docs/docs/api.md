# Alouette API

The Alouette API is meant to be light and simple. It exposes two key functions,
*decay* and *undecay*. Initialisation is automatic, but it can also be done
explicitly if custom settings are needed. Basic [examples][examples] of usage
are shipped with the [source][alouette_source].
{: .justify}

!!! note
    The C and Python API are almost identical. Alouette C library symbols are
    prefixed with `alouette_`, while Python symbols reside in the `alouette`
    package.
    {: .justify}

!!! note
    Alouette uses [PDG IDs][PID] (PIDs) for labeling particles. I.e.  a
    $\tau$-lepton has `pid = 15` while its anti-particle, a $\tau^+$, is
    numbered `-15`.
    {: .justify}

!!! note
    Alouette uses the same natural system of units than TAUOLA. E.g. momenta
    are expressed in $GeV / c$.
    {: .justify}


<div markdown="1" class="shaded-box fancy">
## Decay function

The *decay* function simulates a $\tau$ decay with TAUOLA. An optionnal
polarisation 3-vector can be provided. If the latter is `NULL` (`None`), then
spin effects are not simulated. The format of [decay products](#decay-products)
is described below.
{: .justify}

!!! note
    $\tau$ decay modes are indexed according to TAUOLA, see e.g. Appendix C of
    [Davidson et al. (2012)][Davidson2012] for a description of available modes.
    If *mode* is set to zero, then the decay is randomised over all modes,
    according to branching ratios.
    {: .justify}

### C synopsis

```C
enum alouette_return alouette_decay(int mode, int pid, const double momentum[3],
    const double * polarisation, struct alouette_products * products);
```

### Python synopsis

```Python
alouette.decay(mode=None, pid=None, momentum=None, polarisation=None)
```

!!! note
    Default is `mode=0`, `pid=15` and `momentum=(0,0,0)`, I.e. a center of mass
    $\tau^-$ decay considering all modes.
    {: .justify}
</div>


<div markdown="1" class="shaded-box fancy">
## Undecay function

The *undecay* function simulates a $\tau$ decay from a given decay product using
the [Backward Monte Carlo][BMC] technique (BMC). The spin polarisation of the
primary $\tau$ can be provided a posteriori as a callback function.  Setting
*polarisation* to `NULL` (`None`) results in spin effects to be ignored. The
format of [decay products](#decay-products) is described below.
{: .justify}

!!! note
    $\tau$ decay modes are indexed according to TAUOLA, see e.g. Appendix C of
    [Davidson et al. (2012)][Davidson2012] for a description of available modes.
    If *mode* is set to zero, then the decay is randomised over all modes,
    according to branching ratios.
    {: .justify}

!!! warning
    The daughter *momentum* must be non null in a backward decay. I.e. it is not
    possible to backward simulate center of mass decays.
    {: .justify}

### Configuration parameters

The undecay function has two additional configuration parameters, for advanced
usage. The `alouette_undecay_mother` (`undecay.mother`) parameter allows to
specify the mother PID in a backward decay. Setting it to 0 results in both
$\tau^-$ and $\tau^+$ being considered, which is the default setting.
{: .justify}

The `alouette_undecay_bias` (`undecay.bias`) parameter allows to control the
spin biasing procedure, in backward decays. It must be within $[-1, 1]$.
Default is 1, which assumes left (right) handed $\tau^-$ ($\tau^+$).
{: .justify}

### C synopsis

```C
enum alouette_return alouette_undecay(int mode, int pid,
    const double momentum[3], alouette_polarisation_cb * polarisation,
    struct alouette_products * products);

typedef void alouette_polarisation_cb(
    int pid, const double momentum[3], double * polarisation);
```

### Python synopsis

```Python
alouette.undecay(mode=None, pid=None, momentum=None, polarisation=None)
```

!!! note
    Default is `mode=0` and `pid=16`, a $\nu_\tau$ decay product considering all
    modes.
    {: .justify}
</div>


<div markdown="1" class="shaded-box fancy">
## Decay products

Decay products are stored in a dedicated data structure. The daughters PIDs and
4-momenta, *P*, are accessed as arrays. The decay *polarimeter* vector is also
provided for advanced usage, e.g. in order to simulate spin-spin corellations in
$\tau^-\tau^+$ pairs, following [Jadach et al.  (1991)][Jadach1991].
{: .justify}

!!! note
    In backward mode, the mother particle is contained in the "decay products"
    as first entry, while the initial daughter is suppressed. In addition, the
    *weight* field indicates the corresponding BMC weight.
    {: .justify}

### C synopsis

```C
struct alouette_products {
        /** Number of decay products. */
        int size;
        /** PDG IDs of decay products. */
        int pid[ALOUETTE_MAX_SIZE];
        /** Four momenta (px, py, pz, E) of decay products, in GeV/c. */
        double P[ALOUETTE_MAX_SIZE][4];
        /** Polarimeter vector of the decay. */
        double polarimeter[3];
        /** Monte Carlo weight of the decay. */
        double weight;
```

### Python synopsis

```Python
class alouette.Products:
    size: int
    pid: int
    P: numpy.ndarray
    polarimeter: numpy.ndarray
    weight: float
```

!!! warning
    All attributes are read-only, including the 4-momenta, *P*, and the decay
    polarimeter.
    {: .justify}
</div>


<div markdown="1" class="shaded-box fancy">
## Errors handling

The C library functions return a status code indicating success or failure, as
described [below](#c-return-codes). When an error occurs, the `alouette_message`
function can be used in order to fetch the last error message, as a C string.
{: .justify}

In the Python package, errors are managed automatically. Internal C library
errors result in raising a `ValueError` or a `RuntimeError` exception, with
the corresponding error message.
{: .justify}

### C return codes

```C
enum alouette_return {
        /** Execution was successful. */
        ALOUETTE_RETURN_SUCCESS = 0,
        /** A parameter is out of its validity range. */
        ALOUETTE_RETURN_VALUE_ERROR,
        /** A TAUOLA error occured. */
        ALOUETTE_RETURN_TAUOLA_ERROR
};

/* Get the last (error) message(s). */
const char * alouette_message(void);
```
</div>


<div markdown="1" class="shaded-box fancy">
## Initialisation

Initialisation is done automatically by Alouette using default settings.
However, for advanced usage, the library can also be initialised explicitly
using the `alouette_initialise` (`alouette.initialise`) function. Then, one can
specify an optionnal random seed for TAUOLA's "warm up", i.e. the determination
of $W_\text{max}$ weights, as detailed in [Jadach et al. (1991)][Jadach1991].
In addition, one can override the cutoff, $k_0^\text{decay}$, for radiative
corrections in leptonic decays, see e.g. [Jezabek et al. (1992)][Jezabek1992].
{: .justify}

!!! note
    Providing a `NULL` (`None`) value results in default setting to be used for
    the corresponding parameter.
    {: .justify}

### C synopsis

```C
enum alouette_return alouette_initialise(unsigned long * seed, double * xk0dec);
```

### Python synopsis

```Python
alouette.initialise(seed=None, xk0dec=None)
```
</div>


<div markdown="1" class="shaded-box fancy">
## Random stream

Alouette embeds a Mersenne Twister Pseudo Random Numbers Generator (PRNG). The
next number in the sequence is returned by the `alouette_random`
(`alouette.random`) function. The PRNG can be reset using `alouette_random_set`
(`random.set`), and providing an optionnal seed. The `alouette_random_seed`
function (`random.seed` attribute) returns the current seed value.
{: .justify}

!!! note
    If a `NULL` (`None`) seed is provided, then the PRNG is seeded using the
    OS entropy, i.e. `/dev/urandom`.
    {: .justify}

!!! note
    C users can use their own PRNG by overriding the `alouette_random` function
    pointer. For Python users, the PRNG cannot be modified.
    {: .justify}

### C synopsis

```C
/* The library PRNG, uniform over (0,1). */
extern float (*alouette_random)(void);

/* Get the random seed of the built-in PRNG. */
unsigned long alouette_random_seed(void);

/* (Re)set the built-in PRNG. */
void alouette_random_set(unsigned long * seed);
```

### Python synopsis

```Python
# The library PRNG, uniform over (0,1).
alouette.random()

# Get the random seed of the built-in PRNG.
alouette.random.seed

# (Re)set the built-in PRNG.
alouette.random.set(seed=None)
```
</div>


<div markdown="1" class="shaded-box fancy">
## Tauola data

For more advanced usage, one might need to access TAUOLA's internal data. Those
are exposed to C users with the `tauola_` prefix. E.g., the `/PARMAS/` common
block can be accessed as `tauola_parmas` C structure. For Python users, some of
these data are also wrapped in `alouette.tauola` submodule. They can be
accessed as attributes, using the usual syntax. E.g. `tauola.parmas.amtau`
returns the (read-only) $\tau$-lepton mass, according to TAUOLA.
{: .justify}

!!! note
    The `tauola.h` header file, shipped with alouette, provides C definitions
    for some relevant TAUOLA data.
    {: .justify}
</div>


[alouette_source]: https://github.com/niess/alouette
[BMC]: https://arxiv.org/abs/1705.05636
[Davidson2012]: http://doi.org/10.1016/j.cpc.2011.12.009
[examples]: https://github.com/niess/alouette/tree/master/examples
[Jadach1991]: http://doi.org/10.1016/0010-4655(91)90038-m
[Jezabek1992]: https://doi.org/10.1016/0010-4655(92)90092-D
[PID]: https://pdg.lbl.gov/2021/reviews/rpp2020-rev-monte-carlo-numbering.pdf
