# Alouette [![Linux](https://github.com/niess/alouette/actions/workflows/linux.yml/badge.svg)][workflow_linux] [![OSX](https://github.com/niess/alouette/actions/workflows/osx.yml/badge.svg)][workflow_osx]
(Y-**ET** another encapsula**TE**d TA-**UOLA**, backwards)


## Description

Alouette is a [TAUOLA][TAUOLA] thin wrapper for simulating single &tau; decays.
It can operate in forward or in [backward Monte Carlo][BMC] mode.  Alouette is
built over [Tauola++][tauolapp] Fortran source (version 1.1.8, for LHC). It can
be used as a C library (*libalouette*) or as a Python 3 package
([alouette][alouette_py]).

The Alouette API is meant to be light and simple. It exposes two key functions,
*decay* and *undecay*.  Initialisation is automatic.  A more detailed
description of Alouette is avaible online, from [Read the Docs][alouette_docs].
One can also browse the [examples][examples] shipped with the
[source][alouette_source].


## Portability

Alouette is designed for UNIX systems. It relies on POSIX specifics, e.g.  it
assumes `/dev/urandom` to be available. Alouette is tested on
[Linux][workflow_linux] and [OSX][workflow_osx].


## Thread safety

TAUOLA and the Alouette wrapper do **not** support **multi-threading**.


## License

The Alouette wrapper is  under the **GNU LGPLv3** license. See the
[LICENSE][license] and [COPYING.LESSER][lesser] files.


[alouette_docs]: https://alouette.readthedocs.io/en/latest/
[alouette_py]: https://pypi.org/project/alouette
[alouette_source]: https://github.com/niess/alouette
[BMC]: https://arxiv.org/abs/1705.05636
[examples]: https://github.com/niess/alouette/tree/master/examples
[license]: https://github.com/niess/alouette/blob/master/LICENSE
[lesser]: https://github.com/niess/alouette/blob/master/COPYING.LESSER
[TAUOLA]: https://www.sciencedirect.com/science/article/pii/001046559190038M
[tauolapp]: http://tauolapp.web.cern.ch/tauolapp/
[workflow_linux]: https://github.com/niess/alouette/actions/workflows/linux.yml
[workflow_osx]: https://github.com/niess/alouette/actions/workflows/osx.yml
