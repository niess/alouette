# Introduction

Alouette is a [TAUOLA][TAUOLA] thin wrapper for simulating single $\tau$
decays.  It can operate in forward or in [backward Monte Carlo][BMC] mode.
Alouette is built over [TAUOLA universal interface][tauolapp] (version 1.1.8,
for LHC) Fortran's source. It can be used as a C library (*libalouette*) or as
a Python 3 package ([alouette][alouette_py]).
{: .justify}

## Source and license

Alouette's source is hosted on [GitHub][alouette_source]. It is available under
the terms of the GNU LGPLv3 license. See the provided [LICENSE][license] and
[COPYING.LESSER][lesser] files. The [examples][examples] however have a separate
public domain license allowing them to be copied without any restriction.
{: .justify}

## Issues

Alouette results have been carefully validated through various tests including
comparisons to [Tauola++][tauolapp]. Code updates are unit tested using [GitHub
CI][github_ci].  Yet, whenever you were (un)lucky to come accross a bug, please
consider reporting it as a [GitHub issue][issues].
{: .justify}


[alouette_py]: https://pypi.org/project/alouette
[alouette_source]: https://github.com/niess/alouette
[BMC]: https://arxiv.org/abs/1705.05636
[examples]: https://github.com/niess/alouette/tree/master/examples
[github_ci]: https://github.com/niess/alouette/actions
[issues]: https://github.com/niess/alouette/issues
[license]: https://github.com/niess/alouette/blob/master/LICENSE
[lesser]: https://github.com/niess/alouette/blob/master/COPYING.LESSER
[TAUOLA]: https://www.sciencedirect.com/science/article/pii/001046559190038M
[tauolapp]: http://tauolapp.web.cern.ch/tauolapp/
