# Installation

## The easy way

Binary distributions of Alouette are avaible from [PyPI][alouette_py]. The
Python 3 package and the C library can be installed using pip, e.g. as:
{: .justify}

```bash
pip3 install alouette
```

!!! note
    The C library, *libalouette*, is shipped with the [alouette][alouette_py]
    Python 3 package. Installing the Python package also installs
    *alouette-config*, a configuration script that simplifies linking against
    the C library, e.g. using `$(alouette-config --cflags)` and
    `$(alouette-config --libs)`.
    {: .justify}

## From source

Alouette can be built from [source][alouette_source] using the provided
[Makefile][alouette_make], e.g. as:
{: .justify}

```bash
make CC=gcc
make package CC=gcc PYTHON=python3
```

where the second line builds the Python 3 package. The optionnal `CC` and
`PYTHON` environment variables specify the compiler and the Python interpreter
used for the build process.
{: .justify}

!!! warning
    A Fortran 95 compatible compiler is required in order to build Alouette from
    source, e.g. gfortran. However, there is no runtime dependency on Fortran.
    {: .justify}

[alouette_make]: https://github.com/niess/alouette/blob/master/Makefile
[alouette_py]: https://pypi.org/project/alouette
[alouette_source]: https://github.com/niess/alouette
