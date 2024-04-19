<!--
SPDX-FileCopyrightText: 2024 Helmholtz-Zentrum hereon GmbH
SPDX-FileContributor: Carsten Lemmen <carsten.lemmen@hereon.de>
SPDX-License-Identifier: CC-BY-4.0
-->

[![CI](https://codebase.helmholtz.cloud/kse/generalized-aquatic-ecosystem-model/badges/main/pipeline.svg)](https://codebase.helmholtz.cloud/kse/generalized-aquatic-ecosystem-model/-/pipelines?page=1&scope=all&ref=main)
[![Code coverage](https://codebase.helmholtz.cloud/kse/generalized-aquatic-ecosystem-model/badges/main/coverage.svg)](https://codebase.helmholtz.cloud/kse/generalized-aquatic-ecosystem-model/-/graphs/main/charts)
<!-- TODO: uncomment the following line when the package is registered at https://readthedocs.org -->
<!-- [![Docs](https://readthedocs.org/projects/generalized-aquatic-ecosystem-model/badge/?version=latest)](https://generalized-aquatic-ecosystem-model.readthedocs.io/en/latest/) -->
[![Latest Release](https://codebase.helmholtz.cloud/kse/generalized-aquatic-ecosystem-model/-/badges/release.svg)](https://codebase.helmholtz.cloud/kse/generalized-aquatic-ecosystem-model)
<!-- TODO: uncomment the following line when the package is published at https://pypi.org -->
<!-- [![PyPI version](https://img.shields.io/pypi/v/generalized-aquatic-ecosystem-model.svg)](https://pypi.python.org/pypi/generalized-aquatic-ecosystem-model/) -->
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
[![Imports: isort](https://img.shields.io/badge/%20imports-isort-%231674b1?style=flat&labelColor=ef8336)](https://pycqa.github.io/isort/)
[![PEP8](https://img.shields.io/badge/code%20style-pep8-orange.svg)](https://www.python.org/dev/peps/pep-0008/)
[![Checked with mypy](http://www.mypy-lang.org/static/mypy_badge.svg)](http://mypy-lang.org/)
<!-- TODO: uncomment the following line when the package is registered at https://api.reuse.software -->
<!-- [![REUSE status](https://api.reuse.software/badge/codebase.helmholtz.cloud/kse/generalized-aquatic-ecosystem-model)](https://api.reuse.software/info/codebase.helmholtz.cloud/kse/generalized-aquatic-ecosystem-model) -->


# Trait-based Adaptive Ecosystem Model (TAME)
This is a the development repository of TAME.  TAME is a collection of [FABM](https://fabm.net) models and further utilities


## How to build the parser

The parser is contained in the folder [./parser](./parser).  It is a simple C program and can be compiled with any suitable C compiler.

There is also a `python` template based parser in development in the [./python](./python) folder

## How to build the FABM model

The FABM must be compiled together with FABM. To do this, provide the following additional arguments to cmake [when you build FABM](https://github.com/fabm-model/fabm/wiki/Building-and-installing): `-DFABM_INSTITUTES=tame -DFABM_TAME_BASE=</path/to/tame>`

Here, `</path/to/tame>` is this diretory, i.e. the same directory that contains this ReadMe file. Note that `-DFABM_INSTITUTES=tame` will make FABM compile our models as the *only* available biogeochemical models. If you additionally want to have access to other biogeochemical models included with FABM, you can set `FABM_INSTITUTES` to a semi-colon separated list, e.g., `-DFABM_INSTITUTES="tame;vims;iow"` (to prevent the shell from interpreting the semi-colons, you typically have to enclose this list with quotes).

## Report when it is working/not working for you!

Build test with build chain suggested above

| System | Type     | Compiler | Status |
|--------|----------|----------|--------|
| quoll  | arm64   | gfortran |        |
| femto  | x86-64   | ifort    |      |
| kuro   | x86-64   | intel/openmpi    |      |
| strand | x86-64   | ifort    |        |


## License information

Copyright © 2024 Helmholtz-Zentrum hereon GmbH

Code files in this repository are licensed under the
GPL-3.0-or-later, if not stated otherwise
in the file.

Documentation files in this repository are licensed under CC-BY-4.0, if not stated otherwise in the file.

Supplementary and configuration files in this repository are licensed
under CC0-1.0, if not stated otherwise
in the file.

Please check the header of the individual files for more detailed
information.


### License management

License management is handled with [``reuse``](https://reuse.readthedocs.io/).
If you have any questions on this, please have a look into the
[contributing guide][contributing] or contact the maintainers of TAME.

[contributing]: https://generalized-aquatic-ecosystem-model.readthedocs.io/en/latest/contributing.html
