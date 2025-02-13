<!--
SPDX-FileCopyrightText: 2024-2025 Helmholtz-Zentrum hereon GmbH
SPDX-FileContributor: Carsten Lemmen <carsten.lemmen@hereon.de>
SPDX-License-Identifier: CC0-1.0
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

This is a the development repository of TAME. TAME is a collection of [FABM](https://fabm.net) models and further utilities.  TAME is currently under development and is not a production-ready product.  

The history of TAME is the Model for Adaptive Ecosystems (MAECS), originally coded by Kai Wirtz and published as Wirtz & Kerimoglu (2016).  It is a variation of a nutrient-phytoplankton-zooplankton-detritutus (NPZD) ecosystem model, but different from other ecosystem models, it featured 
- the physiological adaptation of allocation within phytoplankton
- optional loss terms by virus

MAECS was used extensively in the FABM version 0 framework for Southern North Sea applications in the MOSSCO project (Wirtz 2019; Slavik et al. 2019; Nasermoaddeli et al. 2018; Lemmen 2018).  The necessity to move to a new model resulted from
- difficulties of porting to FABM version 1 and above
- strong reliance on namelists
- complexity of the code 

TAME will try to remedy the above difficulties, and it's concept is from the beginning targeted towards flexibilty, user friendlyness, and compatibility with the current version of the FABM framework.

# Obtaining and operating TAME

All of TAME and its dependencies are free and open source software. To operate it locally on your computer, you need (1) a compilation toolchain and (2) the sources of TAME, FABM, and GOTM.

## Compiler toolchain

TAME is written in the programming language Fortran and uses the NetCDF binary self-documenting output format.  So you need a Fortran compiler and a compatible installation of the NetCDF developer library.  The sources are hosted on a `git` repository, and the build is facilitated by `cmake`.  The easiest way is to install those for your system in a package manager, such as `apt`, `brew`, or `conda`.

Typically on Ubuntu linux (and others) you would use `apt`
```
apt install gfortran libnetcdff-dev cmake git
```

Typically on Linux or Mac with the homebrew package manager you would use `brew`
```
brew install gcc netcdf-fortran cmake git
```

Typically on Linux or Mac with conda package manager you would use `conda` or `mamba`
```
conda create -n tame python=3.12
conda activate tame
conda install gfortran netcdf-fortran cmake git
```

## Obtaining source codes

TAME itself is a model within the FABM framework, which relies on GOTM as a hydrodynamic driver.  So you need to get download three software packages to get this to work.  Before you do, decide on a directory where you install your sources and define the environment variables `GOTM_BASE`, `FABM_BASE`, and `TAME_BASE` to point to these.  We will here assume they are all located in your `$HOME/devel/` folder.

```
export TAME_BASE=$HOME/devel/tame
export FABM_BASE=$HOME/devel/fabm
export GOTM_BASE=$HOME/devel/gotm

git clone https://codebase.helmholtz.cloud/kse/generalized-aquatic-ecosystem-model.git $TAME_BASE
git clone https://github.com/fabm-model/fabm.git $FABM_BASE
git clone https://github.com/gotm-model/code.git $GOTM_BASE
```

As the 0d driver needs a stable version of GOTM, change to this one and also recursively download all
GOTM dependencies
```
(cd $GOTM_BASE && git checkout v6.0 && git submodule update --init --recursive)
```

## How to build TAME

As TAME is a FABM model, we rely on the build structure implemented in FABM and give it separate arguments to also consider TAME in its model hierarchy.  Therefore, the source directory we supply to `CMake` is `$FABM_BASE`.  And it is good style to do an out-of-tree build, i.e., to compile in a directory separate from your sources.  We will define this as `$BUILD_GOTM`, as GOTM is the default FABM host

```
export BUILD_GOTM=$HOME/devel/build-gotm
mkdir $BUILD_GOTM

cmake -B $BUILD_GOTM -S $FABM_BASE -DFABM_HOST=gotm -DFABM_BASE=$FABM_BASE -DFABM_INSTITUTES=tame -DFABM_TAME_BASE=$TAME_BASE -DGOTM_BASE=$GOTM_BASE
make
```

For more information on the FABM `CMake` build, consult their [building and installing](https://github.com/fabm-model/fabm/wiki/Building-and-installing) page. 

Note that `-DFABM_INSTITUTES=tame` will make FABM compile our models as the _only_ available biogeochemical models. If you additionally want to have access to other biogeochemical models included with FABM, you can set `FABM_INSTITUTES` to a semi-colon separated list, e.g., `-DFABM_INSTITUTES="tame;vims;iow"` (to prevent the shell from interpreting the semi-colons, you typically have to enclose this list with quotes).

You may also want to build the 0d model, which you achieve by changing the source directory to `CMake`

```
export BUILD_0d=$HOME/devel/build-0d
mkdir $BUILD_0d

cmake -B $BUILD_0d -S $FABM_BASE/src/drivers/0d -DFABM_HOST=0d -DFABM_BASE=$FABM_BASE -DFABM_INSTITUTES=tame -DFABM_TAME_BASE=$TAME_BASE -DGOTM_BASE=$GOTM_BASE
make
```

## How to build the parser

The parser is contained in the folder [./parser](./parser). It is a simple C program and can be compiled with any suitable C compiler.

There is also a `python` template based parser in development in the [./python](./python) folder

## Report when it is working/not working for you!

Build test with build chain suggested above

| System | Type   | Compiler      | Status |
| ------ | ------ | ------------- | ------ |
| quoll  | arm64  | gfortran      |        |
| femto  | x86-64 | ifort         |        |
| kuro   | x86-64 | intel/openmpi |        |
| strand | x86-64 | ifort         |        |

## License information

Copyright 2024-2025 Helmholtz-Zentrum hereon GmbH

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

License management is handled with [`reuse`](https://reuse.readthedocs.io/).
If you have any questions on this, please have a look into the
[contributing guide][contributing] or contact the maintainers of TAME.

[contributing]: https://generalized-aquatic-ecosystem-model.readthedocs.io/en/latest/contributing.html
