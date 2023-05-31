# Correa 
[![Build bindings](https://github.com/yossibokor/correa/actions/workflows/bindings.yml/badge.svg?branch=master)](https://github.com/yossibokor/correa/actions/workflows/bindings.yml)
[![Executables](https://github.com/yossibokor/correa/actions/workflows/executables.yml/badge.svg?branch=master)](https://github.com/yossibokor/correa/actions/workflows/executables.yml)
[![Documentation](https://github.com/yossibokor/correa/actions/workflows/documentation.yml/badge.svg?branch=master&event=push)](https://github.com/yossibokor/correa/actions/workflows/documentation.yml)

## UNDER DEVELOPMENT: This is a public beta. Please use [issues](https://github.com/yossibokor/correa/issues) to report any bugs

Correa is a C++ library, you can choose to build Python bindings and/or executables. 

### License

Correa is released under the [BSD 3-Clause License](md_LICENSE.html), which you should have received a copy of.


## Installation

First clone the [repository](https://github.com/yossibokor/correa) with recursive submodules, and ensure you have [CMake](https://cmake.org/), [BLAS](https://www.netlib.org/blas/) and [LAPACK](https://www.netlib.org/lapack/) installed.

Then, use [cmake](https://cmake.org/) from the directory to build Correa with the options you would like. The options are
- `CORREA_BUILD_EXECUTABLES` (default OFF)
- `CORREA_BUILD_PYTHON_BINDINGS` (default ON)


To build the Python bindings:
```bash
git clone --recurse-submodules https://github.com/yossibokor/correa.git
cd correa
cmake -S . -B build
cd build
make
```

and then ether launch python from `.../correa/build` or add `.../correa/build` to your `PYTHONPATH`.

## UDER DEVELOPMENT: This is a public beta. Please use [issues](https://github.com/yossibokor/correa/issues) to report any bugs.
