# StatFEM_EUCLID

This package provides tooling for unsupervised constitutive model discovery based on statFEM an EUCLID.
Instead of 'hard-wiring' a single finite element software, the package is able to use any solver as a black box
through a [UM-Bridge](https://um-bridge-benchmarks.readthedocs.io/en/docs/) interface.

## Features
- A high(er)-level interface to UM-Bridge for sampling an FEM black box
- Construction of PCE surrogates from samples
- Data assimilation based on the statFEM method (WIP: in refactoring)
- Parameter identification based on the virtual fields method or EUCLID (WIP: in refactoring)


## License
StatFEM_EUCLID.jl is licensed under the MIT License. 