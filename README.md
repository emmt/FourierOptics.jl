# FourierOptics [![Build Status](https://github.com/emmt/FourierOptics.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/emmt/FourierOptics.jl/actions/workflows/CI.yml?query=branch%3Amain) [![Build Status](https://ci.appveyor.com/api/projects/status/github/emmt/FourierOptics.jl?svg=true)](https://ci.appveyor.com/project/emmt/FourierOptics-jl) [![Coverage](https://codecov.io/gh/emmt/FourierOptics.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/emmt/FourierOptics.jl)

Package `FourierOptics` provides tools based on Fourier optics to simulate the
propagation of a complex amplitude in an optical system.

`FourierOptics` was inspired by the [`PROPER`
library](https://sourceforge.net/projects/proper-library/) by John E. Krist
(see References below).

## Usage

The (scalar) field at a given position along the direction `z` of propagation is stored in a `Field` structure.

Methods with a band prefix `!` modify the field in-place, methods without a `!`
suffix yield a modified independent copy of the input `Field` structure.

## Links

* [PROPER](https://sourceforge.net/projects/proper-library/), IDL, Matlab, and
  Python versions.
* [LightPipes](http://www.okotech.com/lightpipes), [Python version](https://github.com/opticspy/lightpipes).

## References

1. John E. Krist, [*PROPER: An Optical Propagation Library for
   IDL*](doi:10.1117/12.731179), in *Optical Modeling and Performance
   Predictions III*, edited by Mark A. Kahan, Proc. of SPIE Vol. **6675**,
   66750P, (2007)·

2. J. W. Goodman, **Fourier Optics**, 2nd edition (1996).
