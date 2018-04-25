# Julia package for Fourier optics simulations

[![License](http://img.shields.io/badge/license-MIT-brightgreen.svg?style=flat)](LICENSE.md)
[![Build Status](https://travis-ci.org/emmt/FourierOptics.jl.svg?branch=master)](https://travis-ci.org/emmt/FourierOptics.jl)

This package provides the building blocks for numerical simulations of an
optical system using Fourier optics approximations.

This file provides some documentation:

* [Installation](#installation) of the package;
* [Simple coordinate transforms](#simple_coordinate_transforms);
* [Physical units](#physical_units) used in the package;


## Installation

[**FourierOptics**](https://github.com/emmt/FourierOptics.jl) is not yet an
[offical Julia package](https://pkg.julialang.org/) but it is easy to install
from the [julia REPL](https://docs.julialang.org/en/stable/manual/interacting-with-julia/) with one of the following commands (depending which of *https* or *ssh* is the most suitable for you):

```julia
# for https:
Pkg.clone("https://github.com/emmt/FourierOptics.jl.git")
# or for ssh:
Pkg.clone("git@github.com:emmt/FourierOptics.jl.git")
```

note that there is no needs to call `Pkg.build("FourierOptics")` after.

To check whether the **FourierOptics** package works correctly:

```julia
Pkg.test("FourierOptics")
```

To update to the last version, just type:

```julia
Pkg.update()
```

and perhaps test again...


## Simple coordinate transforms

A simple 2D coordinate transform is created by:

```julia
using FourierOptics
C = CoordinateTransform(stp, x0, y0)
```

which yields `C` such that:

```julia
C(u, v) -> (x0 + stp*u, y0 + stp*v)
```

Called without arguments, the identity is returned, that is
`CoordinateTransform()` is the same as `CoordinateTransform(1,0,0)`.

The parameters of the coordinate transform `C` can be retrieved by:

```julia
step(C)   -> stp
origin(C) -> (x0, y0)
```

Input coordinates can be a 2-component tuple or vector (note that the
result is simlar to the argument):

```julia
C((u, v)) -> (x, y)
C([u, v]) -> [x, y]
```

Coordinate transforms may be inverted, pre/post-scaled, pre/post-translated:

```julia
A = inv(C)        # A(C(u,v)) = (u,v) and C(A(x,y)) = (x,y)
A = ρ*C           # A(u, v) = ρ*C(u, v)
A = C*ρ           # A(u, v) = C(ρ*u, ρ*v)
A = C + (tu, tv)  # A(u, v) = C(u + tu, v + tv)
A = C + [tu, tv]  # A(u, v) = C(u + tu, v + tv)
A = (tx, ty) + C  # A(u, v) = C(u, v) +. (tx, ty)
A = [tx, ty] + C  # A(u, v) = C(u, v) +. [tx, ty]
```

where `ρ` is some scalar factor.  The `-` operator can be used like `+` to
apply the opposite pre/post translation.

Two coordinate transforms `A` and `B` can be composed as follows:

```julia
C = A ∘ B
```

to form a new coordinate transform, `C` such that `C(u,v) = A(B(u,v))`.
Operators `\` and  `/` can be used for left/right *division*:

```julia
A\B = inv(A) ∘ B
A/B = A ∘ inv(B)
```

Two coordinate transforms `A` and `B` can be approximately compared:

```julia
A ≈ B
```

Using all this, coordinate transforms can be built starting with the identity:

```julia
I = CoordinateTransform()
C = (x0, y0) + stp*I
```

which may be more readable in the code.


## Physical units

In the [FourierOptics](https://github.com/emmt/FourierOptics.jl) module, all
physical units (only lengths and angles are needed) are assumed to be in SI
units (hence lengths are in meters, angles are in radians, *etc.*).

To use short symbolic units names, just do:

```julia
using FourierOptics.Units
```

to be able to define physical parameters with a readable syntax.  For instance
`λ = 500nm` or `λ = 0.5µm` to define a wavelength, or `f = 15.2mm` to define a
focal length.

Exported units are:

* lengths: `m` (meters), `cm` (centimeters), `mm` (millimeters), `µm`
  (micrometers), `nm` (nanometers);

* angles: `rd` (radians), `deg` or `arcdeg` (degrees of arc), `arcmin` (minutes
  of arc), `arcsec` (seconds of arc), `marcsec` or `mas` (milliarcseconds),
  `µarcsec` or `µas` (microarcseconds).
