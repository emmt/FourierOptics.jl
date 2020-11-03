# Julia package for Fourier optics simulations

[![License](http://img.shields.io/badge/license-MIT-brightgreen.svg?style=flat)](LICENSE.md)
[![Build Status](https://travis-ci.com/emmt/FourierOptics.jl.svg?branch=master)](https://travis-ci.com/emmt/FourierOptics.jl)

This package provides the building blocks for numerical simulations of an
optical system using Fourier optics approximations.

This file provides some documentation:

* [Installation](#installation) of the package;
* [Simple coordinate transforms](#simple-coordinate-transforms);
* [Physical units](#physical-units) used in the package;
* [Physical regions](#physical-regions);


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


## Physical regions

A `Region` defines a sampled rectangular piece of physical plane perpendicular
to the direction of propagation.  The sampling step is the same in all
direction and the grid axes are aligned with the `(x,y)` axes of the world
coordinate system (axis `z` is assumed to be the direction of propagation).

A `Region` can be defined by:

```julia
Region(width, height, step=1m)
```

where `width`, `height` and `step` are the physical dimensions and sampling
step (all lengths in SI units, that is meters) and which yields a region
centered at the origin.  The bounding box (see below) of the resulting region
will be at least of size `width` by `height`.  Another possibility is to use
two instances of `Range` to define the coordinates of the grid of samples along
the two axes:

```julia
Region(X, Y)
```

Because of rounding errors, the effective grid of samples embedded into the
resulting region may have slightly different coordinates.

A region can be *translated* or *recentered*:

```julia
R = Region(X, Y)          # define a region
Rt = R + (tx, ty)         # translate the region
Rc = recenter(R)          # recenter the region (0,0)
Rc = recenter(R, (x0,y0)) # recenter the region at (x0, y0)
```

A number of coordinate transformations can be performed:

```julia
R = Region(X, Y)              # define a region
grid2world(R, (i,j)) -> x, y  # get the world coordinates of a grid node
world2grid(R, (x,y)) -> i, j  # get the grid indices from world coordinates
```

Above, `(x,y)` denotes physical (*world*) coordinates and `(i,j)` denotes,
possibly fractional, grid indices.

Some information can be retrieved from the region instance:

```julia
R = Region(X, Y)  # define a region
length(R)         # yields the number of samples
size(R)           # yields the dimensions (in number of samples) of the region
size(R, i)        # yields the i-th dimension
step(R)           # yields the sampling step (in meters)
grid2world(R)     # yields the grid to world coordinate transformation
world2grid(R)     # yields the world to grid coordinate transformation
extrema(R)        # yields the extreme positions `(xmin,xmax,ymin,ymax)`
boundingbox(R)    # yields the bounding box of the region
center(R)         # yields the physical coordinates of the center of R
```

Compared to `extrema(R)` which yields the extreme world positions `(xmin, xmax,
ymin, ymax)` of the embedded grid nodes, `boundingbox(R)` accounts for an
additional margin of 1/2 sample on all sides.
