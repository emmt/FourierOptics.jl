#
# regions.jl -
#
# Implment physical regions.
#

module Regions

export
    Region,
    boundingbox,
    grid2world,
    world2grid

import ...FourierOptics: center, recenter
using ..CoordinateTransforms
using ..Units

"""

A `Region` defines a rectangular sampled piece of physical plane perpendicular
to the direction of propagation.  The sampling step is the same in all
direction and the grid axes are aligned with the `(x,y)` axes of the world
coordinate system (axis `z` is assumed to be the direction of propagation).

A `Region` can be defined by:

```julia
Region(width, height, step=1mm)
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

Finally, the following keyword-only constructor is available:

```julia
Region(; size=(nx,ny), step=stp, center=(xc,yc))
```

where `nx` and `ny` are the number of grid nodes along the dimensions of the
region (both are `100` by default), `stp` is the sampling step size (`1mm` by
default) and `(xc,yc)` are the coordinates of the geometrical center of the
region (`(0mm,0mm)` by default).

Because of rounding errors, the effective grid of samples embedded into the
resulting region may have slightly different coordinates.

A region can be translated or recentered:

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

Here, `(x,y)` denotes physical (*world*) coordinates and `(i,j)` denotes,
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

"""
struct Region
    sz::NTuple{2,Int}        # number of samples in each direction
    i2w::CoordinateTransform # indices to world coordinate transform
    w2i::CoordinateTransform # world to indices coordinate transform
    Region(sz::NTuple{2,Integer}, i2w::CoordinateTransform) =
        (isfinite(step(i2w)) && step(i2w) != 0 ||
         throw(ArgumentError("invalid sampling step"));
         new(sz, i2w, inv(i2w)))
end
Base.size(R::Region) = R.sz
Base.size(R::Region, i) = R.sz[i]
Base.length(R::Region) = prod(size(R))
grid2world(R::Region) = R.i2w
grid2world(R::Region, args...) = R.i2w(args...)
world2grid(R::Region) = R.w2i
world2grid(R::Region, args...) = R.w2i(args...)
Base.step(R::Region) = step(grid2world(R))

Base.show(io::IO, ::MIME"text/plain", R::Region) = show(io, R)

function Base.show(io::IO, R::Region)
    print(io, "Region(size = $(size(R)), step = ")
    showlength(io, step(R))
    print(io, ", center = ")
    showposition(io, center(R))
    print(io, ")")
end

# Keyword-only constructor (to look like the output of the show method).
function Region(;
                size::NTuple{2,Integer} = (100,100),
                step::Real = 1mm,
                center::NTuple{2,Union{Real,Rational,Irrational}} = (0mm, 0mm))
    nx, ny = size
    xc, yc = center
    dx = step*(nx - 1)/2
    dy = step*(ny - 1)/2
    return Region(linspace(xc - dx, xc + dx, nx),
                  linspace(yc - dy, yc + dy, ny))
end

function Region(width::Real, height::Real, step::Real=1mm)
    # Build a large enough region centered at (0,0).
    n1 = ceil(Int, width/step)
    n2 = ceil(Int, height/step)
    x0 = -(n1 + 1)*step/2
    y0 = -(n2 + 1)*step/2
    return Region((n1, n2), CoordinateTransform(step, x0, y0))
end

function Region(X::Range{Tx}, Y::Range{Ty}) where {Tx<:Real,Ty<:Real}
    # Get the sampling step.
    xstp, ystp = step(X), step(Y)
    T = float(promote_type(Tx, Ty))
    if abs(xstp - ystp) > 4*eps(T)*max(abs(xstp), abs(ystp))
        throw(ArgumentError("steps must be (approximately) the same"))
    end
    stp = (xstp + ystp)/2

    # Get the number of samples along each side.
    nx, ny = length(X), length(Y)

    # Compute the offsets so as to keep the central position.
    x0 = ((first(X) - stp) + (last(X) - nx*stp))/2
    y0 = ((first(Y) - stp) + (last(Y) - ny*stp))/2

    # Build the region using the inner constructor.
    return Region((nx, ny), CoordinateTransform(stp, x0, y0))
end

"""
```julia
center(R) -> xc, yc
```

yields the physical coordinates `(xc,yc)` of the geometrical center of the
region `R`.

"""
center(R::Region) = ((n1, n2) = size(R);
                     grid2world(R, (1 + n1)/2, (1 + n2)/2))

"""
```julia
recenter(R, (x0, y0) = (0,0))
```

yields a new physical region of same physical size and sampling step as `R` but
whose center is at `(x0,y0)`.

"""
function recenter(R::Region, c::NTuple{2,Real} = (0,0))
    n1, n2 = size(R)
    stp = step(R)
    x0 = c[1] - (1 + n1)*stp/2
    y0 = c[2] - (1 + n2)*stp/2
    return Region((n1, n2), CoordinateTransform(stp, x0, y0))
end

recenter(R::Region, x0::Real, y0::Real) =
    recenter(R, (x0, y0))

# Translate the region.
Base.:+(R::Region, t::NTuple{2,Real}) =
    ((x0, y0) = origin(grid2world(R));
     Region(size(R), CoordinateTransform(step(R), x0 + t[1], y0 + t[2])))

Base.:-(R::Region, t::NTuple{2,Real}) =
    ((x0, y0) = origin(grid2world(R));
     Region(size(R), CoordinateTransform(step(R), x0 - t[1], y0 - t[2])))

"""
```julia
boundingbox(R) -> xmin, xmax, ymin, ymax
```

yields the physical size of the region `R`.  Compared to `extrema(R)`, an
extra margin of 1/2 sample is taken into account.

"""
function boundingbox(R::Region)
    n1, n2 = size(R)
    x0, y0 = grid2world(R, 0.5, 0.5)
    x1, y1 = grid2world(R, n1 + 0.5, n2 + 0.5)
    return min(x0, x1), max(x0, x1), min(y0, y1), max(y0, y1)
end

function Base.extrema(R::Region)
    x0, y0 = grid2world(R, 1, 1)
    x1, y1 = grid2world(R, size(R))
    return min(x0, x1), max(x0, x1), min(y0, y1), max(y0, y1)
end

end # module
