const default_antialiasing = 11

"""
    FourierOptics.insert_mask(Fin, args...; kwds...) -> Fout

multiplies the complex amplitude of the input field `Fin` by a mask defined by
arguments `args...` and keywords `kwds...` and returns the resulting output
field `Fout`. The input field `Fin` is left unmodified, method
[`FourierOptics.insert_mask!`](@ref) may be used for in-place operation.
See [`FourierOptics.forge_mask`](@ref) for how to define a mask.

"""
insert_mask(F::Field, args...; kwds...) = insert_mask!(copy(F), args...; kwds...)

"""
    FourierOptics.insert_mask!(F, args...) -> F

multiplies in-place the complex amplitude of the field `F` by a mask defined by
arguments `args...` and keywords `kwds...` and returns `F`. See
[`FourierOptics.insert_mask`](@ref) for an out-of-place version and for
details.

"""
insert_mask!(F::Field, args...; kwds...) = multiply!(F, forge_mask(F, args...; kwds...))

# Alias for a 2-tuple representing a pair of coordinates.
const TwoTuple{A,B} = Tuple{A,B}

Base.eltype(x::GeometricalObject) = Base.eltype(typeof(x))
Base.eltype(::Type{<:GeometricalObject{T}}) where {T} = T

"""
    FourierOptics.Point{T}(x, y)
    FourierOptics.Point{T}((x, y))
    FourierOptics.Point{T}(; x=..., y=...)
    FourierOptics.Point{T}(; r=..., θ=...)

construct a point given its Cartesian coordinates `(x,y)` (in the 3 first
examples) or its polar coordinates (in the last example) with `r` the distance
to the origin and `θ` the counterclockwise angle relative to `x`-axis.
Parameter `T` is the type used to store coordinates, it may be omitted.

Note that, in `FourierOptics`, `x` and `y` respectively correspond to the 1st
and 2nd dimensions of 2-dimensional arrays.

"""
Point{T}(xy::TwoTuple) where {T} = Point{T}(xy...)
Point(xy::TwoTuple) = Point(xy...)
Point(x, y) = Point{promote_type(typeof(x), typeof(y))}(x, y)
@inline Point(; kwds...) = build(Point; kwds...)
@inline Point{T}(; kwds...) where {T} = build(Point{T}; kwds...)
@inline build(::Type{T}; x=nothing, y=nothing, r=nothing, θ=nothing) where {T<:Point} =
    if (x !== nothing) & (y !== nothing) & (r === nothing) & (θ === nothing)
        return T(x, y)
    elseif (x === nothing) & (y === nothing) & (r !== nothing) & (θ !== nothing)
        sinθ, cosθ = sincos(θ)
        return T(r*cosθ, r*sinθ)
    else
        throw(ArgumentError(
            "exclusively keywords `x` and `y` or `r` and `θ` must be specified"))
    end

"""
    FourierOptics.Box{T}(ll::Point, ur::Point)
    FourierOptics.Box{T}((llx,lly), (urx,ury))

construct a rectangular box with edges aligned with the Cartesian axes and
given the coordinates of its lower-left, `ll`, and upper-right, `ur`, corners.
Parameter `T` is the type used to store coordinates, it may be omitted.

Corner coordinates may be specified as points or as 2-tuples: `(llx,lly) =
(ll.x,ll.y)` `(urx,ury) = (ur.x,ur.y)`.

If `llx > urx` or `lly > ury`, the box is considered as **empty**.

Boxes are used to represent grid cells and bounding-boxes of other geometrical
shape. Use [`FourierOptics.Rectangle`](@ref) if you want to define rectangular
masks.

"""
Box(ll::Point, ur::Point) = Box{promote_type(eltype(ll),eltype(ur))}(ll, ur)
Box(ll::TwoTuple, ur::TwoTuple) = Box(Point(ll), Point(ur))
Box{T}(ll::TwoTuple, ur::TwoTuple) where {T} = Box{T}(Point(ll), Point(ur))
@inline Box(; kwds...) = build(Box; kwds...)
@inline Box{T}(; kwds...) where {T} = build(Box{T}; kwds...)
@inline build(::Type{T}; xmin, xmax, ymin, ymax) where{T<:Box} = T((xmin,ymin), (xmax,ymax))

"""
    FourierOptics.Rectangle{T}(p0::Point, p1::Point)
    FourierOptics.Rectangle{T}((x0,y0), (x1,y1))

construct a rectangle with edges aligned with the Cartesian axes and given the
coordinates of 2 opposite corners. Parameter `T` is the type used to store
coordinates, it may be omitted.

Corner coordinates may be specified as points or as 2-tuples: `(x0,y0) =
(p0.x,p0.y)` `(x1,y1) = (p1.x,p1.y)`.

Compared to a [`FourierOptics.Box`](@ref), a rectangle is never empty (it
contains at least a point).

"""
Rectangle(p0::Point, p1::Point) = Rectangle{promote_type(eltype(p0),eltype(p1))}(p0, p1)
Rectangle(p0::TwoTuple, p1::TwoTuple) = Rectangle(Point(p0), Point(p1))
Rectangle{T}(p0::TwoTuple, p1::TwoTuple) where {T} = Rectangle{T}(Point(p0), Point(p1))
@inline Rectangle(; kwds...) = build(Rectangle; kwds...)
@inline Rectangle{T}(; kwds...) where {T} = build(Rectangle{T}; kwds...)
@inline build(::Type{T}; x0, y0, x1, y1) where{T<:Rectangle} = T((x0,y0), (x1,y1))

"""
    FourierOptics.Circle{T}(center::Point, radius)
    FourierOptics.Circle{T}((x, y), radius)
    FourierOptics.Circle{T}(; center=..., radius=...)
    FourierOptics.Circle{T}(; center=..., diameter=...)

construct a circle of given `center` and `radius` or `diameter`. The center may
be specified by its coordinates `(x,y)` along the Cartesian axes. Parameter `T`
is the type used to store coordinates, it may be omitted.

"""
Circle{T}(center::TwoTuple, radius) where {T} = Circle{T}(Point(center), radius)
Circle(center::TwoTuple, radius) = Circle(Point(center), radius)
Circle(center::Point, radius) =
    Circle{promote_type(eltype(center),typeof(radius))}(center, radius)
@inline Circle(; kwds...) = build(Circle; kwds...)
@inline Circle{T}(; kwds...) where {T} = build(Circle{T}; kwds...)
@inline build(::Type{T}; center, radius=nothing, diameter=nothing) where{T<:Circle} =
    if (radius !== nothing) & (diameter === nothing)
        return T(center, radius)
    elseif (radius === nothing) & (diameter !== nothing)
        return T(center, diameter/2)
    else
        throw(ArgumentError(
            "exclusively keywords `radius` or `diameter` must be specified"))
    end

Polygon() = throw(ArgumentError("polygon must have at least 3 vertices"))
Polygon(points::TwoTuple...) = Polygon(points)
Polygon(points::Point...) = Polygon(points)
Polygon(points::Tuple{Vararg{Point{T}}}) where {T} = Polygon{T}(points)
Polygon(points::Tuple{Vararg{TwoTuple}}) = Polygon(map(Point, points))
Polygon(points::Tuple{Vararg{Point}}) = Polygon{promote_eltype(points...)}(points)
Polygon{T}() where {T} = throw(ArgumentError("polygon must have at least 3 vertices"))
Polygon{T}(points::TwoTuple...) where {T} = Polygon{T}(points)
Polygon{T}(points::Point...) where {T} = Polygon{T}(points)
Polygon{T}(points::Tuple{Vararg{TwoTuple}}) where {T} = Polygon{T}(map(Point{T}, points))
Polygon{T}(points::Tuple{Vararg{Point}}) where {T} = Polygon(map(Point{T}, points))

# Type conversion and copy constructors for geometrical objects.
for type in (:GeometricalObject, :Point, :Box,
             :ShapeObject, :Rectangle, :Circle, :Polygon,
             :MaskObject, :OpaqueMask, :TransparentMask)
    @eval begin
        $type(obj::$type) = obj
        $type{T}(obj::$type{T}) where {T} = obj
        $type{T}(obj::$type) where {T} = convert_eltype(T, obj)
    end
    if type ∈ (:Point, :Box, :Rectangle, :Circle, :Polygon)
        @eval begin
            convert_eltype(::Type{T}, obj::$type{T}) where {T} = obj
            convert_eltype(::Type{T}, obj::$type) where {T} = $type{T}(parts(obj)...)
        end
    elseif type ∈ (:OpaqueMask, :TransparentMask)
        @eval begin
            $type{T}(obj::ShapeObject) where {T} = $type(convert_eltype(T, obj))
            convert_eltype(::Type{T}, obj::$type{T}) where {T} = obj
            convert_eltype(::Type{T}, obj::$type) where {T} =
                $type(convert_eltype(T, shape(obj)))
        end
    end
end
Base.convert(::Type{T}, obj::T) where {T<:GeometricalObject} = obj
Base.convert(::Type{T}, obj) where {T<:GeometricalObject} = T(obj)

function Base.collect(obj::MaskObject, objs::MaskObject...)
    T = promote_eltype(obj, objs...)
    return MaskObject{T}[obj, objs...]
end

Base.Tuple(obj::Union{Point,Box,Rectangle,Circle}) = parts(obj)

# Iterating over a basic geometrical object is like iteration over the
# arguments of the most basic constructor.
@inline Base.iterate(obj::Point, i::Int=1) =
    i == 1 ? (obj.x, 2) :
    i == 2 ? (obj.y, 3) : nothing
@inline Base.iterate(obj::RectangularObject, i::Int=1) =
    i == 1 ? (first(obj), 2) :
    i == 2 ? (last( obj), 3) : nothing
@inline Base.iterate(obj::Circle, i::Int=1) =
    i == 1 ? (center(obj), 2) :
    i == 2 ? (radius(obj), 3) : nothing

Base.show(io::IO, obj::Point) =
    print(io, "Point(", obj.x, ", ", obj.y, ')')
Base.show(io::IO, obj::Box) =
    print(io, "Box((", obj.xmin, ", ", obj.ymin, "), (", obj.xmax, ", ", obj.ymax, "))")
function Base.show(io::IO, obj::RectangularObject)
    (x0, y0), (x1, y1) = obj
    print(io, showntype(obj), "((", x0, ", ", y0, "), (", x1, ", ", y1, "))")
end
function Base.show(io::IO, obj::CircularObject)
    c = center(obj)
    r = radius(obj)
    print(io, showntype(obj), "(center = (", c.x, ", ", c.y, "), radius = ", r, ')')
end
function Base.show(io::IO, obj::PolygonalObject)
    print(io, showntype(obj), "(")
    flag = false
    for pnt in vertices(obj)
        print(io, (flag ? "), (" : "("), pnt.x, ", ", pnt.y)
        flag = true
    end
    print(io, (flag ? "))" : ")"))
end

Base.:(==)(a::Point, b::Point) =
    (a.x == b.x) &
    (a.y == b.y)
Base.:(==)(a::Box, b::Box) =
    (a.xmin == b.xmin) &
    (a.ymin == b.ymin) &
    (a.xmax == b.xmax) &
    (a.ymax == b.ymax)
Base.:(==)(a::Rectangle, b::Rectangle) =
    (a.x0 == b.x0) &
    (a.y0 == b.y0) &
    (a.x1 == b.x1) &
    (a.y1 == b.y1)
Base.:(==)(a::Circle, b::Circle) =
    (a.x == b.x) &
    (a.y == b.y) &
    (a.r == b.r)
function Base.:(==)(a::Polygon, b::Polygon)
    va, vb = vertices(a), vertices(b)
    len = length(va)
    length(vertices(b)) == len || return false
    @inbounds for i in 1:len
        va[i] == vb[i] || return false
    end
    return true
end

Base.isapprox(a::Point, b::Point; kwds...) =
    isapprox(a.x, b.x; kwds...) &&
    isapprox(a.y, b.y; kwds...)
Base.isapprox(a::Box, b::Box; kwds...) =
    isapprox(a.xmin, b.xmin; kwds...) &&
    isapprox(a.ymin, b.ymin; kwds...) &&
    isapprox(a.xmax, b.xmax; kwds...) &&
    isapprox(a.ymax, b.ymax; kwds...)
Base.isapprox(a::Rectangle, b::Rectangle; kwds...) =
    isapprox(a.x0, b.x0; kwds...) &&
    isapprox(a.y0, b.y0; kwds...) &&
    isapprox(a.x1, b.x1; kwds...) &&
    isapprox(a.y1, b.y1; kwds...)
Base.isapprox(a::Circle, b::Circle; kwds...) =
    isapprox(a.x, b.x; kwds...) &&
    isapprox(a.y, b.y; kwds...) &&
    isapprox(a.r, b.r; kwds...)
function Base.isapprox(a::Polygon, b::Polygon; kwds...)
    va, vb = vertices(a), vertices(b)
    len = length(va)
    length(vertices(b)) == len || return false
    @inbounds for i in 1:len
        isapprox(va[i], vb[i]; kwds...) || return false
    end
    return true
end

for type in (:Point, :Box, :Rectangle, :Circle, :Polygon,
             :RectangularAperture, :RectangularObscuration,
             :CircularAperture, :CircularObscuration,
             :PolygonalAperture, :PolygonalObscuration)
    @eval begin
        showntype(::Union{$type,Type{<:$type}}) = $(string(type))
    end
end

"""
    FourierOptics.grow(box, dx, dy=dx)

yields a new box object corresponding to the input `box` object with 1st and
2nd dimensions respectively grown by `dx` and `dy`.

Note that the algebraic (not absolute) values are applied. Hence, if `dx` and
`dy` are both negative, the box is effectively shrunk by `abs(dx)` and
`abs(dy)`.

See also [`FourierOptics.shrink`](@ref).

"""
grow(box::Box, dx, dy=dx) = Box((box.xmin - dx, box.ymin - dy),
                                (box.xmax + dx, box.ymax + dy))

"""
    FourierOptics.shrink(box, dx, dy=dx)

yields a new box object corresponding to the input `box` object with 1st and
2nd dimensions respectively shrunk by `dx` and `dy`.

Note that the algebraic (not absolute) values are applied. Hence, if `dx` and
`dy` are both negative, the box is effectively grown by `abs(dx)` and
`abs(dy)`.

See also [`FourierOptics.grow`](@ref).

"""
shrink(box::Box, dx, dy=dx) = Box((box.xmin + dx, box.ymin + dy),
                                  (box.xmax - dx, box.ymax - dy))

Base.intersect(a::Box, b::Box) = Box((fastmax(a.xmin, b.xmin), fastmax(a.ymin, b.ymin)),
                                     (fastmin(a.xmax, b.xmax), fastmin(a.ymax, b.ymax)))
Base.isempty(box::Box) = (box.xmin > box.xmax)|(box.ymin > box.ymax)
Base.first(obj::Box) = Point(obj.xmin, obj.ymin)
Base.last(obj::Box) = Point(obj.xmax, obj.ymax)

Base.first(obj::Rectangle) = Point(obj.x0, obj.y0)
Base.last(obj::Rectangle) = Point(obj.x1, obj.y1)

Base.first(obj::RectangularMask) = first(shape(obj))
Base.last(obj::RectangularMask) = last(shape(obj))

# Methods for circular objects.
diameter(obj::CircularObject) = twice(radius(obj))
radius(obj::Circle) = obj.r
radius(obj::CircularMask) = radius(shape(obj))
center(obj::Circle) = Point(obj.x, obj.y)
center(obj::CircularMask) = center(shape(obj))

twice(x) = x + x

# Methods for polygonal objects.
vertices(obj::Polygon) = obj.vertices
vertices(obj::PolygonalMask) = vertices(shape(obj))
function vertices(obj::Union{Box,Rectangle})
    (x0, y0), (x1, y1) = obj
    return Point(x0,y0), Point(x1,y0), Point(x1,y1), Point(x0,y1)
end

shape(obj::MaskObject) = obj.shape
shape(obj::ShapeObject) = obj

is_opaque(obj::OpaqueMask) = true
is_opaque(obj::MaskObject) = false

is_transparent(obj::MaskObject) = !is_opaque(obj)

is_convex(mask::MaskObject) = is_convex(shape(mask))
is_convex(rect::Rectangle) = true
is_convex(circle::Circle) = true
is_convex(poly::Polygon) = poly.convex

# Geometrical objects can be used to do some math...
Base.abs2(obj::Point) = abs2(obj.x) + abs2(obj.y)
Base.abs(obj::Point) = sqrt(abs2(obj))
Base.Math.atan(obj::Point) = atan(obj.y, obj.x)
LinearAlgebra.norm(obj::Point) = abs(obj)
inner(a::Point, b::Point) = a.x*b.x + a.y*b.y # scalar/inner product
outer(a::Point, b::Point) = a.x*b.y - a.y*b.x # outer product
Base.:(*)(a::Point, b::Point) = inner(a, b) # scalar product

# Just permute operands for some common operations.
Base.:(+)(a::Point, b::GeometricalObject) = b + a
Base.:(*)(a::GeometricalObject, b::Number) = b*a
Base.:(\)(a::Number, b::GeometricalObject) = b/a

# Adding a point to a shape amounts to shifting the shape. Unary plus does
# nothing. Unary minus behaves as multiplying by -1. Multiplying/dividing a
# geometrical object by a number amounts to scaling around origin.
Base.:(+)(a::GeometricalObject) = a

Base.:(-)(a::Point) = Point(-a.x, -a.y)
Base.:(+)(a::Point,  b::Point) = Point(a.x + b.x, a.y + b.y)
Base.:(-)(a::Point,  b::Point) = Point(a.x - b.x, a.y - b.y)
Base.:(*)(a::Number, b::Point) = Point(a*b.x, a*b.y)
Base.:(/)(a::Point,  b::Number) = Point(a.x/b, a.y/b)

Base.:(-)(a::Circle) = Circle(-center(a), radius(a))
Base.:(+)(a::Circle, b::Point) = Circle(center(a) + b, radius(a))
Base.:(-)(a::Circle, b::Point) = Circle(center(a) - b, radius(a))
Base.:(-)(a::Point,  b::Circle) = Circle(a - center(b), radius(b))
Base.:(*)(a::Number, b::Circle) = Circle(a*center(b), abs(a)*radius(b))
Base.:(/)(a::Circle, b::Number) = Circle(center(a)/b, radius(a)/abs(b))

Base.:(-)(a::Polygon) = Polygon(map(-, vertices(a)))
Base.:(+)(a::Polygon, b::Point) = Polygon(map(Fix2(+, b), vertices(a)))
Base.:(-)(a::Polygon, b::Point) = Polygon(map(Fix2(-, b), vertices(a)))
Base.:(-)(a::Point,   b::Polygon) = Polygon(map(Fix1(-, a), vertices(b)))
Base.:(*)(a::Number,  b::Polygon) = Polygon(map(Fix1(*, a), vertices(b)))
Base.:(/)(a::Polygon, b::Number) = Polygon(map(Fix2(/, b), vertices(a)))

for type in (:Box, :Rectangle)
    @eval begin
        Base.:(-)(a::$type) = $type(-first(a), -last(a))
        Base.:(+)(a::$type,  b::Point) = $type(first(a) + b, last(a) + b)
        Base.:(-)(a::$type,  b::Point) = $type(first(a) - b, last(a) - b)
        Base.:(-)(a::Point,  b::$type) = $type(a - first(b), a - last(b))
        Base.:(*)(a::Number, b::$type) = $type(a*first(b), a*last(b))
        Base.:(/)(a::$type,  b::Number) = $type(first(a)/b, last(a)/b)
    end
end

# Performing a geometrical operation on a mask amounts to performing the
# operation on the embedded shape.
for type in (:OpaqueMask, :TransparentMask)
    @eval begin
        Base.:(-)(a::$type) = $type(-shape(a))
        Base.:(+)(a::$type,  b::Point) = $type(shape(a) + b)
        Base.:(-)(a::$type,  b::Point) = $type(shape(a) - b)
        Base.:(-)(a::Point,  b::$type) = $type(a - shape(b))
        Base.:(*)(a::Number, b::$type) = $type(a*shape(b))
        Base.:(/)(a::$type,  b::Number) = $type(shape(a)/b)
    end
end

"""
    FourierOptics.Box{T=eltype(obj)}(obj)

yields the bounding-box of the geometrical object `obj`.

"""
Box(obj::ShapeObject) = Box{eltype(obj)}(obj)
Box{T}(obj::MaskObject) where {T} = Box{T}(shape(obj))
Box{T}(obj::Rectangle) where {T} = Box{T}(first(obj), last(obj))
function Box{T}(obj::Circle) where {T}
    c = center(obj)
    r = radius(obj)
    return Box{T}((c.x - r, c.y - r), (c.x + r, c.y + r))
end
Box{T}(obj::Polygon) where {T} = convert_eltype(T, obj.box)

# `parts(obj)` yields the individual elements of `obj` from which it can be
# re-built without ambiguities.
parts(obj::Point) = (obj.x, obj.y)
parts(obj::Union{Box,Rectangle}) = (first(obj), last(obj))
parts(obj::Circle) = (center(obj), radius(obj))
parts(obj::Polygon) = vertices(obj)
parts(obj::MaskObject) = parts(shape(obj))

# Fast versions of `min` and `max` which return their first argument if any
# argument is a NaN. The fast version of `minmax` cannot warrants this
# property.
fastmin(x::T, y::T) where {T} = (x > y ? y : x)
fastmax(x::T, y::T) where {T} = (x < y ? y : x)
fastminmax(x::T, y::T) where {T} = (x > y ? (y,x) : (x,y))
for func in (:fastmax, :fasmin, :fastminmax)
    @eval $func(x, y) = $func(promote(x, y)...)
end

"""
    FourierOptics.RectangularAperture(args...; kwds...)

yields a mask object representing a rectangular aperture defined by given
arguments `args...` and keywords `kwds...` and whose edges are aligned with the
Cartesian axes. See [`FourierOptics.Rectangle`](@ref) constructor for possible
arguments and keywords. A rectangular aperture is a transparent rectangular
mask.

"""
RectangularAperture(args...; kwds...) = TransparentMask(Rectangle(args...; kwds...))

"""
    FourierOptics.RectangularObscuration(args...; kwds...)

yields a mask object representing a rectangular obscuration defined by given
arguments `args...` and keywords `kwds...` and whose edges are aligned with the
Cartesian axes. See [`FourierOptics.Rectangle`](@ref) constructor for possible
arguments and keywords. A rectangular obscuration is an opaque rectangular
mask.

"""
RectangularObscuration(args...; kwds...) = OpaqueMask(Rectangle(args...; kwds...))

"""
    FourierOptics.CircularAperture(args...; kwds...)

yields a mask object representing a circular aperture defined by given
arguments `args...` and keywords `kwds...`. See [`FourierOptics.Circle`](@ref)
constructor for possible arguments and keywords. A circular aperture is a
transparent circular mask.

"""
CircularAperture(args...; kwds...) = TransparentMask(Circle(args...; kwds...))

"""
    FourierOptics.CircularObscuration(args...; kwds...)

yields a mask object representing a circular obscuration defined by given
arguments `args...` and keywords `kwds...`. See [`FourierOptics.Circle`](@ref)
constructor for possible arguments and keywords. A circular obscuration is an
opaque circular mask.

"""
CircularObscuration(args...; kwds...) = OpaqueMask(Circle(args...; kwds...))

"""
    FourierOptics.PolygonalAperture(args...; kwds...)

yields a mask object representing a polygonal aperture defined by given
arguments `args...` and keywords `kwds...`. See [`FourierOptics.Polygon`](@ref)
constructor for possible arguments and keywords. A polygonal aperture is a
transparent polygonal mask.

"""
PolygonalAperture(args...; kwds...) = TransparentMask(Polygon(args...; kwds...))

"""
    FourierOptics.PolygonalObscuration(args...; kwds...)

yields a mask object representing a polygonal obscurationdefined by given
arguments `args...` and keywords `kwds...`. See [`FourierOptics.Polygon`](@ref)
constructor for possible arguments and keywords. A polygonal obscuration is an
opaque polygonal mask.

"""
PolygonalObscuration(args...; kwds...) = OpaqueMask(Polygon(args...; kwds...))

# Yields a range with `n` points centered around zero and step `s/n`. This
# method is used to index a sub-cell.
function subrange(n::Int, s::Number)
    r = (1:n) .* (s/n)
    return r .- (first(r) + last(r))/2
end

"""
    FourierOptics.forge_mask(F::Field, objs...; kwds...) -> msk

yields a 2-dimensional transmission mask for the field `F` and combining
elementary mask objects `objs...`.

"""
forge_mask(F::Field, args...; kwds...) =
    forge_mask!(Matrix{floating_point_type(F)}(undef, size(F)), F, args...; kwds...)

function forge_mask!(dst::AbstractMatrix, F::Field, args...; kwds...)
    X = grid_step(F)*RolledCoordinates(grid_size(F))
    Y = X
    return forge_mask!(dst, X, Y, args...; kwds...)
end

"""
    FourierOptics.forge_mask(X, Y, objs...; kwds...) -> msk

yields a 2-dimensional transmission mask with coordinates given by `X` and `Y`
along the 1st and 2nd dimensions and combining aperture/obscuration objects
`objs...`. The following *painting* algorithm is used:

- The mask is initially filled with the transparent or opaque value depending
  on whether the first component is opaque or transparent.

- Then, for each component in turn, the cells of the mask that are inside the
  component are painted with the opaque or transparent value depending on
  whether the component is opaque or transparent.

- The parts of the mask overlapping the boundaries of the topmost components
  are set to an intermediate value between the opaque and transparent ones and
  (approximately) proportionally to the transparent fraction of the cell area.

Note that the order of the components of the mask is relevant: an aperture
component drills holes in the previously opaque parts while an obscuration
masks previously transparent parts.

Keyword `antialiasing` can be set to specify the number of sub-cells (per side)
to determine the transmission of grid cells partially overlapping the boundary
delimiting the mask components. By default, `antialiasing =
$default_antialiasing`. If `antialiasing ≤ 1`, a 50% transmission is assumed
for partially overlapping cells (sharp edges); otherwise, overlapping cells are
subdivided in `antialiasing × antialiasing` sub-cells to estimate their partial
transmission.

Keywords `opaque` and `transparent` can be used to specify the values of the
the respectively opaque and transparent parts of the mask. Values of partially
opaque/transparent parts will be interpolated between these.

Example:

    using FourierOptics
    using FourierOptics: Point
    δ = 0.07mm # grid sampling step
    N = 512   # grid size
    X = RolledCoordinates(N)*δ # coordinates
    Y = X
    center = Point(0mm,0mm)
    outer_radius = 3.2mm
    inner_radius = outer_radius/3
    width = 0.23mm # width of spider arms
    height = 2*(outer_radius + δ) # length of spider arms
    voff = Point(width, height)
    hoff = Point(height, width)
    mask = forge_mask(
        X, Y,
        CircularAperture(center, outer_radius), # aperture
        CircularObscuration(center, inner_radius), # central obscuration
        RectangularObscuration(center - voff/2, center + voff/2),
        RectangularObscuration(center - hoff/2, center + hoff/2),
        ...)

"""
function forge_mask(X::AbstractVector,
                    Y::AbstractVector,
                    args...; kwds...)
    T = floating_point_type(promote_eltype(X, Y))
    return forge_mask!(Matrix{T}(undef, length(X), length(Y)), X, Y, args...; kwds...)
end

function forge_mask!(dst::AbstractMatrix,
                     X::AbstractVector,
                     Y::AbstractVector,
                     args::Tuple{Vararg{MaskObject}};
                     kwds...)
    return forge_mask!(dst, X, Y, args...; kwds...)
end

function forge_mask!(dst::AbstractMatrix,
                     X::AbstractVector,
                     Y::AbstractVector,
                     args::MaskObject...;
                     antialiasing::Integer = 11,
                     opaque = zero(eltype(dst)),
                     transparent = oneunit(eltype(dst)))
    # Check that coordinate units are compatible and determine a suitable unit
    # for all coordinates.
    T = float(promote_eltype(X, Y, args...))

    # Check arguments.
    I, J = axes(dst)
    axes(X) == (I,) || throw(DimensionMismatch("coordinates `X` have incompatible indices"))
    axes(Y) == (J,) || throw(DimensionMismatch("coordinates `Y` have incompatible indices"))
    antialiasing = as(Int , antialiasing)
    opaque = as(eltype(dst), opaque)
    transparent = as(eltype(dst), transparent)

    # Call the real method.
    return unsafe_forge_mask!(dst, convert_eltype(T, X), convert_eltype(T, Y),
                              map(A -> convert_eltype(T, A), args)...;
                              antialiasing, opaque, transparent)
end

function unsafe_forge_mask!(dst::AbstractMatrix{T},
                            X::AbstractVector{C},
                            Y::AbstractVector{C},
                            objs::MaskObject{C}...;
                            antialiasing::Int,
                            opaque::T,
                            transparent::T) where {T,C}
    partial = interpolate(opaque, transparent, 1//2)
    δx = grid_step(X)
    δy = grid_step(Y)

    # Determine whether a cell of the mask is fully or partially opaque of
    # transparent.
    initial = true
    for obj in objs
        unsafe_forge_mask!(dst, X, δx, Y, δy, obj, opaque, partial, transparent, initial)
        initial = false
    end
    if initial
        # No objects, assume mask is transparent.
        fill!(mask, transparent)
    end

    if antialiasing > 1
        # Refine mask at partially opaque cells.
        I, J = axes(dst)
        DX = subrange(antialiasing, δx)
        DY = subrange(antialiasing, δy)
        state = Array{Bool}(undef, antialiasing, antialiasing)
        @inbounds for j in J
            y = Y[j]
            for i in I
                if dst[i,j] == partial
                    # Cell is partially opaque/transparent.
                    x = X[i]
                    count = -1 # to trigger resetting of state array
                    for obj in objs
                        count = unsafe_forge_mask!(state, x, DX, y, DY, obj, count)
                    end
                    #@assert max(count, 0) == Base.count(state)
                    fraction = max(count, 0)//antialiasing^2
                    dst[i,j] = interpolate(opaque, transparent, fraction)
                end
            end
        end
    end

    return dst
end

function unsafe_forge_mask!(dst::AbstractMatrix{V},
                            X::AbstractVector{T}, δx::T,
                            Y::AbstractVector{T}, δy::T,
                            mask::MaskObject{T},
                            opaque::V,
                            partial::V,
                            transparent::V,
                            initial::Bool) where {T,V}
    if initial
        fill!(dst, is_opaque(mask) ? transparent : opaque)
    end
    return unsafe_forge_mask!(dst, X, δx, Y, δy, shape(mask),
                              is_opaque(mask) ? opaque : transparent,
                              partial)
end

function unsafe_forge_mask!(dst::AbstractMatrix{V},
                            X::AbstractVector{T}, δx::T,
                            Y::AbstractVector{T}, δy::T,
                            obj::ShapeObject{T},
                            value::V,
                            partial::V) where {T,V}
    box = grow(Box(obj), δx, δy)
    I, J = axes(dst)
    @inbounds for j in J
        y = Y[j]
        if (y < box.ymin)|(y > box.ymax)
            # Cell is outside the bounding-box of the mask.
            continue
        end
        for i in I
            x = X[i]
            if (x < box.xmin)|(x > box.xmax)
                # Cell is outside the bounding-box of the mask.
                continue
            end
            cell = Box((x - δx/2, y - δy/2), (x + δx/2, y + δy/2))
            overlap = Overlap(cell, obj)
            if overlap == INSIDE
                dst[i,j] = value
            elseif overlap == PARTIAL
                dst[i,j] = partial
            end
        end
    end
    return dst
end

function unsafe_forge_mask!(state::AbstractMatrix{Bool},
                            x0::T, dx::AbstractVector{T},
                            y0::T, dy::AbstractVector{T},
                            obj::MaskObject{T},
                            count::Int) where {T}
    # FIXME: Skip all this if pixel is outside boundaries.
    I, J = axes(state)
    if count < 0
        # This is the first mask applied to this cell.
        count = 0
        @inbounds for j in J
            y = y0 + dy[j]
            for i in I
                x = x0 + dx[i]
                inside = Overlap(Point(x,y), obj) != OUTSIDE
                transp = xor(inside, is_opaque(obj))
                state[i,j] = transp
                count += transp
            end
        end
    else
        @inbounds for j in J
            y = y0 + dy[j]
            for i in I
                x = x0 + dx[i]
                if Overlap(Point(x,y), obj) != OUTSIDE
                    # Point is considered as being inside boundaries of mask.
                    if is_opaque(obj) # NOTE: This branch optimized out by compiler.
                        # Sub-cell must be opaque.
                        if state[i,j]
                            # Sub-cell previously counted as transparent.
                            state[i,j] = false
                            count -= 1
                        end
                    else
                        # Sub-cell must be transparent.
                        if !state[i,j]
                            # Sub-cell previously considered as opaque.
                            state[i,j] = true
                            count += 1
                        end
                    end
                end
            end
        end
    end
    return count
end

"""
    FourierOptics.grid_step(x::AbstractVector)

yields the step along a vector of coordinates throwing an error if the
increment between successive values of `x` is not positive or not uniform.

"""
function grid_step(rng::Union{AbstractRange,AbstractCoordinates})
    stp = step(rng)
    stp > zero(stp) || throw(ArgumentError("grid step must be positive"))
    return stp
end

function grid_step(x::AbstractVector)
    len = length(x)
    len > 1 || throw(ArgumentError("zero-length vector of coordinates has no defined step"))
    stp = (last(x) - first(x))/(len - 1)
    stp > zero(stp) || throw(ArgumentError("grid step must be positive"))
    tol = sqrt(eps(floating_point_type(stp)))
    stpmin = stp*(one(tol) - tol)
    stpmax = stp*(one(tol) + tol)
    i = firstindex(x)
    i_last = lastindex(x)
    flag = true
    @inbounds while i < i_last
        i += 1
        flag &= (stpmin ≤ x[i] - x[i-1] ≤ stpmax)
    end
    flag || throw(ArgumentError("vector of coordinates has non-uniform steps"))
    return stp
end

"""
    FourierOptics.Overlap(pxl, obj)

yields the overlapping of pixel `pxl` with shape object `obj`. `pxl` may be a
point or a box. Returned value is:

- `INSIDE` if `pxl` is fully inside the boundaries of `obj`;

- `OUTSIDE` if `pxl` is fully outside the boundaries of `obj`;

- `PARTIAL` if `pxl` straddles some boundaries of `obj`.

"""
Overlap(pxl, obj::MaskObject) = Overlap(pxl, shape(obj))

function Overlap(cell::Box{T}, rect::Rectangle{T}) where {T}
    (x0, y0), (x1, y1) = cell
    (xmin, ymin), (xmax, ymax) = rect
    if (x0 ≥ xmin) & (x1 ≤ xmax) & (y0 ≥ ymin) & (y1 ≤ ymax)
        return INSIDE
    elseif (x0 > xmax) | (x1 < xmin) | (y0 > ymax) | (y1 < ymin)
        return OUTSIDE
    else
        return PARTIAL
    end
end

function Overlap(point::Point{T}, rect::Rectangle{T}) where {T}
    x, y = point
    (xmin, ymin), (xmax, ymax) = rect
    if (xmin < x < xmax) & (ymin < y < ymax)
        return INSIDE
    elseif (x < xmin) | (x > xmax) | (y < ymin) | (y > ymax)
        return OUTSIDE
    else
        return PARTIAL
    end
end

function Overlap(cell::Box{T}, circle::Circle{T}) where {T}
    # Get coordinates of cell corners relative to circle center and quickly
    # check whether cell is certainly outside circle.
    x0 = cell.xmin - circle.x
    x1 = cell.xmax - circle.x
    y0 = cell.ymin - circle.y
    y1 = cell.ymax - circle.y
    r = radius(circle)
    if (x0 > r) | (x1 < -r) | (y0 > r) | (y1 < -r)
        return OUTSIDE
    end
    r² = r^2
    # Most distant point in cell from center of circle is one of the corners.
    # If this point is inside the circle, then the cell is fully inside the
    # circle.
    x²max = fastmax(x0^2, x1^2)
    y²max = fastmax(y0^2, y1^2)
    if x²max + y²max ≤ r²
        return INSIDE
    end
    # Otherwise, if closest point of one of the cell edge is inside the circle,
    # overlap is partial.
    if clamp(zero(T), x0, x1)^2 + y²max ≤ r²
        return PARTIAL
    end
    if x²max + clamp(zero(T), y0, y1)^2 ≤ r²
        return PARTIAL
    end
    # Otherwise, there is no overlap.
    return OUTSIDE
end

function Overlap(point::Point{T}, circle::Circle{T}) where {T}
    d² = abs2(point - center(circle))
    r² = abs2(radius(circle))
    if d² < r²
        return INSIDE
    elseif d² == r²
        return PARTIAL
    else
        return OUTSIDE
    end
end

function Overlap(cell::Box{T}, polygon::Polygon{T,N}) where {T,N}
    (x0, y0), (x1, y1) = cell # retrieve coordinates of cell vertices
    V = vertices(polygon) # retrieve the list of vertices of the polygon

    # If any polygon edges (strictly) crosses one of the 4 edges of the cell, then
    # there is some partial overlapping.
    A = last(V) # initialize A, the 1st point of edges in cyclic list of points
    @inbounds for B in V # loop over B, the 2nd point of edges
        xmin, xmax = fastminmax(A.x, B.x)
        ymin, ymax = fastminmax(A.y, B.y)
        if xmin < x0 < xmax
            # Left edge of cell may cross the segment AB.
            y = interpolate(A.y, B.y, (x0 - A.x)/(B.x - A.x))
            if y0 < y < y1
                return PARTIAL
            end
        end
        if xmin < x1 < xmax
            # Right edge of cell may cross the segment AB.
            y = interpolate(A.y, B.y, (x1 - A.x)/(B.x - A.x))
            if y0 < y < y1
                return PARTIAL
            end
        end
        if ymin < y0 < ymax
            # Bottom edge of cell may cross the segment AB.
            x = interpolate(A.x, B.x, (y0 - A.y)/(B.y - A.y))
            if x0 < x < x1
                return PARTIAL
            end
        end
        if ymin < y1 < ymax
            # Top edge of cell may cross the segment AB.
            x = interpolate(A.x, B.x, (y1 - A.y)/(B.y - A.y))
            if x0 < x < x1
                return PARTIAL
            end
        end
        A = B # update 1st point of next edge
    end

    # Otherwise, the cell is either fully inside or fully outside the polygon.
    # This can be tested for any point of the cell, its center for example.
    return Overlap(Point((x0 + x1)/2, (y0 + y1)/2), polygon)
end

Overlap(point::Point{T}, polygon::Polygon{T,N}) where {T,N} =
    winding_number_test(point, vertices(polygon)) ? INSIDE : OUTSIDE
