#
# coords.jl -
#
# Implement coordinate transforms.
#

module CoordinateTransforms

export
    CoordinateTransform,
    compose,
    jacobian,
    origin

"""
```julia
CoordinateTransform(stp, x0, y0) -> C
```

yields a simple 2D coordinate transform `C` such that:

```julia
C(u, v) -> (x0 + stp*u, y0 + stp*v)
```

The parameters of the coordinate transform `C` can be retrieved as follows:

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

Coordinate transforms may be inverted, pre/post-scaled, pre-post-translated:

```julia
A = inv(C)        # A(C(u,v)) = (u,v) and C(A(x,y)) = (x,y)
A = ρ*C           # A(u, v) = ρ*C(u, v)
A = C*ρ           # A(u, v) = C(ρ*u, ρ*v)
A = C + (tu, tv)  # A(u, v) = C(u + tu, v + tv)
A = C + [tu, tv]  # A(u, v) = C(u + tu, v + tv)
A = (tx, ty) + C  # A(u, v) = C(u, v) +. (tx, ty)
A = [tx, ty] + C  # A(u, v) = C(u, v) +. [tx, ty]
```

the `-` operator can also be used to apply the opposite pre/post translation.

Two coordinate transforms `A` and `B` can be composed as follows:

```julia
C = A ∘ B
```

to form a new coordinate transform, `C` such that `C(u,v) = A(B(u,v))`.
Operators `\\` and  `/` can be used for left/right "division":

```julia
A\\B = inv(A) ∘ B
A/B = A ∘ inv(B)
```

Two coordinate transforms `A` and `B` can be approximately compared:

```julia
A ≈ B
```


See also: [`compose`](@ref), [`jacobian`](@ref), [`origin`](@ref),
          [`step`](@ref).

"""
struct CoordinateTransform
    stp::Float64  # step size
    x0::Float64   # abscissa at (u=0,v)
    y0::Float64   # ordinate at (u,v=0)
end

# Yields the identity without arguments.
CoordinateTransform() = CoordinateTransform(1, 0, 0)

"""
```julia
step(C)
```

yields the step size the coordinate transform `C` (see
[`CoordinateTransform`](@ref)).

"""
Base.step(C::CoordinateTransform) = C.stp

"""
```julia
jacobian(C) -> jac
```

yields the Jacobian of the coordinate transform `C`.

See also: [`CoordinateTransform`](@ref), [`jacobian`](@ref), [`origin`](@ref).

"""
jacobian(C::CoordinateTransform) = step(C)^2

"""
```julia
origin(C) -> (x0, y0)
```

yields the origin of coordinates for the coordinate transform `C`; that is
the coordinates given by `C(0,0)`.

See also: [`CoordinateTransform`](@ref), [`jacobian`](@ref), [`step`](@ref).

"""
origin(C::CoordinateTransform) = (C.x0, C.y0)

compose(A::CoordinateTransform, B::CoordinateTransform) =
    CoordinateTransform(A.stp*B.stp,
                        A.x0 + A.stp*B.x0,
                        A.y0 + A.stp*B.y0)

Base.:∘(A::CoordinateTransform, B::CoordinateTransform) = compose(A, B)

Base.isapprox(A::CoordinateTransform, B::CoordinateTransform; kwds...) =
    (isapprox(A.stp, B.stp; kwds...) &&
     isapprox(A.x0,  B.x0;  kwds...) &&
     isapprox(A.y0,  B.y0;  kwds...))

# Union of possible coordinate types.
const Coordinates = Union{NTuple{2,Real},AbstractVector{<:Real}}

@inline __checkcoordinates(c::AbstractVector{Real}) =
    (length(c) == 2 || throw(ArgumentError("expecting 2 coordinates"));
     nothing)

@inline __checkcoordinates(c::NTuple{2,Real}) = nothing

(C::CoordinateTransform)(u::Real, v::Real) = (u*C.stp + C.x0,
                                              v*C.stp + C.y0)

(C::CoordinateTransform)(uv::NTuple{2,Real}) = C(uv[1], uv[2])

(C::CoordinateTransform)(uv::AbstractVector{Real}) =
    (__checkcoordinates(uv); return Float64[uv[1]*C.stp + C.x0,
                                            uv[2]*C.stp + C.y0])

Base.inv(C::CoordinateTransform) =
    (stp = 1/C.stp; return CoordinateTransform(stp, -stp*C.x0, -stp*C.y0))

Base.:\(A::CoordinateTransform, B::CoordinateTransform) =
    compose(inv(A), B)

Base.:/(A::CoordinateTransform, B::CoordinateTransform) =
    compose(A, inv(B))

#for op in (:(*), :(⋅), :(∘))
#    @eval begin
#
#        Base.$op(C::CoordinateTransform, uv::Coordinates) = C(uv)
#
#        Base.$op(A::CoordinateTransform, B::CoordinateTransform) =
#            compose(A, B)
#
#    end
#end

Base.:*(C::CoordinateTransform, ρ::Real) =
    CoordinateTransform(C.stp*ρ, C.x0, C.y0)

Base.:*(ρ::Real, C::CoordinateTransform) =
    CoordinateTransform(ρ*C.stp, ρ*C.x0, ρ*C.y0)

Base.:+(C::CoordinateTransform, t::Coordinates) =
    (__checkcoordinates(t);
     return CoordinateTransform(C.stp,
                                C.x0 + C.stp*t[1],
                                C.y0 + C.stp*t[2]))

Base.:-(C::CoordinateTransform) =
    CoordinateTransform(-C.stp, -C.x0, -C.y0)

Base.:-(C::CoordinateTransform, t::Coordinates) =
    (__checkcoordinates(t);
     return CoordinateTransform(C.stp,
                                C.x0 - C.stp*t[1],
                                C.y0 - C.stp*t[2]))

Base.:+(t::Coordinates, C::CoordinateTransform) =
    (__checkcoordinates(t);
     return CoordinateTransform(C.stp,
                                C.x0 + t[1],
                                C.y0 + t[2]))

Base.:-(t::Coordinates, C::CoordinateTransform) =
    (__checkcoordinates(t);
     return CoordinateTransform(-C.stp,
                                t[1] - C.x0,
                                t[2] - C.y0))

end # module
