#
# utils.jl -
#
# General purpose functions.
#

"""

```julia
circularmask(R, [rmin = 0,] rmax) -> M
```

yields an array `M` of same dimensions (and element type) as `R` which is
set with a circular mask at radii `R`.  Argument `rmax` is the outer radius
of the mask while optional argument `rmin` is the inner radius of the mask
(0 by default).

The value of the result is 1 inside the mask, 0 ouside the mask and in the
range `[0,1]` at the mask borders.  Keyword `rev` can be set true to
reverse the mask, that is to set the result to 0 inside the mask and 1
outside the mask.

Keyword `shape` can be set with a non-zero value to use simple
anti-aliasing rules for the edges.  The thickness (in same units as `R`) of
the mask border is given by the absolute value of `shape`.  Thus, the
support of the mask has a radius: `rmax + abs(shape)/2`.  With `shape > 0`,
a linear ramp is used; with `shape = 0` (the default), a rectangular
profile is used; with `shape < 0`, a sine profile is used.

The destination mask may be specified:

```julia
circularmask!(M, R, [rmin = 0,] rmax) -> M
```

"""
circularmask(R::AbstractArray{<:AbstractFloat}, args...; kwds...) =
     circularmask!(similar(R), R, args...; kwds...)

function circularmask!(M::AbstractArray{Tm,N},
                       R::AbstractArray{Tr,N},
                       r1::Real, r2::Real=0;
                       shape::Real = 0,
                       rev::Bool = false) where {Tm<:AbstractFloat,
                                                 Tr<:AbstractFloat, N}
    @assert size(M) == size(R)
    const rmin = convert(Tr, min(r1, r2))
    const rmax = convert(Tr, max(r1, r2))
    const inside = (rev ? zero(Tm) : one(Tm))  # value inside mask
    const outside = (rev ? one(Tm) : zero(Tm)) # value ouside mask
    const edge = (inside + outside)/2          # value at middle of edge

    if shape == 0
        # Points exactly at the edges are set to 1/2.
        if rmin > 0
            # Account for a central obscuration.
            @inbounds for i in eachindex(M, R)
                r = R[i]
                M[i] = (r < rmin || r > rmax ? outside :
                        rmin < r < rmax ? inside : edge)
            end
        else
            # No central obscuration.
            @inbounds for i in eachindex(M, R)
                r = R[i]
                M[i] = (r > rmax ? outside :
                        r < rmax ? inside  : edge)
            end
        end
    else
        # Points near the edges are set to a value between 0 and 1
        # following a given profile.  Note that the expressions at the
        # edges are centered to minimize rounding errors.
        const h = abs(convert(Tr, shape))/2 # half border thickness
        const lim1 = rmin - h
        const lim2 = (rmin + h ≤ rmax - h ? rmin + h : (rmin + rmax)/2)
        const lim3 = (rmin + h ≤ rmax - h ? rmax - h : (rmin + rmax)/2)
        const lim4 = rmax + h
        if shape > 0
            # Use a linear B-spline for borders.
            const a = (inside - outside)/(h + h)
            if lim2 > 0
                # Account for a central obscuration.
                @inbounds for i in eachindex(M, R)
                    r = R[i]
                    M[i] = (r ≥ lim4 ? outside             :
                            r > lim3 ? (rmax - r)*a + edge :
                            r ≥ lim2 ? inside              :
                            r > lim1 ? (r - rmin)*a + edge : outside)
                end
            else
                # No central obscuration.
                @inbounds for i in eachindex(M, R)
                    r = R[i]
                    M[i] = (r ≥ lim4 ? outside             :
                            r > lim3 ? (rmax - r)*a + edge : inside)
                end
            end
        else # shape < 0
            # Use a sine window for borders.
            const a = convert(Tr, π)*(inside - outside)/(h + h)
            if lim2 > 0
                # Account for a central obscuration.
                @inbounds for i in eachindex(M, R)
                    r = R[i]
                    M[i] = (r ≥ lim4 ? outside                       :
                            r > lim3 ? sin((rmax - r)*a)*edge + edge :
                            r ≥ lim2 ? inside                        :
                            r > lim1 ? sin((r - rmin)*a)*edge + edge : outside)
                end
            else
                # No central obscuration.
                @inbounds for i in eachindex(M, R)
                    r = R[i]
                    M[i] = (r ≥ lim4 ? outside                       :
                            r > lim3 ? sin((rmax - r)*a)*edge + edge : inside)
                end
            end
        end
    end
    return M
end

@doc @doc(circularmask) circularmask!

"""
```julia
fftshiftphasor([T,] s, n)
```

yields a vector of `n` complexes whose elements have values `exp(2*i*π*k*s/n)`
for `k ∈ fftfreq(n)`.  Optional argument `T <: Complex{<:AbstractFloat}` is the
type of the elements of the result.  By default, `T =
Complex{float(typeof(s))}`.

"""
fftshiftphasor(s::T, n::Integer) where {T<:Real} =
    fftshiftphasor(Complex{float(T)}, s, n)

function fftshiftphasor(::Type{Complex{T}}, s::Real,
                        n::Integer) where {T<:AbstractFloat}
    fftshiftphasor!(Array{Complex{T}}(n), s)
end

function fftshiftphasor!(dst::AbstractVector{Complex{T}},
                         s::Real) where {T<:AbstractFloat}
    # For fractional shifts, it is important to use the signed discrete
    # frequency index.
    len = length(dst)
    n = div(len, 2)
    α = convert(T, 2π*s/len)
    @inbounds @simd for k in 1:len-n
        ϕ = α*(k - 1)
        dst[k] = complex(cos(ϕ), sin(ϕ))
    end
    l = len + 1
    @inbounds @simd for k in len-n+1:len
        ϕ = α*(k - l)
        dst[k] = complex(cos(ϕ), sin(ϕ))
    end
    if iseven(len)
        # Treat the Nyquist frequency specifically for a slight gain in
        # precision.
        k = n + 1
        dst[k] = real(dst[k])
    end
    return dst
end

phasor(phi::AbstractArray{T,N}) where {T<:Real,N} =
    phasor(Complex{float(T)}, phi)

function phasor(::Type{Complex{T}},
                phi::AbstractArray{<:Real,N}) where {T<:AbstractFloat,N}
    phasor!(Array{Complex{T}}(size(phi)), phi)
end

function phasor!(dst::AbstractArray{Complex{T},N},
                 phi::AbstractArray{<:Real,N}) where {T<:AbstractFloat,N}
    @assert size(dst) == size(phi)
    @inbounds @simd for k in eachindex(dst, phi)
        ϕ = convert(T, phi[k])
        dst[k] = complex(cos(ϕ), sin(ϕ))
    end
    return dst
end

phasor(alpha::Ta, X::AbstractArray{Tx,N}) where {Ta<:Real,Tx<:Real,N} =
    phasor(Complex{float(promote_type(Ta, Tx))}, alpha, X)

function phasor(::Type{Complex{T}},
                alpha::Real,
                X::AbstractArray{<:Real,N}) where {T<:AbstractFloat,N}
    phasor!(Array{Complex{T}}(size(X)), alpha, X)
end

function phasor!(dst::AbstractArray{Complex{T},N},
                 alpha::Real,
                 X::AbstractArray{<:Real,N}) where {T<:AbstractFloat,N}
    @assert size(dst) == size(X)
    α = convert(T, alpha)
    @inbounds @simd for k in eachindex(dst, X)
        ϕ = α*convert(T, X[k])
        dst[k] = complex(cos(ϕ), sin(ϕ))
    end
    return dst
end
