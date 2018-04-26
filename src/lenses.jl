#
# lenses.jl -
#
# Implement lenses model.
#

Base.show(io::IO, ::MIME"text/plain", L::Lens) = show(io, L)
Base.show(io::IO, ::MIME"text/plain", A::LensOperator) = show(io, A)

function Base.show(io::IO, L::Lens)
    print(io, "Lens(f = ")
    showlength(io, focal_length(L))
    print(io, ", D = ")
    showlength(io, diameter(L))
    print(io, ", center = ")
    showposition(io, center(L))
    print(io, ")")
end

function Base.show(io::IO, A::LensOperator)
    print(io, "LensOperator(")
    show(io, A.lens)
    print(io, ", λ = ")
    showlength(io, A.lambda)
    print(io, ", inp = ")
    show(io, A.Rinp)
    print(io, ", out = ")
    show(io, A.Rout)
    print(io, ")")

end

"""

An instance of `Lens` stores the parameters of a circular lens defined by its
focal length `f`, its diameter `D` and the coordinates `(x0,y0)` of the center
of the lens in a plane perpendicular to the direction of propagation.

A lens can be created by a keyword-only constructor:

```julia
Lens(; f=focal_length; D=diameter; center=(x0,y0))
```

the default center is `(0mm,0mm)`, the default focal length is `30mm` and the
default diameter is `20mm`.

"""
Lens(; f::Real=30mm, D::Real=20mm, center::NTuple{2,Real} = (0mm, 0mm)) =
    Lens(f, D, center[1], center[2])

"""
```julia
focal_length(L) -> f
```

yields the focal length `f` of the lens `L` in SI units (meters).

"""
focal_length(L::Lens) = L.f

"""
```julia
diameter(L) -> D
```

yields the diameter `D` of the circular lens `L` in SI units (meters).

"""
diameter(L::Lens) = L.D

"""
```julia
center(L) -> (x0, y0)
```

yields the coordinates `(x0, y0)` of the center of the lens `L` in its pupil
plane (in SI units, that is meters).

"""
center(L::Lens) = (L.x0, L.y0)

"""
```julia
LensOperator(L, P, λ, n)
```

yields a linear operator to propagate the complex amplitude from a pupil plane
`P` through a lens `L` to its focal plane.  `λ` is the wavelength (in SI units)
and `n` is the size of the output region in the focal plane and centered at the
lens center.

"""
function LensOperator(L::Lens, Rinp::Region, lambda::Real, n::Integer)
    # Get the bounding box of the part of the input which encompass the pupil.
    # FIXME: There seem to be 1-2 extra pixels.
    dr = step(Rinp)
    if diameter(L) > n*dr/2
        error("there will be aliasing")
    end
    n1, n2 = size(Rinp)
    x0, y0 = center(L)
    r1 = (diameter(L) - dr)/2
    r2 = (diameter(L) + dr)/2
    i0, j0 = floor.(Int, world2grid(Rinp, x0 - r2, y0 - r2))
    i1, j1 =  ceil.(Int, world2grid(Rinp, x0 + r2, y0 + r2))
    if ! (1 ≤ i0 ≤ i1 ≤ n1 && 1 ≤ j0 ≤ j1 ≤ n2)
        error("there is some vignetting!")
    end

    # Compute the offset (in fractional samples) of the lens center w.r.t. the
    # origin of the extracted part (which is also the origin for the FFT).
    inpoff1, inpoff2 = world2grid(Rinp, x0, y0) .- (i0, j0)

    # Define the output region so that it is centered on the window computed by
    # FFT and has the same size.
    ds = convert(Float64, focal_length(L)*lambda/(n*dr)) # ouput step
    outoff1 = outoff2 = (n - 1)/2
    s = ds*(n - 1)/2 # half-width of output region between first and last nodes
    Xout = linspace(-s, +s, n)
    Yout = linspace(-s, +s, n)
    Rout = Region(Xout, Yout)

    # Pre-phase modulation to translate the output.
    preshift1 = fftshiftphasor(outoff1, n)
    preshift2 = fftshiftphasor(outoff2, n)

    # Post-phase modulation to compensate the off-centering of the pupil.
    η = 2π*dr/(focal_length(L)*lambda)
    postshift1 = phasor(Complex{Float64}, η*inpoff1, Xout)
    postshift2 = phasor(Complex{Float64}, η*inpoff2, Yout)

    alpha = Float64(dr^2)/Float64(focal_length(L)*lambda)
    Qinp = Array{Complex{Float64}}(i1 - i0 + 1, j1 - j0 + 1)
    k = 0
    for j in j0:j1
        for i in i0:i1
            x, y = grid2world(Rinp, i, j)
            r = hypot(x - x0, y - y0)
            k += 1
            if r ≥ r2
                Qinp[k] = 0
            else
                rho = (r ≤ r1 ? alpha : alpha*(r2 - r)/dr)
                Qinp[k] = rho*preshift1[i - i0 + 1]*preshift2[j - j0 + 1]
            end
        end
    end

    Qout = Array{Complex{Float64}}(n, n)
    @inbounds for j in 1:n
        @simd for i in 1:n
            Qout[i,j] = postshift1[i]*postshift2[j]
        end
    end

    # Allocate a workspace and a plan for in-place FFT.
    ws = Array{Complex{Float64}}(n, n)
    F = plan_fft!(ws; flags=FFTW.MEASURE)

    return LensOperator(L, lambda,
                        Rinp, Qinp, IndexBox(i0:i1,j0:j1),
                        Rout, Qout,
                        F, ws)

end

function LazyAlgebra.vcreate(::Type{Direct}, A::LensOperator,
                             x::AbstractMatrix{<:Complex{<:AbstractFloat}})
    return similar(A.ws)
end

function LazyAlgebra.apply!(alpha::Real, ::Type{Direct}, A::LensOperator,
                            x::AbstractMatrix{<:Complex{<:AbstractFloat}},
                            beta::Real,
                            y::AbstractMatrix{<:Complex{<:AbstractFloat}})
    if size(x) != size(A.Rinp)
        throw(DimensionMismatch("source size $(size(x)) must be $(size(A.Rinp))"))
    end
    if size(y) != size(A.Rout)
        throw(DimensionMismatch("destination size $(size(y)) must be $(size(A.Rout))"))
    end
    if alpha == 0
        vscale!(y, beta)
    else
        # Get indices ranges of the relevant input part and offset them
        # w.r.t. workspace arrays ws and Qinp.
        I, J = A.Binp.I, A.Binp.J
        if ! (1 ≤ first(I) ≤ last(I) ≤ size(A.Rinp,1) &&
              1 ≤ first(J) ≤ last(J) ≤ size(A.Rinp,2) &&
              size(A.Qinp) == (length(I), length(J)) &&
              size(A.ws) == size(A.Qout) == size(A.Rout))
            throw(ArgumentError("corrupted LensOperator structure!"))
        end
        i0 = first(I) - 1
        j0 = first(J) - 1
        I -= i0
        J -= j0

        # Compute transmitted complex amplitude in the pupil plane.
        ws = vzero!(A.ws)
        Qinp = A.Qinp
        @inbounds for j in J
            @simd for i in I
                ws[i,j] = Qinp[i,j]*x[i+i0,j+j0]
            end
        end

        # Compute complex amplitude in the focal plane.
        println(size(ws))
        A_mul_B!(ws, A.FFT, ws)
        Qout = A.Qout
        if beta == 0
            if alpha == 1
                @inbounds @simd for k in eachindex(Qout, y, ws)
                    y[k] = Qout[k]*ws[k]
                end
            else
                α = convert(T, alpha)
                @inbounds @simd for k in eachindex(Qout, y, ws)
                    y[k] = α*Qout[k]*ws[k]
                end
            end
        elseif beta == 1
            if alpha == 1
                @inbounds @simd for k in eachindex(Qout, y, ws)
                    y[k] += Qout[k]*ws[k]
                end
            else
                α = convert(T, alpha)
                @inbounds @simd for k in eachindex(Qout, y, ws)
                    y[k] += α*Qout[k]*ws[k]
                end
            end
        else
            β = convert(T, beta)
            if alpha == 1
                @inbounds @simd for k in eachindex(Qout, y, ws)
                    y[k] = Qout[k]*ws[k] + β*y[k]
                end
            else
                α = convert(T, alpha)
                @inbounds @simd for k in eachindex(Qout, y, ws)
                    y[k] = α*Qout[k]*ws[k] + β*y[k]
                end
            end
        end
    end
    return y
end

LazyAlgebra.input_type(A::LensOperator) = AbstractMatrix{Complex{Float64}}
LazyAlgebra.input_ndims(A::LensOperator) = 2
LazyAlgebra.input_size(A::LensOperator) = size(A.Rinp)
LazyAlgebra.input_size(A::LensOperator, d...) = size(A.Rinp, d...)
LazyAlgebra.input_eltype(A::LensOperator) = Complex{Float64}

LazyAlgebra.output_type(A::LensOperator) =  AbstractMatrix{Complex{Float64}}
LazyAlgebra.output_ndims(A::LensOperator) = 2
LazyAlgebra.output_size(A::LensOperator) = size(A.Rout)
LazyAlgebra.output_size(A::LensOperator, d...) = size(A.Rout, d...)
LazyAlgebra.output_eltype(A::LensOperator) = Complex{Float64}
