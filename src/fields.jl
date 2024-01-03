# All allowed planning flags.
const FFTW_PLANNING_FLAGS = (FFTW.ESTIMATE|FFTW.MEASURE|FFTW.PATIENT|FFTW.EXHAUSTIVE|FFTW.WISDOM_ONLY)

"""
    FourierOptics.Field(ampl; refractive_index=1, wavelength, sampling, kwds...)
    FourierOptics.Field(copy(ampl); refractive_index=1, wavelength, sampling, kwds...)

constructs a new field with complex amplitude stored in array `ampl`. This
array is not modified by the constructor and is shared by the object (the
caller can make a copy first to avoid this, as shown in the second example).
Other parameters are specified by keywords (some are mandatory):
`refractive_index` is the refractive index of the propagation medium,
`wavelength` is the wavelength in the vacuum, and `sampling` is the lateral grid
step.

Other possible keywords in `kwds...` are `fftw_flags`, and `fftw_timelimit` the
flags and time-limit for creating the plans for FFTW. If the caller has already
created such plans, the constructor may be called as:

    FourierOptics.Field(ampl, forward, backward; refractive_index=1, wavelength, sampling)

with `forward` and `backward` the FTTW plans for performing the forward and
backward FFT.

---

   FourierOptics.Field{T=Float64}(; refractive_index=1, wavelength, dim, sampling, fill=0, kwds...)

constructs a new field with complex amplitude stored in an array whose elements
have type `Complex{T}`. All parameters are specified by keywords (see above).
Additional keywords are `dim` the number of samples along each dimension of a
tranversal plane and `fill` the initial complex amplitude of the field.

"""
Field(; kwds...) = Field{Float64}(; kwds...)
Field{T}(; fill=zero(Complex{T}), dim::Integer, kwds...) where {T<:AbstractFloat} =
    Field(fill!(Array{Complex{T}}(undef, dim, dim), fill); kwds...)
function Field(ampl::StridedMatrix{<:FFTW.fftwComplex};
               fftw_timelimit::Real = default_fftw_timelimit,
               fftw_flags::Integer = default_fftw_flags,
               kwds...)
    # Create the FFTW plans with suitable flags. For maximum efficiency, the
    # transforms are always applied in-place and thus cannot preserve their
    # inputs.
    flags = check_fftw_flags(fftw_flags) | FFTW.DESTROY_INPUT
    forward = plan_fft!(ampl; flags, timelimit = fftw_timelimit)
    backward = plan_bfft!(ampl; flags, timelimit = fftw_timelimit)

    # Call inner constructor.
    return Field(ampl, forward, backward; kwds...)
end

# Copy constructor.
Base.copy(F::Field) = copyto!(similar(F), F)

function Base.copyto!(dst::Field, src::Field)
    @assert_same_axes dst src
    copyto!(amplitude(dst), amplitude(src))
    copy_struct_field!(dst, src, :n)
    copy_struct_field!(dst, src, :λ₀)
    copy_struct_field!(dst, src, :z)
    copy_struct_field!(dst, src, :fact)
    copy_struct_field!(dst, src, :curv)
    return dst
end

Base.similar(F::Field, ::Type{T}, dims::Dims{2}) where {T} =
    Field(similar(amplitude(F), T, dims),
          F.forward_plan, F.backward_plan;
          refractive_index = refractive_index(F),
          wavelength = wavelength_in_vacuum(F),
          sampling = grid_step(F))

"""
    FourierOptics.reset!(F::Field) -> F

resets the complex amplitude of the field `F` and returns `F`.

"""
function reset!(F::Field)
    ampl = amplitude(F)
    fill!(ampl, one(eltype(ampl)))
    F.curv = zero(F.curv)
    F.fact = one(F.fact)
    return F
end

# Implement the abstract array API for fields.
Base.length(F::Field) = length(amplitude(F))
Base.size(F::Field) = size(amplitude(F))
Base.axes(F::Field) = axes(amplitude(F))
Base.ndims(F::Field) = ndims(typeof(F))
Base.ndims(::Type{<:Field{T,A}}) where {T,A} = ndims(A)
Base.eltype(F::Field) = eltype(typeof(F))
Base.eltype(::Type{<:Field{T,A}}) where {T,A} = eltype(A)
Base.IndexStyle(::Type{<:Field{T,A,true}}) where {T,A} = IndexLinear()
@inline function Base.getindex(F::Field{T,A,true}, i::Int) where {T,A}
    @boundscheck checkbounds(F, i)
    return @inbounds getindex(amplitude(F), i)
end
@inline function Base.setindex!(F::Field{T,A,true}, x, i::Int) where {T,A}
    @boundscheck checkbounds(F, i)
    @inbounds setindex!(amplitude(F), x, i)
    return F
end
Base.IndexStyle(::Type{<:Field{T,A,false}}) where {T,A} = IndexCartesian()
@inline function Base.getindex(F::Field{T,A,false}, I::Vararg{Int,2}) where {T,A}
    @boundscheck checkbounds(F, I...)
    return @inbounds getindex(amplitude(F), I...)
end
@inline function Base.setindex!(F::Field{T,A,false}, x, I::Vararg{Int,2}) where {T,A}
    @boundscheck checkbounds(F, I...)
    @inbounds setindex!(amplitude(F), x, I...)
    return F
end

Base.show(io::IO, ::MIME"text/plain", F::Field) = show(io, F)
function Base.show(io::IO, F::Field{T}) where {T}
    n1, n2 = size(F)
    print(io, n1, "×", n2, " Field{", T,
          "}: n = ", round(refractive_index(F), sigdigits=4),
          ", λ₀ = ", round(ustrip(nm, wavelength_in_vacuum(F)), sigdigits=6),
          "nm,\n    δx = ", round(ustrip(mm, grid_step(F)), sigdigits=4),
          "mm, z = ")
    z = F.z
    if abs(z) < 10m
        print(io, round(ustrip(mm, z), sigdigits=8), "mm, R = ")
    else
        print(io, round(ustrip(m, z), sigdigits=8), "m, R = ")
    end
    R = inv(curvature(F))
   if isinf(R)
        print(io, "+∞")
    else
        print(io, round(ustrip(m, R), sigdigits=4))
    end
    print(io, " m, Iₜₒₜ = ", total_intensity(F))
end

"""
    FourierOptics.curvature(F::Field) -> 1/R

yields the wavefront curvature of the field `F`.

"""
curvature(F::Field) = F.curv

"""
    FourierOptics.multiplier(F::Field) -> 1/R

yields the uniform factor for the field `F`.

"""
multiplier(F::Field) = F.fact

"""
    FourierOptics.amplitude(F::Field) -> amp

yields the complex amplitude of the field `F` sampled in the lateral grid. The
returned array is the complex amplitude as stored in `F`, not accounting for a
possible wavefront curvature and uniform factor. Call
[`FourierOptics.true_amplitude`](@ref) to retrieve the complex amplitude with
all these terms.

"""
amplitude(F::Field) = F.ampl

"""
    FourierOptics.true_amplitude(F::Field) -> amp

yields the complex amplitude of the field `F` sampled in the lateral grid.
Compared to [`FourierOptics.amplitude`](@ref), the returned array accounts for
a possible wavefront curvature and uniform factor.

The true complex amplitude of the field `F` is given by:

    η*amp[j1,j2]*exp_i((π/(λ*R))⋅x[j1]^2)*exp_i((π/(λ*R))⋅y[j2]^2)

with `amp = amplitude(F)` the complex amplitude as stored in `F`, `η =
multipler(F)` a uniform factor, `λ = wavelength_in_medium(F)` the wavelength in
the medium, and `1/R = curvature(F)` the wavefront curvature.

"""
true_amplitude(F::Field) = true_amplitude!(similar(amplitude(F)), F)

function true_amplitude!(dst::AbstractMatrix, F::Field{T}) where {T}
    amp = amplitude(F)
    J₁, J₂ = axes(amp)
    axes(dst) == axes(J₁, J₂) || throw(DimensionMismatch(
        "destination array has incompatible axes"))
    curv = curvature(F)
    η = multiplier(F)
    if iszero(curv)
        if η == oneunit(η)
            copyto!(dst, amp)
        else
            mul!(dst, η, amp)
        end
    else
        δx = grid_step(F)
        λ = wavelength_in_medium(F)
        x = RolledCoordinates(F)
        πρ = π*as(T, δx^2/λ*curv)
        exp_iπρx² = Vector{Complex{T}}(undef, size(x))
        axes1(exp_iπρx²) == J₁ == J₂ || throw(DimensionMismatch(
            "complex amplitude has unexpected axes"))
        @inbounds for j in eachindex(exp_iπρx², x)
            exp_iπρx²[j] = exp_i(πρ*x[j]^2)
        end
        @inbounds for j₂ in J₂
            ξ = η*exp_iπρx²[j₂]
            for j₁ in J₁
                dst[j₁,j₂] = ξ*exp_iπρx²[j₁]
            end
        end
    end
    return dst
end

"""
    FourierOptics.intensity(F::Field) -> I

yields the intensity of the field `F` sampled in the lateral grid.

"""
intensity(F::Field) =
    intensity!(similar(amplitude(F), real(eltype(F))), F)

"""
    FourierOptics.intensity!(dst, F::Field) -> dst

overwrites `dst` with the intensity of the field `F` sampled in the lateral
grid.

"""
function intensity!(dst::AbstractMatrix, F::Field)
    A = amplitude(F)
    @assert_same_axes dst A
    η = abs2(F.fact)
    if iszero(η)
        fill!(dst, zero(eltype(dst)))
    elseif isone(η)
        @inbounds @simd for i in eachindex(dst, A)
            dst[i] = abs2(A[i])
        end
    else
        @inbounds @simd for i in eachindex(dst, A)
            dst[i] = η*abs2(A[i])
        end
    end
    return dst
end

"""
    FourierOptics.total_intensity(F::Field)

yields the total intensity of the field `F` at its current position.

"""
function total_intensity(F::Field{T}) where {T}
    # FIXME: Assume area in SI units.
    δx = as(T, grid_step(F)/oneunit(StdLength{T}))
    η = abs2(F.fact*δx)
    s = zero(η)
    if !iszero(η)
        A = amplitude(F)
        @inbounds @simd for i in eachindex(A)
            s += abs2(A[i])
        end
    end
    return η*s
end

"""
    FourierOptics.refractive_index(F::Field) -> n

yields the refractive index of the medium for the field `F`.

"""
refractive_index(F::Field) = F.n

function check_refractive_index(n::Real)
    (isfinite(n) & (n > zero(n))) || error("refractive index must be positive and finite")
    return n
end

"""
    FourierOptics.wavelength_in_vacuum(F::Field) -> λ₀

yields the wavelength of field `F` in the vacuum.

"""
wavelength_in_vacuum(F::Field) = F.λ₀

"""
    FourierOptics.wavelength_in_medium(F::Field) -> λ

yields the wavelength of field `F` in the medium.

"""
wavelength_in_medium(F::Field) =
    wavelength_in_vacuum(F)/refractive_index(F)

"""
    FourierOptics.wavenumber(F::Field) -> k

yields the wave-number `k = 2⋅π⋅n/λ₀ = 2⋅π/λ` for the field `F`.

"""
wavenumber(F::Field) = 2*refractive_index(F)*π/wavelength_in_vacuum(F)

function check_wavelength(λ::Length)
    (isfinite(λ) & (λ > zero(λ))) || error("wavelength must be positive and finite")
    return λ
end

"""
    FourierOptics.grid_size(F::Field) -> N

yields the number of samples along any lateral dimension of the field `F`.

"""
function grid_size(F::Field)
    @noinline bad_grid_size(n1, n2) =
        AssertionError("invalid non-square grid size $n1 × $n2")
    n1, n2 = size(F)
    n1 == n2 || throw(bad_grid_size(n1, n2))
    return n1
end

"""
    FourierOptics.grid_range(F::Field) -> rng

yields the index range along any lateral dimension of the field `F`.

"""
grid_range(F::Field) = Base.OneTo(grid_size(F))

"""
    FourierOptics.grid_step(F::Field) -> δx

yields the current grid sampling step size in field `F`.

"""
grid_step(F::Field) = F.δx

function check_grid_step(δx::Length)
    (isfinite(δx) & (δx > zero(δx))) || error("grid step must be positive and finite")
    return δx
end

# Private function.
function _set_grid_step!(F::Field, δx::Length)
    F.δx = check_grid_step(δx)
    return nothing
end

"""
    Fout = Fin*η
    Fout = η*Fin
    Fout = FourierOptics.multiply(Fin, η)

multiplies the complex amplitude of the field `Fin` by `η` and returns the
resulting field `Fout` leaving `Fin` unchanged. Argument `η` may be a scalar or
any type of abstract array of the same size as the field. In this latter case,
the multiplication is performed elementwise. This may be used to apply a
transmission mask (in amplitude) to the field.

See [`FourierOptics.multiply!`](@ref) for an in-place version.

"""
multiply(F::Field, η::Union{T,AbstractMatrix{T}}) where {T<:Union{Real,Complex}} =
    multiply!(copy(F), η)
Base.:(*)(F::Field, η::Union{T,AbstractMatrix{T}}) where {T<:Union{Real,Complex}} =
    multiply(F, η)
Base.:(*)(η::Union{T,AbstractMatrix{T}}, F::Field) where {T<:Union{Real,Complex}} =
    multiply(F, η)

"""
    FourierOptics.multiply!(F, η) -> F

multiplies in-place the complex amplitude of the field `F` by `η`. Argument `η`
may be a scalar or any type of abstract array of the same size as the field

See [`FourierOptics.multiply`](@ref) for an out-place version and for more
details.

"""
function multiply!(F::Field, η::Union{Real,Complex})
    F.fact *= η # in general this is sufficient
    if iszero(F.fact)
        # Multiplying by 0 is a shortcut to set the complex amplitude to zero
        # everywhere.
        fill!(F.ampl, zero(eltype(F.ampl)))
    end
    return F
end

function multiply!(F::Field, η::AbstractArray{<:Union{Real,Complex}})
    @assert_same_axes F η
    @inbounds @simd for i in eachindex(F, η)
            F[i] *= η[i]
    end
    return F
end

"""
    FourierOptics.split_beam(F, τ) -> Ft, Fr

splits the input field `F` into a transmitted field `Ft` and a reflected field
`Fr` as if an infinitely thin beam-splitter with intensity transmission
coefficient `τ` is inserted at the curent propagation position of the field
`F`.

"""
split_beam(F::Field, τ::Real) = split_beam!(copy(F), τ)

"""
    FourierOptics.split_beam!(F, τ) -> Ft, Fr

splits in-place the input field `F` into a transmitted field `Ft` and a
reflected field `Fr` as if an infinitely thin beam-splitter with intensity
transmission coefficient transmission `τ` is inserted at the curent propagation
position of the field `F`. The operation is done in-place: the input field `F`
is modified and returned as one the resulting fields.

"""
function split_beam!(F::Field{T}, τ::Real) where {T}
    # Compute the complex amplitude transmission and reflection coefficients so
    # as to preserve energy.
    zero(τ) ≤ τ ≤ one(τ) || throw(ArgumentError("transmission coefficient must be in [0,1]"))
    t = sqrt(as(T, τ))
    r = sqrt(one(T) - as(T, τ))
    # Compute the output fileds by simple multiplication. NOTE in-place
    # multiplication must be the last.
    Ft = multiply(F, t) # transmitted field (out-of-place multiplication)
    Fr = multiply!(F, r) # reflected field (in-place multiplication)
    return Ft, Fr
end
