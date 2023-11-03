# All allowed planning flags.
const FFTW_PLANNING_FLAGS = (FFTW.ESTIMATE|FFTW.MEASURE|FFTW.PATIENT|FFTW.EXHAUSTIVE|FFTW.WISDOM_ONLY)

# Outer constructors.
Field(; kwds...) = Field{Float64}(; kwds...)
Field{T}(; grid_size::Integer, kwds...) where {T<:AbstractFloat} =
    Field(Array{Complex{T}}(undef, grid_size, grid_size); kwds...)
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
    T = real(eltype(ampl))
    A = typeof(ampl)
    F = typeof(forward)
    B = typeof(backward)
    return Field(ampl, forward, backward; kwds...)
end

# Copy constructor.
function Base.copy(Fin::Field)
    Fout = Field(copy(get_amplitude(Fin)), Fin.forward_plan, Fin.backward_plan;
                 wavelength = get_wavelength(Fin),
                 grid_step = get_grid_step(Fin),
                 verbosity = get_verbosity(Fin),
                 Rayleigh_factor = get_Rayleigh_factor(Fin),
                 phase_offset = get_phase_offset(Fin),
                 _preserve_amplitude = true)
    copy_field!(Fout, Fin, :w0)
    copy_field!(Fout, Fin, :z)
    copy_field!(Fout, Fin, :z_w0)
    copy_field!(Fout, Fin, :surface_type)
    return Fout
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
refected field `Fr` as if an infinitely thin beam-splitter with intensity
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

"""
    FourierOptics.get_amplitude(F::Field) -> amp

yields the complex amplitude of the field `F` at its reference surface.

"""
get_amplitude(F::Field) = F.ampl

# Implement the abstract array API for fields.
Base.length(F::Field) = length(get_amplitude(F))
Base.size(F::Field) = size(get_amplitude(F))
Base.axes(F::Field) = axes(get_amplitude(F))
Base.ndims(F::Field) = ndims(typeof(F))
Base.ndims(::Type{<:Field{T,A}}) where {T,A} = ndims(A)
Base.eltype(F::Field) = eltype(typeof(F))
Base.eltype(::Type{<:Field{T,A}}) where {T,A} = eltype(A)
Base.IndexStyle(::Type{<:Field{T,A,true}}) where {T,A} = IndexLinear()
@inline function Base.getindex(F::Field{T,A,true}, i::Int) where {T,A}
    @boundscheck checkbounds(F, i)
    return @inbounds getindex(get_amplitude(F), i)
end
@inline function Base.setindex!(F::Field{T,A,true}, x, i::Int) where {T,A}
    @boundscheck checkbounds(F, i)
    @inbounds setindex!(get_amplitude(F), x, i)
    return F
end
Base.IndexStyle(::Type{<:Field{T,A,false}}) where {T,A} = IndexCartesian()
@inline function Base.getindex(F::Field{T,A,false}, I::Vararg{Int,2}) where {T,A}
    @boundscheck checkbounds(F, I...)
    return @inbounds getindex(get_amplitude(F), I...)
end
@inline function Base.setindex!(F::Field{T,A,false}, x, I::Vararg{Int,2}) where {T,A}
    @boundscheck checkbounds(F, I...)
    @inbounds setindex!(get_amplitude(F), x, I...)
    return F
end

"""
    FourierOptics.get_intensity(F::Field)

yields the intensity of the field `F`.

"""
get_intensity(F::Field) = abs2.(get_amplitude(F))

"""
    FourierOptics.get_total_intensity(F::Field)

yields the total intensity of the field `F` at its reference surface.

"""
get_total_intensity(F::Field) = sum(get_intensity(F)) # FIXME optimize

"""
    FourierOptics.get_wavelength(F::Field) -> λ

yields the wavelength of field `F`.

"""
get_wavelength(F::Field) = F.lambda

function check_wavelength(λ::Length)
    (isfinite(λ) & (λ > zero(λ))) || error("wavelength must be positive and finite")
    return λ
end

"""
    FourierOptics.get_grid_size(F::Field) -> n

yields the number of samples along any dimension of the field `F`.

"""
function get_grid_size(F::Field)
    @noinline bad_grid_size(n1, n2) =
        AssertionError("invalid non-square grid size $n1 × $n2")
    n1, n2 = size(F)
    n1 == n2 || throw(bad_grid_size(n1, n2))
    return n1
end

"""
    FourierOptics.get_grid_range(F::Field) -> rng

yields the index range along any dimension of the field `F`.

"""
get_grid_range(F::Field) = Base.OneTo(get_grid_size(F))

"""
    FourierOptics.get_grid_step(F::Field) -> Δx

yields the current grid sampling step size in field `F`.

"""
get_grid_step(F::Field) = F.dx

function check_grid_step(dx::Length)
    (isfinite(dx) & (dx > zero(dx))) || error("grid step must be positive and finite")
    return dx
end

# Private function.
function _set_grid_step!(F::Field, dx::Length)
    F.dx = check_grid_step(dx)
    return nothing
end

"""
    FourierOptics.get_grid_step(F::Field, dz::Length) -> Δx′

yields the grid step after applying the *spherical-to-waist* or
*waist-to-spherical* propagators to the field `F` over a distance `dz`.
The result is given by:

    Δx′ = λ*abs(dz)/(n*Δx)

with `λ` the wavelength, `n` the number of grid nodes along a dimension, and
`Δx` the grid step **before** the propagation.

"""
get_grid_step(F::Field{T}, dz::Length) where {T} = get_grid_step(F, as(Meters{T}, dz))
function get_grid_step(F::Field{T}, dz::Meters{T}) where {T}
    n = get_grid_size(F)
    Δx = get_grid_step(F)
    λ = get_wavelength(F)
    return λ*abs(dz)/(n*Δx) :: Meters{T}
end

"""
    FourierOptics.get_position(F::Field) -> z

yields the position along the propagation axis of the field `F`.

"""
get_position(F::Field) = F.z

function _set_position!(F::Field, z::Length)
    F.z = z
    return nothing
end

_increment_position!(F::Field, dz::Length) = _set_position!(F, get_position(F) + dz)

"""
    FourierOptics.get_minimum_waist_radius(F::Field) -> w

yields the radius of the minimum beam waist for field `F`.

"""
get_minimum_waist_radius(F::Field) = F.w0

"""
    FourierOptics.get_minimum_waist_position(F::Field) -> z

yields the position along `z`-axis of the minimum beam waist for field `F`.

"""
get_minimum_waist_position(F::Field) = F.z_w0

"""
    FourierOptics.get_Rayleigh_distance(F::Field) -> d

yields the Rayleigh distance from minimum beam waist of field `F` and defined by:

    π⋅w₀^2/λ

with `w₀` the radius of the beam waist and `λ` the wavelength.

The Rayleigh distance from the beam waist (times some factor) specifies the
boundary of the near and far fields, and thus which reference surface type
(planar or spherical) is best used to minimize aliasing. Within the Rayleigh
distance of the waist, the wavefront can usually be fit best with a planar
reference surface, while outside it is best fit with a curved surface. The
radius of a reference sphere is equal to the distance from the current position
to the beam waist.

"""
function get_Rayleigh_distance(F::Field)
    λ  = get_wavelength(F)
    w0 = get_minimum_waist_radius(F)
    return π*w0^2/λ
end

"""
    FourierOptics.get_Rayleigh_factor(F::Field) -> f

yields the Rayleigh factor for the field `F`.

"""
get_Rayleigh_factor(F::Field) = F.Rayleigh_factor

function check_Rayleigh_factor(η::Real)
    (isfinite(η) & (η > zero(η))) || error("Rayleigh factor must be positive and finite")
    return η
end

function set_Rayleigh_factor!(F::Field, x::Real)
    F.Rayleigh_factor = check_Rayleigh_factor(x)
    return nothing
end

"""
    FourierOptics.get_waist_radius_at_surface(F::Field) -> w

yields the waist radius at the current surface of field `F`.

"""
function get_waist_radius_at_surface(F::Field{T}) where {T}
    z          = get_surface_position(F)
    z_w0       = get_minimum_waist_position(F)
    w0         = get_minimum_waist_radius(F)
    d_Rayleigh = get_Rayleigh_distance(F)
    return w0*sqrt(one(T) + ((z - z_w0)/d_Rayleigh)^2)
end

"""
    FourierOptics.get_focal_ratio(F::Field) -> fratio

yields the focal ratio of the current beam in field `F`:

    fratio = abs(z_w0 - z)/(2w)

with `z_w0` and `z` the positions along the propagation direction of the
minimum waist and of the current reference surface, and `w` the waist radius at
the reference surface.

"""
function get_focal_ratio(F::Field{T}) where {T}
    z    = get_surface_position(F)
    z_w0 = get_minimum_waist_position(F)
    w    = get_waist_radius_at_surface(F)
    return as(T, abs(z_w0 - z)/(2w))
end

"""
    FourierOptics.get_surface_type(F::Field)

yields the type (`:PLANAR` or `:SPHERICAL`) of the current reference surface of
the field `F`.

"""
get_surface_type(F::Field) = F.surface_type

function _set_surface_type!(F::Field, type::Symbol)
    ((type === :SPHERICAL) | (type === :PLANAR)) || throw(ArgumentError(
        "invalid reference surface type "))
    F.surface_type = type
    return nothing
end

"""
    FourierOptics.get_field_type(F::Field) -> rng

yields the last field range approximation used to propagate field `F` to its
current position. Returned value is `:NEAR` if the current reference surface is
spherical, or `:FAR` if the current surface is planar.

"""
get_field_type(F::Field) = (get_surface_type(F) === :PLANAR ? :NEAR : :FAR)

"""
    FourierOptics.get_field_type(F::Field, dz::Length) -> rng

yields the field range condition to propagate the field `F` by a distance `dz`.
Returned value is `:NEAR` for near field or `:FAR` for far field conditions.

"""
function get_field_type(F::Field, dz::Length)
    # Distance after propagation from the waist.
    distance = abs(get_surface_position(F) - get_minimum_waist_position(F) + dz)
    # Boundary between the near and far fields.
    boundary = get_Rayleigh_factor(F)*get_Rayleigh_distance(F)
    return distance < boundary ? :NEAR : :FAR
end

"""
    FourierOptics.get_verbosity(F::Field) -> verb

yields the verbosity lebel for computations with field `F`.

"""
get_verbosity(F::Field) = F.verbosity

function set_verbosity!(F::Field, verb::Integer)
    F.verbosity = verb
    return nothing
end

is_verbose(F::Field) = get_verbosity(F) > 0

propagate_phase_offset(F::Field) = get_phase_offset(F)

get_phase_offset(F::Field) = F.phase_offset

function set_phase_offset!(F::Field, val::Bool)
    F.phase_offset = val
    return nothing
end

"""
    FourierOptics.multiply(Fin, α) -> Fout

multiplies the complex amplitude of the field `Fin` by `α` and returns the
resulting field `Fout`. Argument may also be any type of abstract array.

"""
multiply(A::AbstractArray, α::Union{Real,Complex}) = multiply!(copy(F), α)

"""
    FourierOptics.multiply!(F, α) -> F

multiplies in-place the complex amplitude of the field `F` by `α`.

"""
function multiply!(F::Field, α::Union{Real,Complex})
    multiply!(get_amplitude(F), α)
    return F
end

function multiply!(A::AbstractArray, α::Union{Real,Complex})
    if iszero(α)
        fill!(A, zero(eltype(A)))
    elseif !isone(α)
        α = convert_multiplier(α, A)
        @inbounds @simd for i in eachindex(A)
            A[i] *= α
        end
    end
    return A
end
