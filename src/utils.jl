"""
    FourierOptics.default_coordinates_center(inds)

yields the central index in the range of indices `inds` following the same
conventions as in `fftshift` and `ifftshift`.

"""
default_coordinates_center(inds::AbstractUnitRange{<:Integer}) =
    return first(inds) + (length(inds)÷2)

# Constructors for array coordinates and frequencies.
Coordinates(F::Field) = Coordinates(axes(F)[1])
Frequencies(F::Field) = Frequencies(axes(F)[1])
Coordinates(len::Integer) = Coordinates(Base.OneTo{Int}(len))
Coordinates(len::Integer, center::Integer) = Coordinates(Base.OneTo{Int}(len), center)
Frequencies(len::Integer) = Frequencies(Base.OneTo{Int}(len))

# Implement abstract array API for coordinates.
Base.length(A::AbstractCoordinates) = length(A.indices)
Base.axes(A::AbstractCoordinates) = (A.indices,)
Base.size(A::AbstractCoordinates) = (length(A),)
Base.IndexStyle(::Type{<:AbstractCoordinates}) = IndexLinear()
Base.firstindex(A::AbstractCoordinates) = first(A.indices)
Base.lastindex(A::AbstractCoordinates) = last(A.indices)
@inline function Base.getindex(A::Coordinates, i::Int)
    @boundscheck checkbounds(A, i)
    return i - A.center
end
@inline function Base.getindex(A::Frequencies, i::Int)
    @boundscheck checkbounds(A, i)
    return ifelse(i ≤ A.half, i - firstindex(A), i - lastindex(A) - 1)
end

"""
    FourierOptics.infinity(x)

yields the value corresponding to the positive infinity for the numerical
floating-point type or value `x`. Quantities are legitimate as far as their
bare numerical type is floating-point.

"""
infinity(x) = infinity(typeof(x))
infinity(::Type{T}) where {T<:AbstractFloat} = typemax(T)
infinity(::Type{Q}) where {Q<:Unitful.AbstractQuantity{T}} where {T<:AbstractFloat} = typemax(Q)
@noinline infinity(::Type{T}) where {T} = error("infinity is not defined for type `$T`")

"""
    FourierOptics.exp_i(ϕ) -> exp(i⋅ϕ)

yields `exp(i⋅ϕ)` computed efficiently. if Phase `ϕ` is not an angular
quantity, it must be a real or a dimensionless quantity assumed to be in
radians.

"""
function exp_i(ϕ::Union{Dimensionless{Real},Angle})
    sinϕ, cosϕ = sincos(ϕ)
    return complex(cosϕ, sinϕ)
end

function reset!(F::Field)
    ampl = get_amplitude(F)
    fill!(ampl, one(eltype(ampl)))
    return F
end

floating_point_type(A::Field) = floating_point_type(typeof(A))
floating_point_type(::Type{<:Field{T}}) where {T} = T

"""
    FourierOptics.standard_length([T = Float64,] x::Length)

yields a unitless floating-point value of type `T` corresponding to the length
`x` in meters.

"""
standard_length(x::Length) = standard_length(Float64, x)
standard_length(::Type{T}, x::Length) where {T<:AbstractFloat} = ustrip(T, m, x)

"""
    FourierOptics.standard_angle([T = Float64,] x::Angle)

yields a unitless floating-point value of type `T` corresponding to the angle
`x` in radians.

"""
standard_angle(x::Angle) = standard_angle(Float64, x)
standard_angle(::Type{T}, x::Angle) where {T<:AbstractFloat} = ustrip(T, rad, x)

"""
    FourierOptics.convert_multiplier(α, A)

converts multiplier `α` to the same floating-point type as the elements of the
array or field `A`.

"""
convert_multiplier(α::Number, A::AbstractArray) = convert_multiplier(α, typeof(A))
convert_multiplier(α::Real, ::Type{<:AbstractArray{T}}) where {T<:AbstractFloat} =
    as(T, α)
convert_multiplier(α::Real, ::Type{<:AbstractArray{Complex{T}}}) where {T<:AbstractFloat} =
    as(T, α)
convert_multiplier(α::Complex, ::Type{<:AbstractArray{Complex{T}}}) where {T<:AbstractFloat} =
    as(Complex{T}, α)
@noinline convert_multiplier(α::Number, A::Type{<:AbstractArray}) =
    throw(ArgumentError("invalid multiplier of type `$(typeof(α))` for array with elements of type `$(eltype(A))`"))

# Copy a field between two similar structures.
@inline copy_field!(dst::T, src, key::Symbol) where {T} =
    setfield!(dst, key, convert(fieldtype(T, key), getfield(src, key)))

"""
    FourierOptics.check_fftw_flags(flags)

checks whether `flags` is an allowed bitwise-or combination of FFTW planner
flags (see http://www.fftw.org/doc/Planner-Flags.html) and returns the filtered
flags.

"""
function check_fftw_flags(flags::Integer)
    iszero(flags & ~FFTW_PLANNING_FLAGS) || throw(ArgumentError(
        "only FFTW planning flags can be specified"))
    return as(UInt32, flags & FFTW_PLANNING_FLAGS)
end

get_input_size(P::FFTW.FFTWPlan) = P.sz
get_output_size(P::FFTW.FFTWPlan) = P.osz
get_flags(P::FFTW.FFTWPlan) = P.flags
