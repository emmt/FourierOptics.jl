"""
    FourierOptics.default_coordinates_center(inds)

yields the central index in the range of indices `inds` following the same
conventions as in `fftshift` and `ifftshift`.

"""
default_coordinates_center(inds::AbstractUnitRange{<:Integer}) =
    return first(inds) + (length(inds)÷2)

# Constructors for array coordinates and frequencies.
Coordinates(F::Field) = Coordinates(axes(F)[1])
RolledCoordinates(F::Field) = RolledCoordinates(axes(F)[1])
Coordinates(len::Integer) = Coordinates(Base.OneTo{Int}(len))
Coordinates(len::Integer, center::Integer) = Coordinates(Base.OneTo{Int}(len), center)
RolledCoordinates(len::Integer) = RolledCoordinates(Base.OneTo{Int}(len))

Base.show(io::IO, ::MIME"text/plain", x::AbstractCoordinates) = show(io, x)
Base.show(io::IO, x::Coordinates) = print(io, "Coordinates(", first(x), ":", last(x), ")")
Base.show(io::IO, x::RolledCoordinates) = print(io, "RolledCoordinates(", length(x), ")")
Base.show(io::IO, x::ScaledCoordinates) = begin
    #=
    print(io, "ScaledCoordinates(", x.step, "*")
    show(io, x.coords)
    print(io, ")")
    =#
    print(io, x.step, "*")
    show(io, x.coords)
end

Base.:(*)(α::Number, x::ScaledCoordinates) = (α*x.step)*x.coords
Base.:(*)(α::Number, x::AbstractCoordinates) = ScaledCoordinates(α, x)
Base.:(*)(x::AbstractCoordinates, α::Number) = α*x
Base.:(\)(α::Number, x::ScaledCoordinates) = (x.step/α)*x.coords
Base.:(\)(α::Number, x::AbstractCoordinates) = inv(α)*x
Base.:(/)(x::AbstractCoordinates, α::Number) = α\x

Base.broadcasted(::typeof(*), α::Number, x::AbstractCoordinates) = α*x
Base.broadcasted(::typeof(*), x::AbstractCoordinates, α::Number) = α*x
Base.broadcasted(::typeof(\), α::Number, x::AbstractCoordinates) = α\x
Base.broadcasted(::typeof(/), x::AbstractCoordinates, α::Number) = α\x

for type in (:Coordinates, :RolledCoordinates, :ScaledCoordinates)
    if type === :ScaledCoordinates
        @eval begin
            convert_eltype(::Type{T}, obj::$type{T}) where {T} = obj
            convert_eltype(::Type{T}, obj::$type) where {T} =
                (as(T, obj.step)*obj.coords) :: AbstractVector{T}
        end
    else
        @eval begin
            convert_eltype(::Type{Int}, obj::$type) = obj
            convert_eltype(::Type{T}, obj::$type) where {T} = as_eltype(T, obj)
        end
    end
end

# Implement abstract array API for coordinates.
Base.length(x::AbstractCoordinates) = length(indices(x))
Base.axes(x::AbstractCoordinates) = (indices(x),)
Base.size(x::AbstractCoordinates) = (length(x),)
Base.IndexStyle(::Type{<:AbstractCoordinates}) = IndexLinear()
Base.firstindex(x::AbstractCoordinates) = first(indices(x))
Base.lastindex(x::AbstractCoordinates) = last(indices(x))
@inline function Base.getindex(x::AbstractCoordinates, i::Int)
    @boundscheck checkbounds(x, i)
    return unsafe_getindex(x, i)
end
@inline function Base.getindex(x::ScaledCoordinates, i::Int)
    @boundscheck checkbounds(x, i)
    return x.step*unsafe_getindex(x.coords, i)
end
Base.step(x::AbstractCoordinates) = 1 # FIXME: nearly correct
Base.step(x::ScaledCoordinates) = x.step

indices(x::AbstractCoordinates) = x.indices
indices(x::ScaledCoordinates) = indices(x.coords)
center(x::Coordinates) = x.center

unsafe_getindex(x::Coordinates, i::Int) = i - center(x)
unsafe_getindex(x::RolledCoordinates, i::Int) =
    i - ifelse(i ≤ x.half, firstindex(x), lastindex(x) + 1)

Base.iterate(x::AbstractCoordinates, i::Int = firstindex(x)) =
    i > lastindex(x) ? nothing : (unsafe_getindex(x, i), i + 1)

Base.minimum(x::Coordinates) = unsafe_getindex(x, lastindex(x))
Base.maximum(x::Coordinates) = unsafe_getindex(x, firstindex(x))

Base.minimum(x::RolledCoordinates) = unsafe_getindex(x, x.half + 1)
Base.maximum(x::RolledCoordinates) = unsafe_getindex(x, x.half)

Base.copy(x::AbstractCoordinates) = x

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
    FourierOptics.standardize(T, val)
    FourierOptics.standardize(obj, val)

yield value `val` with standard units and floating-point type `T`. First
argument may also be an object `obj` (e.g., a field) or a type of such an
object to use the floating-point type suitable for this object.

"""
standardize(obj::AbstractArray, val) = standardize(typeof(obj), val)
standardize(::Type{<:Field{T}}, val) where {T} = standardize(T, val)
standardize(::Type{<:AbstractArray{T}}, x) where {T} =
    standardize(floating_point_type(T), val)

standardize(::Type{T}, val::Real) where {T<:AbstractFloat} = as(T, val)
standardize(::Type{T}, val::Complex) where {T<:AbstractFloat} = as(Complex{T}, val)
standardize(::Type{T}, val::Angle) where {T<:AbstractFloat} = as(StdAngle{T}, val)
standardize(::Type{T}, val::Length) where {T<:AbstractFloat} = as(StdLength{T}, val)
standardize(::Type{T}, val::ReciprocalLength) where {T<:AbstractFloat} =
    as(StdReciprocalLength{T}, val)

# catch errors
@noinline standardize(::Type{T}, x::X) where {T,X} =
    error("don't know how to standardize value of type `$X` to use with elements of type `$T`")

"""
    FourierOptics.standard_length([T = floating_point_type(x),] x::Length)

yields a unitless floating-point value of type `T` corresponding to the length
`x` in meters.

"""
standard_length(x::Length) = standard_length(floating_point_type(x), x)
standard_length(::Type{T}, x::Length) where {T<:AbstractFloat} = ustrip(T, m, x)

"""
    FourierOptics.standard_angle([T = floating_point_type(x),] x::Angle)

yields a unitless floating-point value of type `T` corresponding to the angle
`x` in radians.

"""
standard_angle(x::Angle) = standard_angle(floating_point_type(x), x)
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
