"""
    FourierOptics.Float

is the default floating-point type for computations.

"""
const Float = Float64

"""
    FourierOptics.Angle

is the union of quantities suitable to specifiy an angle.

"""
const Angle = Union{typeof(°),typeof(rad)}

"""
    FourierOptics.Meters{F}

is the type of a length expressed in meters and with floating-point type `F`.

"""
const Meters{F<:AbstractFloat} = Unitful.Quantity{F, Unitful.𝐋, Unitful.FreeUnits{(Unitful.Unit{:Meter,Unitful.𝐋}(0, 1//1),), Unitful.𝐋, nothing}}

"""
    FourierOptics.Dimensionless{Real}
    FourierOptics.Dimensionless{Complex}

are unions of types suitable for dimensionless real and complex values. These
aliases are introduced to cope with dimensionless quantities that may result
from expressions involving `Unitful` quantities and that have not yet been
converted to a bare numerical type.

"""
const Dimensionless{T} = Union{T,Unitful.Quantity{<:T,Unitful.NoDims}}

# Default settings.
const default_fftw_flags = FFTW.MEASURE
const default_fftw_timelimit = FFTW.NO_TIMELIMIT
const default_Rayleigh_factor = 1
const default_phase_offset = false
const default_verbosity = 0

"""
    FourierOptics.Field(ampl; wavelength, grid_step, kwds...)

builds a new field with complex amplitude stored in array `ampl`. The
wavelength and grid sampling step size must be specified by the keywords
`wavelength` and `grid_step`.

Other optional keywords are:

- `fftw_flags`, and `fftw_timelimit` the flags and time-limit for creating
  the plans for FFTW.

- `Rayleigh_factor` is a factor for adjusting the Rayleigh distance to
  determine the boundary of the near and far fields. By default,
  `Rayleigh_factor = $default_Rayleigh_factor`.

- `phase_offset` specifies whether to propagate the global phase offset of the
  field. This may be important when modeling the separate arms of an
  interferometer with a path difference between the two. By default,
  `phase_offset = $default_phase_offset`.

- `verbosity` speicifies the level of verbosity. By default, `verbosity =
  $default_verbosity`.

"""
mutable struct Field{T<:AbstractFloat,               # floating-point type
                     A<:AbstractMatrix{Complex{T}},  # complex amplitude
                     L,                              # linear index style?
                     F<:Plan,                        # forward FFT plan
                     B<:Plan,                        # backward FFT plan
                     } <: AbstractMatrix{Complex{T}}
    ampl::A               # complex amplitude
    forward_plan::F       # forward FFT plan
    backward_plan::B      # backward FFT plan
    Rayleigh_factor::T    # factor to determine the boundary of the near and far fields
    lambda::Meters{T}     # wavelength
    dx::Meters{T}         # sampling per pixel
    w0::Meters{T}         # minimal radius of beam waist
    z::Meters{T}          # position along propagation axis of surface
    z_w0::Meters{T}       # position along propagation axis of minimal beam waist
    surface_type::Symbol  # type of reference surface: `:PLANAR` or `:SPHERICAL`
    verbosity::Int        # verbosity level
    phase_offset::Bool    # propagate global phase offset?

    # Inner constructor.
    function Field(ampl::A,
                   forward_plan::F,
                   backward_plan::B;
                   wavelength::Length,
                   grid_step::Length,
                   verbosity::Integer = 0,
                   Rayleigh_factor::Real = default_Rayleigh_factor,
                   phase_offset::Bool = default_phase_offset,
                   _preserve_amplitude::Bool = false, # private option
                   ) where {T<:AbstractFloat,
                            A<:AbstractMatrix{Complex{T}},
                            F<:Plan, B<:Plan}
        # Check arguments.
        n = size(ampl, 1)
        size(ampl, 2) == n || error("expecting a square complex amplitude array")
        Base.has_offset_axes(ampl) && error("complex amplitude array must have 1-based indices")
        eltype(forward_plan) == Complex{T} || error("incompatible element type for FFT forward plan")
        axes(forward_plan) == axes(ampl) || error("incompatible indices for FFT forward plan")
        eltype(backward_plan) == Complex{T} || error("incompatible element type for FFT backward plan")
        axes(backward_plan) == axes(ampl) || error("incompatible indices for FFT backward plan")

        # Assume the beam diameter is one-half of the surface size, hence its
        # radius is one-fourth of the surface size.
        w0 = n*grid_step/4

        # Fill amplitude with ones.
        _preserve_amplitude || fill!(ampl, one(eltype(ampl)))

        # Create and instantiate structure.
        L = IndexStyle(ampl) === IndexLinear()
        obj = new{T,A,L,F,B}()
        obj.ampl = ampl
        obj.forward_plan = forward_plan
        obj.backward_plan = backward_plan
        obj.Rayleigh_factor = check_Rayleigh_factor(Rayleigh_factor)
        obj.lambda = check_wavelength(wavelength)
        obj.dx = check_grid_step(grid_step)
        obj.w0 = w0
        obj.z = 0.0m
        obj.z_w0 = 0.0m
        obj.surface_type = :PLANAR
        obj.verbosity = verbosity
        obj.phase_offset = phase_offset
        return obj
    end
end

"""
    FourierOptics.AbstractCoordinates{T} <: AbstractVector{T}

is the super-type of vector-like objects whose elements give the coordinates
along an array axis. Usually such objects are built so as to have the same
indices as the corresponding array axis.

"""
abstract type AbstractCoordinates{T}  <: AbstractVector{T} end

"""
    FourierOptics.Coordinates(inds[, cen])

yields a lightweight vector-like object whose elements give the coordinates
along an array axis with indices `inds`. If argument `inds` is an integer, it
is assumed to be the length of an array axis with 1-based indices. Optional
argument `cen` is the *central* index whose coordinate is zero. If `cen` is
unspecified, it is assumed to be at the geometrical center of the axis
(following the same conventions as in `fftshift` and `ifftshift`).

"""
struct Coordinates{I<:AbstractUnitRange{Int}} <: AbstractCoordinates{Int}
    indices::I  # range of valid indices
    center::Int # index of center
    function Coordinates(indices::I,
                        center::Integer = default_coordinates_center(indices),
                        ) where {I<:AbstractUnitRange{Int}}
        return new{I}(indices, center)
    end
end

"""
    FourierOptics.Frequencies(inds)

yields a lightweight vector-like object whose elements give the discrete FFT
frequencies along an array axis with indices `inds`. If argument `inds` is an
integer, it is assumed to be the length of an array axis with 1-based indices.

The same conventions as in `fftshift` and `ifftshift` are assumed. In other
words, if `n` is the length of the considered axis, then the discrete
frequencies are:

    [0, 1, ..., n÷2,     -(n÷2), 1 - (n÷2), ..., -2, -1]    if `n` is odd

    [0, 1, ..., n÷2 - 1, -(n÷2), 1 - (n÷2), ..., -2, -1]    if `n` is even

"""
struct Frequencies{I<:AbstractUnitRange{Int}} <: AbstractCoordinates{Int}
    indices::I # range of valid indices
    half::Int  # index of positive Nyquist frequency
    function Frequencies(indices::I) where {I<:AbstractUnitRange{Int}}
        n = length(indices)
        half = first(indices) + (isodd(n) ? n÷2 : n÷2 - 1)
        return new{I}(indices, half)
    end
end
