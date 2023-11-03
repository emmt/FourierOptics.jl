#
# types.jl -
#
# Types used in Fourier optics.
#

"""

Positions in planes perpendicular to the direction of propagation are specified
as tuples of two reals `(x,y)`.

"""
const Position = NTuple{2,Real}

abstract type AbstractLens end

struct Lens <: AbstractLens
    f::Float64   # focal length of lens
    D::Float64   # diameter of lens
    x0::Float64  # abscissa of center of lens
    y0::Float64  # ordinate of center of lens
end

struct IndexBox
    I::UnitRange{Int64}            # range of first index
    J::UnitRange{Int64}            # range of second index
end

struct LensOperator{Tl<:AbstractLens,
                    Tf<:Base.DFT.Plan{Complex{Float64}}} <: LinearMapping
    lens::Tl                       # lens parameters
    lambda::Float64                # wavelength
    Rinp::Region                   # pupil plane
    Qinp::Matrix{Complex{Cdouble}} # complex transmission in pupil plane
    Binp::IndexBox                 # part to consider the input region
    Rout::Region                   # focal plane
    Qout::Matrix{Complex{Cdouble}} # complex modulation in focal plane
    FFT::Tf                        # FFT plan
    ws::Matrix{Complex{Cdouble}}   # workspace for the FFT
end

struct ComplexAmplitude{T<:Complexes,M<:AbstractMatrix{<:Complexes}} <: AbstractMatrix{T}
    amp::M          # the complex amplitude
    R::Region       # physical plane
    lambda::Float64 # wavelength (in meters)
end
