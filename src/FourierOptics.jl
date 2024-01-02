module FourierOptics

export
    # Exports from this package.
    Angle,
    Field,
    Coordinates, RolledCoordinates,
    RectangularAperture, RectangularObscuration,
    CircularAperture, CircularObscuration,
    PolygonalAperture, PolygonalObscuration,
    propagate, propagate!,
    forge_mask, forge_mask!,

    # Exports from Unitful.
    Length,
    °, rad,
    km, m, cm, mm, µm, nm

using Unitful
using Unitful:
    Length,
    °, rad,
    km, m, cm, mm, µm, nm

using AbstractFFTs, FFTW
import AbstractFFTs: Plan, fftshift, ifftshift

using ArrayTools, TypeUtils, Unitless
using LinearAlgebra

import Unitless: floating_point_type
import TypeUtils: as_eltype, convert_eltype

using Base: axes1, Fix1, Fix2

include("types.jl")
include("utils.jl")
include("fields.jl")
include("propagate.jl")
include("winding.jl")
include("masks.jl")
include("lenses.jl")

end
