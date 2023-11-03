module FourierOptics

export
    # Exports from this package.
    Angle,
    Field,
    propagate,
    propagate!,

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

using TypeUtils
using LinearAlgebra

include("types.jl")
include("utils.jl")
include("fields.jl")
include("propagate.jl")
include("lenses.jl")

end
