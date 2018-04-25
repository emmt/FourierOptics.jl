#
# FourierOptics.jl -
#
# Tools for Fourier optics computations in Julia.
#

__precompile__(true)

module FourierOptics

export
    CoordinateTransform,
    circularmask,
    circularmask!,
    fftshiftphasor,
    fftshiftphasor!

include("units.jl")
using .Units

include("coords.jl")
import .CoordinateTransforms:
    CoordinateTransform,
    compose,
    jacobian,
    origin

include("utils.jl")

end # module
