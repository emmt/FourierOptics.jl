#
# FourierOptics.jl -
#
# Tools for Fourier optics computations in Julia.
#

module FourierOptics

export
    CoordinateTransform,
    circularmask,
    circularmask!,
    fftshiftphasor,
    fftshiftphasor!

include("coords.jl")
import .CoordinateTransforms:
    CoordinateTransform,
    compose,
    jacobian,
    origin

include("utils.jl")

end # module
