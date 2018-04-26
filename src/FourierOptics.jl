#
# FourierOptics.jl -
#
# Tools for Fourier optics computations in Julia.
#

__precompile__(true)

module FourierOptics

export
    CoordinateTransform,
    Region,
    center,
    circularmask!,
    circularmask,
    fftshiftphasor!,
    fftshiftphasor,
    recenter

using LazyAlgebra
import LazyAlgebra: Complexes, Reals

using AbstractFFTs
import FFTW

# Declaration of functions to be extended by sub-modules.
function center end
function recenter end

include("units.jl")
using .Units

include("coords.jl")
import .CoordinateTransforms:
    CoordinateTransform,
    compose,
    jacobian,
    origin

include("regions.jl")
using .Regions

include("types.jl")
include("utils.jl")
include("lenses.jl")

end # module
