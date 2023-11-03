#
# FourierOptics.jl -
#
# Tools for Fourier optics computations in Julia.
#

__precompile__(true)

module FourierOptics

export
    CoordinateTransform,
    Lens,
    LensOperator,
    Region,
    center,
    diameter,
    circularmask!,
    circularmask,
    fftshiftphasor!,
    fftshiftphasor,
    focal_length,
    grid2world,
    recenter,
    world2grid

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
include("basics.jl")
include("utils.jl")
include("lenses.jl")

end # module
