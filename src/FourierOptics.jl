#
# FourierOptics.jl -
#
# Tools for Fourier optics computations in Julia.
#

module FourierOptics

export
    circularmask,
    circularmask!,
    fftshiftphasor,
    fftshiftphasor!

include("utils.jl")

end # module
