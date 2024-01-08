module FourierOptics

export
    # Exports from this package.
    Angle,
    Field,
    Coordinates, RolledCoordinates,
    rectangular_aperture, rectangular_obscuration,
    circular_aperture, circular_obscuration,
    polygonal_aperture, polygonal_obscuration,
    propagate, propagate!,
    forge_mask, forge_mask!,
    insert_lens, insert_lens!,
    insert_mask, insert_mask!,
    pupil,
    split_beam, split_beam!,
    intensity, total_intensity, amplitude, true_amplitude, curvature,
    refractive_index, wavenumber, wavelength_in_medium, wavelength_in_vacuum,
    grid_size, grid_step,

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
include("pupils.jl")
include("lenses.jl")

end
