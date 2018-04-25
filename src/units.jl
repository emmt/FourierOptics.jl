#
# units.jl -
#
# Define constants for physical units.  All physical units (only lengths and
# angles) are assumed to be in SI units (hence lengths are in meters, angles
# are in radians, etc.).
#

"""
Usage:

``julia
using FourierOptics.Units
```

provides short names to be used as units.  For instance `λ = 500nm` or `λ =
0.5µm` can be used to define the wavelength, or `f = 15.2mm` to define a focal
length.

"""
module Units

export
    m, cm, mm, µm, nm,
    rd, deg, arcdeg, arcmin, arcsec, marcsec, µarcsec, mas, µas

# Lengths are given in SI units (meters).
const m  = 1.0
const cm = 1E-2
const mm = 1E-3
const µm = 1E-6
const nm = 1E-9

# Angles are given in SI units (radians).
const rd = 1.0
const arcdeg = π/180
const arcmin = π/(180*60)
const arcsec = π/(180*60*60)
const marcsec = 1E-3arcsec
const µarcsec = 1E-6arcsec
const deg = arcdeg
const mas = marcsec
const µas = µarcsec

end
