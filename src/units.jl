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
    showlength,
    showposition,
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

"""
```julia
showlength(io, x)
```

prints length `x` (assumed to be in SI units, that is meters) to `io`
in a human readable way.

"""
function showlength(io::IO, x::Real)
    a = abs(x)
    if 1mm ≤ a < 1m
        print(io, @sprintf("%.6gmm", x/1mm))
    elseif 1µm ≤ a < 1mm
        print(io, @sprintf("%.6gµm", x/1µm))
    elseif 1nm ≤ a < 1µm
        print(io, @sprintf("%.6gnm", x/1nm))
    else
        print(io, @sprintf("%.6gm", x/1m))
    end
end

showlength(io::IO, x::Union{Rational,Irrational}) =
    showlength(io, convert(Float64, x))

"""
```julia
showposition(io, (x,y))
```

prints position `(x,y)` (both coordinates are assumed to be in SI units, that
is meters) to `io` in a human readable way.

"""
function showposition(io::IO, c::NTuple{2,Union{Real,Rational,Irrational}})
    print(io, "(")
    showlength(io, c[1])
    print(io, ",")
    showlength(io, c[2])
    print(io, ")")
end

end
