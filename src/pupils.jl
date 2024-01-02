"""
    FourierOptics.pupil([T=$Float,] id; kwds...) -> mask

yields a composite geometrical object describing the entrance pupil of the
instrument identified by `id`. Optional argument `T` is the floating-point type
of the result.

Instrument identifier `id` may be a symbolic name or a `Val` singleton
parameterized by this symbolic name.

Keyword `orient` can be used to specify a rotation angle of the pupil.

Keyword `center` can be used to specify a central position of the pupil other
than the origin.

Keyword `scale` can be used to specify a scaling factor for the pupil.

Some instruments have additional keywords. For instance, keyword `focus` can be
used for a VLT's Unit Telescope (UT) to specify which focus to consider.

"""
pupil(id::Symbol, args...; kwds...) = pupil(Val(id), args...; kwds...)
pupil(::Type{T}, id::Symbol, args...; kwds...) where {T} = pupil(T, Val(id), args...; kwds...)
pupil(id::Val, args...; kwds...) = pupil(Float, id, args...; kwds...)

"""
    FourierOptics.segment(p0, p1, w) -> rect

yields a polygonal object corresponding to a rectangular segment of thickness
`w` joining points `p0` and `p1`. The result is similar to cap style *butt* in
computer graphics. Schematically, for an horizontal segment:

     ┌───────────────────┐  ┬
     │                   │  │
     ┼ p0             p1 ┼  │w
     │                   │  │
     └───────────────────┘  ┴

"""
function segment(p0::Point, p1::Point, w)
    u = p1 - p0 # vector along "length" side
    v = (w/(2*norm(u)))*Point(-u.y, u.x) # vector along "half-width" side
    return Polygon(p0 - v, p1 - v, p1 + v, p0 + v)
end

# Parameters for the VLT Unit Telescopes (source: https://www.eso.org/sci/facilities/develop/documents/VLT-SPE-ESO-10000-2723_is1.pdf):
#
# VLT Cassegrain focus optical parameters:
#   pupil diameter = 8.115 m
#   focal ratio = 13.4106
#   focal length = 108.827 mm
#   object field of view = 15 arcmin (total)
#                        = 2.68 arcmin (unvignetted)
#   image field of view = 474.4 mm (total)
#                       =  85.0 mm (unvignetted)
#   image scale = 528 μm/arcsec
#   radius of image curvature = 1981.4 mm (concave towards M2)
#
# VLT Nasmyth focus optical parameters:
#   pupil diameter = 8.000 m
#   focal ratio = 15.0
#   focal length = 120.000 m
#   object field of view = 30 arcmin (total)
#                        = 7.15 arcmin (unvignetted)
#   image field of view = 1043.8 mm (total)
#                       =  249.6 mm (unvignetted)
#   image scale = 582 μm/arcsec
#   radius of image curvature = 2089.6 mm (concave towards M2)
#
function pupil(::Type{T}, ::Val{:VLT};
               focus::Symbol = :Nasmyth,
               center::Union{Point,TwoTuple} = Point(0mm,0mm),
               orient::Union{Angle,Real} = zero(T),
               scale::Real = one(T)) where {T<:AbstractFloat}
    focus ∈ (:Cassegrain, :Nasmyth) || throw(ArgumentError("invalid focus"))
    L = typeof(one(T)*mm)  # length units
    A = typeof(one(T)*°)   # angle units
    center = Point{L}(center)
    orient = convert(A, orient)

    # Diamneters of pupil aperture and obscuration.
    outer_diameter = convert(L, (focus === :Cassegrain ? 8.115m : 8.000m))
    inner_diameter = convert(L, 1.116m)

    # Spider vannes of a VLT's Unit Telescope (UT) consist in 4 segments
    # forming 2 V's with an angle of 101° between the branches of each V. The
    # submits of the V's are attached at the inner radius of the mount, where
    # sit the M2 mirror, at angles 45° and 45° + 180° = 225°. The endpoints of
    # the V's are attached at the outer radius of the mount at angles 0°, 90°,
    # 180°, and 270°.
    #
    # S = (r/√2, r/√2)        # submit at inner radius `r` and 45°
    # A = (R,0)               # 1st endpoint at outer radius `R` and 0°
    # SA: A + t*(cosβ,-sinβ)  # equation of SA line
    # β = 101°/2 - 45° = 5.5° # SA angle with Ox
    #
    # S along SA is such that
    #
    #     Ax + t*cosβ = Ay - t*sinβ
    #     => t = -R/(sinβ + cosβ)
    #     => Sx = Sy = R*sinβ/(sinβ + cosβ)
    #     => r = √2*|sinβ/(sinβ + cosβ)|*R
    #
    spider_width = convert(L, 40mm) # thickness of spider vannes
    spider_outer_radius = convert(L, 4.2197m)
    spider_V_angle = 101°
    β = convert(A, spider_V_angle)/2 - 45°
    sinβ, cosβ = sincos(β)
    spider_inner_radius = abs(sqrt(as(T, 2))*sinβ/(sinβ + cosβ))*spider_outer_radius

    # Submits of the V's:
    S1 = center + Point{L}(r = spider_inner_radius, θ = orient + 45°)
    S2 = center + Point{L}(r = spider_inner_radius, θ = orient + (45° + 180°))

    # Endpoints of the V's:
    E1 = center + Point{L}(r = spider_outer_radius, θ = orient +   0°)
    E2 = center + Point{L}(r = spider_outer_radius, θ = orient +  90°)
    E3 = center + Point{L}(r = spider_outer_radius, θ = orient + 180°)
    E4 = center + Point{L}(r = spider_outer_radius, θ = orient + 270°)

    # Spider vannes are segments:
    V1 = segment(S1, E1, spider_width)
    V2 = segment(S1, E2, spider_width)
    V3 = segment(S2, E3, spider_width)
    V4 = segment(S2, E4, spider_width)

    return (
        # FIXME: specify type parameter
        CircularAperture(; center, diameter = outer_diameter), # aperture
        CircularObscuration(; center, diameter = inner_diameter), # central obscuration
        PolygonalObscuration(V1),
        PolygonalObscuration(V2),
        PolygonalObscuration(V3),
        PolygonalObscuration(V4))
end
