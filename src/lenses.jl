"""
    FourierOptics.insert_lens(Fin::Field, args...; kwds...) -> Fout

out-of-place version of [`FourierOptics.insert_lens!`](@ref).

"""
insert_lens(F::Field, args...; kwds...) =
    insert_lens!(copy(F), args...; kwds...)

"""
    FourierOptics.insert_lens!(F::Field, fl::Length) -> F

inserts a thin lens of focal length `fl` in the field `F` at its current
position. A positive focal length corresponds to a convex lens or concave
mirror; a negative length corresponds to a concave lens or convex mirror.

See [`FourierOptics.insert_lens`](@ref) for an out-of-place version.

"""
function insert_lens!(F::Field, fl::Length)
    # Simply add reciprocal of focal length to current wavefront curvature.
    F.curv -= inv(standardize(F, fl))
    return F
end
