"""
    FourierOptics.insert_lens(Fin::Field, args...; kwds...) -> Fout

out-of-place version of [`FourierOptics.insert_lens!`](@ref).

"""
insert_lens(F::Field, args...; kwds...) =
    insert_lens!(copy(F), as(Meters{T}, fl),args...; kwds...)

"""
    FourierOptics.insert_lens!(F::Field, fl::Length) -> F

inserts a lens of focal length `fl` at the current position of the field
`F`.

This method computes the phase change cause by a perfect lens that has a focal
length specified by the user. A positive focal length corresponds to a convex
lens or concave mirror; a negative length corresponds to a concave lens or
convex mirror. This method updates the new beam waist position.

See [`FourierOptics.insert_lens`](@ref) for an out-of-place version.

"""
insert_lens!(F::Field{T}, lens_fl::Length, args...; kwds...) where {T} =
    # Make sure focal length has correct units and precision.
    insert_lens!(F, as(Meters{T}, lens_fl),args...; kwds...)

function insert_lens!(F::Field{T}, lens_fl::Meters{T}) where {T}
    # Originally written by John Krist (JPL), February 2005. Adapted for Julia
    # by Éric Thiébaut, November 2023.
    #
    is_verbose(F) && println("Inserting lens")
    Δx         = get_grid_step(F)
    z          = get_surface_position(F)
    z_w0       = get_minimum_waist_position(F)
    λ          = get_wavelength(F)
    d_Rayleigh = get_Rayleigh_distance(F)
    w          = get_waist_radius_at_surface(F)

    # Determine the radius of the Gaussian beam before and after the lens. If
    # the radius is finite, then the beam is spherical; otherwise, the radius
    # is infinite and the beam is planar.
    if iszero(z - z_w0)
        # The lens is at the focus or at the entrance pupil.
        gR_beam_old = infinity(Meters{T}) # input beam is planar
        gR_beam_new = -lens_fl            # output bean is spherical
    else
        # The lens is not at the focus or at the entrance pupil.
        gR_beam_old = (z - z_w0) + d_Rayleigh^2/(z - z_w0) :: Meters{T}
        if gR_beam_old != lens_fl
            # Output beam is spherical.
            gR_beam_new = inv(inv(gR_beam_old) - inv(lens_fl)) :: Meters{T}
        else
            # Output beam is planar.
            gR_beam_new = infinity(Meters{T})
        end
    end
    is_verbose(F) && println(
        "LENS: Gaussian beam radius: input = ", gR_beam_old, ", output = ", gR_beam_new)

    # Determine current beam radius and type (near/far field).
    beam_type_old = get_field_type(F) # FIXME only use PLANAR and SPHERICAL
    if beam_type_old === :NEAR
        # Reference surface is planar.
        R_beam_old = zero(Meters{T})
    else
        # Reference surface is spherical.
        R_beam_old = z - z_w0
    end

    # Update minimal waist position and radius and recompute new Rayleigh
    # distance from focus.
    if isfinite(gR_beam_new)
        # Output beam is spherical.
        α = λ*gR_beam_new
        β = π*w^2
        z_w0 = z - gR_beam_new/(one(T) + (α/β)^2)
        w0 = w/sqrt(one(T) + (β/α)^2)
    else
        # Output beam is planar.
        z_w0 = z
        w0 = w
    end
    F.z_w0 = z_w0 # FIXME use accessor
    F.w0 = w0     # FIXME use accessor
    d_Rayleigh = get_Rayleigh_distance(F)

    # If new Rayleigh distance from focus is currently inside this,
    # then output beam is planar.
    if abs(z_w0 - z) < get_Rayleigh_factor(F)*d_Rayleigh
        beam_type_new = :NEAR
        R_beam_new = zero(Meters{T})
    else
        beam_type_new = :FAR
        R_beam_new = z - z_w0
    end


    # -- apply phase changes as needed, but don't apply if the phase is going to be
    # -- similarly altered during propagation

    if is_verbose(F)
        println("  LENS: propagator type = ", lowercase(string(beam_type_old)),
                "-to-", lowercase(string(beam_type_new)))
        println("  LENS: R_beam_old = ", R_beam_old, ", R_beam = ", R_beam,
                ",  lens_fl = ", lens_fl)
        println("  LENS: Beam diameter at lens = ", 2w)
    end

    if beam_type_old === :NEAR
        if beam_type_new === :NEAR
            lens_phase = inv(lens_fl)
        else # beam_type_new === :FAR
            lens_phase = inv(lens_fl) + inv(R_beam_new)
        end
    else # beam_type_old === :FAR
        if beam_type_new === :NEAR
            lens_phase = inv(lens_fl) - inv(R_beam_old)
        else # beam_type_new === :FAR
            if iszero(R_beam_old)
                lens_phase = inv(lens_fl) + inv(R_beam_new)
            elseif iszero(R_beam_new)
                lens_phase = inv(lens_fl) - inv(R_beam_old)
            else
                lens_phase = inv(lens_fl) - inv(R_beam_old) + inv(R_beam_new)
                if is_verbose(F)
                    println("  LENS: 1/lens_fl = ", inv(lens_fl))
                    println("  LENS: 1/R_beam_old  = ", inv(R_beam_old))
                    println("  LENS: 1/R_beam_new = ", inv(R_beam_new))
                    println("  LENS: lens_phase = ", lens_phase)
                end
            end
        end
    end
    _add_quadratic_phase!(F, -(2π/λ)*Δx^2*(lens_phase/2), Frequency(F))
    _set_surface_type!(F, beam_type_new === :NEAR ? :PLANAR : :SPHERICAL)
    is_verbose(F) && println("  LENS: Rayleigh distance = ", d_Rayleigh)
    return F
end
