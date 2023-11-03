const default_to_plane = false

"""
    FourierOptics.propagate(Fin::Field, dz::Length; to_plane=$default_to_plane) -> Fout

propagates the input field `Fin` over a distance `dz` and returns the resulting
field `Fout`. This propagation method automatically chooses among near-to-near,
near-to-far, far-to-near, and far-to-far field propagators.

Keyword `to_plane` specifies whether to force propagation to a planar reference
surface.

"""
propagate(F::Field, args...; kwds...) = propagate!(copy(F), args...; kwds...)

"""
    FourierOptics.propagate!(F::Field, dz::Length; to_plane=$default_to_plane) -> F

propagates in-place the field `F` over a distance `dz`. This propagation method
automatically chooses among near-to-near, near-to-far, far-to-near, and
far-to-far field propagators.

Keyword `to_plane` specifies whether to force propagation to a planar reference
surface.

"""
function propagate!(F::Field, dz::Length;
                    to_plane::Bool = default_to_plane)
    # Code originally written in IDL by John Krist (JPL), February 2005.
    # Adapted to Julia by Éric Thiébaut, November 2023.
    field_type_old = get_field_type(F)
    field_type_new = to_plane ? :NEAR : get_field_type(F, dz)
    is_verbose(F) && println(
        "Propagating by ", dz, " with a ", lowercase(string(field_type_old)),
        "-to-", lowercase(string(field_type_new)), " propagator.")
    z_w0 = get_minimal_waist_position(F)
    z_old = get_surface_position(F)
    z_new = z_old + dz
    if field_type_old === :NEAR
        # Current reference surface is planar, previous propagation was to the
        # near field.
        if field_type_new === :NEAR
            # Apply near-to-near field propagator.
            _apply_ptp!(F, dz)
        else
            # Apply near-to-far field propagator. First from current position
            # to the waist, then from the waist to the new position.
            _apply_ptp!(F, z_w0 - z_old)
            _apply_wts!(F, z_new - z_w0)
        end
    else
        # Current reference surface is spherical, previous propagation was to
        # the far field.
        if field_type_new === :NEAR
            # Apply far-to-near field propagator.
            _apply_stw!(F, z_w0 - z_old)
            _apply_ptp!(F, z_new - z_w0)
        else
            # Apply far-to-far field propagator.
            _apply_stw!(F, z_w0 - z_old)
            _apply_wts!(F, z_new - z_w0)
        end
    end
    @assert get_surface_type(F) === (field_type_new === :NEAR ? :PLANAR : :SPHERICAL)
    @assert get_surface_position(F) ≈ z_new
    is_verbose(F) && println("Total intensity = ", get_total_intensity(F))
    return F
end

"""
    FourierOptics._apply_ptp!(F::Field, dz::Length) -> F

propagates in-place a planar input wavefront `F` over distance `dz`, keeping it
planar.

This private method is used to propagate a planar input wavefront over some
distance to produce an planar output wavefront. This occurs when both the start
and end point are both within the Rayleigh distance of focus.

"""
function _apply_ptp!(F::Field{T}, dz::Length) where {T}
    #
    # Near-to-near field propagation amounts to convolving by a kernel whose
    # Fourier transform is a phasor with a quadratic phase given by:
    #
    #     exp(-i*π*λ*dz*(u₁^2 + u₂^2))
    #
    # where
    #
    #     (u₁,u₂) = (k₁/(n*Δx),k₂/(n*Δx))
    #
    # is the 2-dimensional spatial frequency and `(k₁,k₂)` the 2-dimensional
    # discrete spatial frequency. As an optimiation, expanding the exponential
    # yields:
    #
    #     exp(-i*π*λ*dz*ρ^2) = exp(i*α*k₁^2)*exp(i*α*k₂^2)
    #
    # with `α = -π*λ*dz/(n*Δx)^2`.
    #
    # Originally written by John Krist (JPL), February 2005. Updated by John
    # Krist, (JPL) Oct 2013. Adapted to Julia by Éric Thiébaut, November 2023.
    #
    @assert F.reference_surface === :PLANAR
    is_verbose(F) && println("PTP: dz = ", dz)
    if !iszero(dz) # FIXME use some small threshold
        dz = as(Meters{T}, dz) # to avoid multiple conversions
        λ = get_wavelength(F)
        n = get_grid_size(F)
        Δx = get_grid_step(F)
        _apply_fft!(F, -1) # apply forward orthogonal FFT
        _apply_quadratic_phase!(F, -π*λ*dz/(n*Δx)^2, Frequency(F))
        _apply_fft!(F, +1) # apply backward orthogonal FFT
        if propagate_phase_offset(F); multiply!(F, exp_i(2π*dz/λ)); end
        _increment_position!(F, dz)
    end
    is_verbose(F) && println("PTP: z = ", get_position(F), "  dx = ", get_grid_step(F))
    return F
end

"""
    FourierOptics._apply_stw!(F::Field, dz::Length) -> F

applies in-place the *spherical-to-waist* propagator to the field `F` over the
distance `dz`.

This private method propagates from a spherical reference surface that is
outside the Rayleigh limit from focus to a planar one that is inside.

"""
function _apply_stw!(F::Field{T}, dz::Length) where {T}
    # Code originally written in IDL by by John Krist (JPL), February 2005.
    # Adapted to Julia by Éric Thiébaut, November 2023.
    @assert get_surface_type(F) === :SPHERICAL
    @assert !iszero(dz)
    is_verbose(F) && println("STW: dz = ", dz)
    dz = as(Meters{T}, dz) # to avoid multiple conversions
    λ = get_wavelength(F)
    n = get_grid_size(F)
    Δx = get_grid_step(F, dz) # grid step *after* propagation
    _apply_fft!(F, (dz ≥ zero(dz) ? -1 : +1))
    _apply_quadratic_phase!(F, (π/(λ*dz))*Δx^2, Frequency(F))
    if propagate_phase_offset(F); multiply!(F, exp_i(2π*dz/λ)); end
    _set_grid_step!(F, Δx)
    _increment_position!(F, dz)
    _set_surface_type(F, :PLANAR)
    is_verbose(F) && println("STW: z = ", get_position(F), "  dx = ", get_grid_step(F))
    return F
end

"""
    FourierOptics._apply_wts!(F::Field, dz::Length) -> F

applies in-place the *waist-to-spherical* propagator to the field `F` over the
distance `dz`. This operation is done in-place.

This private method propagates from a planar reference surface that is inside
the Rayleigh distance from focus to a spherical reference surface that is
outside.

"""
function _apply_wts!(F::Field{T}, dz::Length) where {T}
    # Code originally written in IDL by by John Krist (JPL), February 2005.
    # Adapted to Julia by Éric Thiébaut, November 2023.
    @assert get_surface_type(F) === :PLANAR
    @assert !iszero(dz)
    is_verbose(F) && println("WTS: dz = ", dz)
    dz = as(Meters{T}, dz) # to avoid multiple conversions
    λ = get_wavelength(F)
    n = get_grid_size(F)
    Δx = get_grid_step(F) # grid step *before* propagation
    _apply_quadratic_phase!(F, (π/(λ*dz))*Δx^2, Frequency(F))
    _apply_fft!(F, (dz ≥ zero(dz) ? -1 : +1))
    if propagate_phase_offset(F); multiply!(F, exp_i(2π*dz/λ)); end
    Δx = get_grid_step(F, dz) # grid step *after* propagation
    _set_grid_step!(F, Δx)
    _increment_position!(F, dz)
    _set_surface_type(F, :SPHERICAL)
    is_verbose(F) && println("WTS: z = ", get_position(F), "  dx = ", get_grid_step(F))
    return F
end

# FIXME for func in (_apply_ptp!, _apply_stw!, _apply_wts!)
# FIXME     @eval $func(F::Field{T}, dz::Length) where {T} = $func(F, as(Meters{T}, dz))
# FIXME end

"""
    FourierOptics._apply_quadratic_phase!(F::Field, α, x) -> F

multiplies in-place the complex amplitude at every 2-dimensional index
`(j₁,j₂)` of the field `F` by:

    exp(i*α*(x[j₁]^2 + x[j₂]^2)) = q[j₁]*q[j₂]

with:

    q[j] = exp(i*α*x[j]^2)

Argument `α` is a real factor and argument `x` is the vector of coordinates
(assumed to be the same for the 2 axes of the field).

"""
function _apply_quadratic_phase!(F::Field{T},
                                 α::Dimensionless{Real},
                                 x::AbstractVector{<:Real}) where {T}
    # Convert `α` to a bare real of the correct floating-point type.
    return _apply_quadratic_phase!(F, as(T, α), x)
end

# FIXME: In the future, use x and y.
function _apply_quadratic_phase!(F::Field{T}, α::T, x::AbstractVector{<:Real}) where {T}
    ampl = get_amplitude(F)
    J₁, J₂ = axes(ampl)
    J = axes(x,1)
    @assert J == J₁ == J₂
    q = similar(ampl, Complex{T}, axes(x))
    @assert axes(q) == (J,)
    @inbounds @simd for j in J
        xⱼ = as(T, x[j])
        q[j] = exp_i(α*xⱼ^2)
    end
    @inbounds for j₂ in J₂
        @simd for j₁ in J₁
            ampl[j₁,j₂] *= q[j₁]*q[j₂]
        end
    end
    return F
end

"""
    FourierOptics._apply_fft!(F::Field, dir) -> F

applies the orthogonal Fast Fourier Transform (FFT) to the complex amplitude in
the field `F`. Argument `dir` specifies the sign of the argument in the complex
exponential of the transform: if `dir < 0`, the direct orthogonal FFT is
applied; otherwise, `dir > 0` and the inverse orthogonal FFT is applied.

"""
function _apply_fft!(F::Field{T}, dir::Integer) where {T}
    # Apply in-place complex-complex forward/backward FFT transform, then multiply
    # by 1/sqrt(length(F)).
    ampl = get_amplitude(F)
    if dir < zero(dir)
        mul!(ampl, F.forward_plan, ampl)
    elseif dir > zero(dir)
        mul!(ampl, F.backward_plan, ampl)
    else
        throw(ArgumentError("FFT direction has undefinite sign"))
    end
    α = inv(sqrt(as(T, length(ampl)))) # α = 1/√(n₁*n₂)
    multiply!(ampl, α)
    return F
end
