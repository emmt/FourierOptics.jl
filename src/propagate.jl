const default_propagation_method = :Fresnel2

"""
    FourierOptics.propagate(Fin::Field, Δz::Length; by=:$default_propagation_method, kwds...) -> Fout

propagates the input field `Fin` over a distance `Δz` and returns the resulting
field `Fout`. The input field `Fin` is left unmodified, method
[`FourierOptics.propagate!`](@ref) may be used for in-place propagation.

Keyword `by` specifies the method to perform the propagation:

* `by = :Fresnel1` to apply a single-FFT Fresnel propagation method. This
  propagation method costs 1 FTT and imposes that the output lateral sampling
  step be:

  ```julia
  δx1 = λ⋅|Δz|/(N⋅δx0)
  ```

  with `λ` the wavelength in the propagation medium, `N` the number of samples
  along a dimension of the transverse plane, and `δx0` the lateral sampling
  step before propagation. Single-FFT Fresnel propagation method is adapted for
  *far field* propagation, that is for:

  ```julia
  |Δz| ≥ min(D,N⋅δx0)⋅δx0/λ
  ```

  with `D` the diameter of the beam.

* `by = :Fresnel2` to apply a 2-step Fresnel propagation method. Keyword
  `sampling = ±δx1` can be used to specify the sampling step `δx1` after
  propagation (up to a change of sign). The value of the `sampling` keyword may
  be negative to use a negative algebraic magnification, the output sampling
  step is `δx1 = abs(sampling)`. By default, the sampling step of `F` is
  retained and the method is equivalent to an angular spectrum method under the
  paraxial approximation. The 2-step Fresnel propagation method costs 2 FTTs
  and is equivalent to the fractional Fourier transform method.

* `by = :RayleighSommerfeld` to propagate the field by an angular spectrum
  method which amounts to convolving the complex amplitude by the
  Rayleigh-Sommerfeld propagation kernel. Keyword `no_evanescent_waves`
  (`false` by default) indicates whether to explicitly filter out evanescent
  waves.

*Fresnel propagation* amounts to assuming *paraxial conditions*, that is small
diffraction angles or, equivalently, spatial frequencies much smaller than
`1/λ` with `λ` the wavelength in the propagation medium.

*Angular spectrum* method means that the propagation is performed by convolving
the complex amplitude by the propagation kernel and using the sampled Fourier
transform of the kernel, the so-called propagation transfer function. Angular
spectrum propagation costs 2 FFTs and leaves the lateral sampling step `δx`
unchanged. Angular spectrum method is adapted for *near field* propagation,
that is for:

    |Δz| ≤ min(D,N⋅δx)⋅δx/λ

with `D` the diameter of the beam, `N` the number of samples along a dimension
of the transverse plane, and `λ` the wavelength in the propagation medium.

"""
propagate(F::Field, args...; kwds...) = propagate!(copy(F), args...; kwds...)

"""
    FourierOptics.propagate!(F::Field, Δz::Length; kwds...) -> F

propagates in-place the field `F` over a distance `Δz`. See
[`FourierOptics.propagate`](@ref) for a description.

"""
propagate!(F::Field, Δz::Length; by=default_propagation_method, kwds...) =
    propagate!(by, F, Δz; kwds...)

propagate!(by::Symbol, F::Field, Δz::Length; kwds...) =
    propagate!(Val(by), F, Δz; kwds...)

# Catch errors.
@noinline propagate!(by::Val{M}, F::Field, Δz::Length; kwds...) where {M} =
    throw(ArgumentError("unknown propagation method `$M`"))

function propagate!(by::Val{:Fresnel1}, F::Field, Δz::Length)
    # Retrieve current field parameters.
    Δz = standardize(F, Δz)
    iszero(Δz) && return F
    δx = grid_step(F)
    λ = wavelength_in_medium(F)
    k = wavenumber(F)
    N = grid_size(F)

    # Propagate the complex amplitude.
    ρ = (δx^2/λ)*(curvature(F) + inv(Δz))
    println("ρ = $ρ")
    apply_quadratic_phase_factor!(F, ρ; planar=true)
    apply_fft!(F, (Δz ≥ zero(Δz) ? -1 : 1))
    F.curv = inv(Δz) # set resulting wavefront curvature

    # Update field parameters other than the complex amplitude and the
    # wavefront curvature.
    F.fact *= δx^2*exp_i(k*Δz)/(i*λ*Δz) # scale complex amplitude
    F.δx = λ*abs(Δz)/(N*δx) # update lateral sampling step
    F.z += Δz # update longitudinal position
    return F
end

function propagate!(by::Val{:Fresnel2}, F::Field, Δz::Length;
                    sampling::Length = grid_step(F))
    # Retrieve current field parameters.
    Δz = standardize(F, Δz)
    iszero(Δz) && return F
    δx1 = standardize(F, sampling)
    δx0 = grid_step(F)
    gamma = δx1/δx0 # algebraic magnification
    δx1 = abs(δx1) # output sampling step is always positive
    λ = wavelength_in_medium(F)
    k = wavenumber(F)
    N = grid_size(F)

    # Parameters of the quadratic phase factors.
    ρ1 = (δx0^2/λ)*(curvature(F) + (one(gamma) - gamma)/Δz)
    ρ2 = -λ*Δz/(gamma*N^2*δx0^2)
    ρ3 = gamma*(gamma - one(gamma))*δx0^2/(λ*Δz)

    # Propagate the complex amplitude.
    apply_quadratic_phase_factor!(F, ρ1; planar=true)
    apply_fft!(F, -1) # complex amplitude -> angular spectrum
    apply_quadratic_phase_factor!(F, ρ2)
    apply_fft!(F, +1) # angular spectrum -> complex amplitude
    F.curv = (λ/δx1^2)*ρ3 # set resulting wavefront curvature

    # Update field parameters other than the complex amplitude and the
    # wavefront curvature.
    F.fact *= exp_i(k*Δz)/(gamma*N^2) # scale complex amplitude
    F.δx = δx1 # update sampling
    F.z += Δz # update longitudinal position
    return F
end

function propagate!(by::Val{:RayleighSommerfeld}, F::Field, Δz::Length;
                    no_evanescent_waves::Bool = false)
    # Retrieve current field parameters.
    Δz = standardize(F, Δz)
    iszero(Δz) && return F
    δx = grid_step(F)
    λ = wavelength_in_medium(F)
    k = wavenumber(F)
    N = grid_size(F)
    δα = inv(N*δx) # spatial frequency sampling step

    # Parameters for computing the angular spectrum propagation transfer
    # function.
    T = floating_point_type(F)
    kΔz = as(T, k*Δz)
    λα = as(T, λ*δα)*RolledCoordinates(F) # pseudo-vector of λ*α values
    λβ = λα                               # pseudo-vector of λ*β values
    amp = amplitude(F)
    J₁, J₂ = axes(amp)
    (axes1(λα) == J₁ && axes1(λβ) == J₂) || throw(DimensionMismatch(
        "complex amplitude has invalid axes"))

    # Propagate the complex amplitude.
    apply_quadratic_phase_factor!(F) # make sure wavefront is planar
    apply_fft!(F, -1) # complex amplitude -> angular spectrum
    @inbounds for j₂ in J₂
        λ²β² = λβ[j₂]^2
        for j₁ in J₁
            λ²α² = λα[j₁]^2
            ξ = one(T) - (λ²α² + λ²β²)
            if ξ > zero(ξ)
                # Low frequency part: propagating waves.
                amp[j₁,j₂] *= exp_i(kΔz*sqrt(ξ))
            elseif ξ < zero(ξ)
                # High frequency part: evanescent waves.
                if no_evanescent_waves
                    amp[j₁,j₂] = zero(eltype(amp))
                else
                    amp[j₁,j₂] *= exp(-kΔz*sqrt(-ξ))
                end
            end
        end
    end
    apply_fft!(F, +1) # angular spectrum -> complex amplitude

    # Update field parameters other than the complex amplitude and the
    # wavefront curvature.
    F.fact *= exp_i(k*Δz)/N^2 # scale complex amplitude
    F.z += Δz # update longitudinal position
    return F
end

"""
    FourierOptics.propagate!(Val(:Mas1999), F::Field, Δz::Length) -> F

propagates the field `F` over a distance `Δz` by the propagation method of Mas
et al. (Optics Communications, vol. 164, pp. 233–245, 1999) which amounts to
performing a 2-step Fresnel propagation with a magnification automatically
chosen.

The operation is done in-place. On output, the grid sampling step is updated.

"""
function propagate!(by::Val{:Mas1999}, F::Field, Δz::Length; other::Bool=false)
    Δz = standardize(F, Δz)
    iszero(Δz) && return F
    λ = wavelength_in_medium(F)
    δx₀ = grid_step(F)
    N = grid_size(F)
    f₁ = N*δx₀^2/λ # far-/near-fields limit
    #Δx = N*δx
    #ϕ = atan(λ*Δz/Δx^2)
    η = sqrt(1 + (Δz/f₁)^2) # magnification factor
    # displacement to intermediate plane:
    Δzₘ = Δz/(other ? η - one(η) : one(η) + η)
    propagate!(:Fresnel1, F, Δzₘ)
    propagate!(:Fresnel1, F, Δz - Δzₘ)
    δx₁ = δx₀*η # output sampling step
    @assert grid_step(F) ≈ δx₁
    return F
end

"""
    FourierOptics.range(F::Field, D::Length=N⋅δx) -> Δzₘ

yields:

    Δzₘ = D⋅δx/λ

with `δx` the lateral sampling step for the field `F`, and `λ` the wavelength
in the propagation medium. The range `Δzₘ` is the boundary of the near and the
far fields for the most simple propagation methods. For a given propagation
distance `Δz`:

- if `|Δz| < Δzₘ`, then *angular spectrum* propagation method is the most
  appropriate;

- if `|Δz| > Δzₘ`, then single-FFT *Fresnel* propagation method is the most
  appropriate.

Optional argument `D` is the largest lateral width of the beam, the full
lateral size of the field is used by default (`N` is the number of samples
along a dimension of the transverse plane).

"""
function range(F::Field, D::Length = grid_size(F)*grid_step(F))
    D = standardize(F, D)
    λ = wavelength_in_medium(F)
    δx = grid_step(F)
    return D*δx/λ
end

"""
    FourierOptics.apply_quadratic_phase_factor!(F) -> F

applies in-place the quadratic phase factor corresponding to the current
wavefront curvature for the field `F`. If the wavefront curvature is non-zero,
the complex amplitude stored by `F` is modified and the wavefront curvature is
set to zero.

"""
function apply_quadratic_phase_factor!(F::Field)
    curv = curvature(F)
    if !iszero(curv)
        λ = wavelength_in_medium(F)
        δx = grid_step(F)
        ρ = (δx^2/λ)*curv
        apply_quadratic_phase_factor!(F, ρ; planar=true)
    end
    return F
end

"""
    FourierOptics.apply_quadratic_phase_factor!(F, ρ; planar=false) -> F

applies in-place the quadratic phase factor of parameter `ρ` to the complex
amplitude stored by the field `F`.

If keyword `planar` is true, it is assumed that `ρ` takes into account the
wavefront curvature and the wavefront curvature of `F` is set to zero;
otherwise, the wavefront curvature of `F` is left unchanged.

"""
function apply_quadratic_phase_factor!(F::Field{T},
                                       ρ::Dimensionless{Real};
                                       planar::Bool=false) where {T}
    if !iszero(ρ)
        apply_quadratic_phase_factor!(amplitude(F), ρ, RolledCoordinates(F))
    end
    if planar
        # Wavefront curvature has been taken into account by ρ.
        F.curv = zero(F.curv)
    end
    return F
end

"""
    FourierOptics.apply_quadratic_phase_factor!(A, ρ, x) -> A

multiplies in-place the complex amplitude `A[j₁,j₂]` by:

    exp(i⋅π⋅ρ⋅(x[j₁]^2 + x[j₂]^2)) = q[j₁]⋅q[j₂]

for all 2-dimensional indices `(j₁,j₂)` and with:

    q[j] = exp(i⋅π⋅ρ⋅x[j]^2)

Argument `ρ` is a real factor and argument `x` is the vector of coordinates
(assumed to be the same for the 2 axes of `A`).

"""
function apply_quadratic_phase_factor!(A::AbstractMatrix{Complex{T}},
                                       ρ::Dimensionless{Real},
                                       x::AbstractVector{<:Real}) where {T<:AbstractFloat}
    if !iszero(ρ)
        # Check indices.
        J = axes(x, 1)
        axes(A) == (J, J) || throw(DimensionMismatch("complex amplitude has incompatible axes"))

        # Allocate workspace.
        q = similar(A, Complex{T}, (J,))

        # Call unsafe version (which can assume @inbounds) with `ρ` converted
        # to a bare real of suitable floating-point type.
        unsafe_apply_quadratic_phase_factor!(A, as(T, ρ), x, q)
    end
    return A
end

function unsafe_apply_quadratic_phase_factor!(A::AbstractMatrix{Complex{T}},
                                              ρ::T,
                                              x::AbstractVector{<:Real},
                                              q::AbstractVector{Complex{T}}) where {T}
    # NOTE: This method is not called if ρ = 0.
    πρ = π*ρ
    @inbounds @simd for j in eachindex(q, x)
        xⱼ = as(T, x[j])
        q[j] = exp_i(πρ*xⱼ^2)
    end
    J₁, J₂ = axes(A)
    @inbounds for j₂ in J₂
        @simd for j₁ in J₁
            A[j₁,j₂] *= q[j₁]*q[j₂]
        end
    end
    return nothing
end

"""
    FourierOptics.apply_fft!(F::Field, dir) -> F

applies the Fast Fourier Transform (FFT) to the complex amplitude array stored
by the field `F`. Argument `dir` specifies the sign of the argument in the
complex exponential of the transform: if `dir < 0`, the forward FFT is applied;
otherwise, `dir > 0` and the backward FFT is applied.

An error is thrown if the wavefront has a non-zero curvature.

The sampling step is left unchanged.

"""
function apply_fft!(F::Field{T}, dir::Signed) where {T}
    iszero(curvature(F)) || error("FFT not applicable with non-zero wavefront curvature")
    ampl = amplitude(F)
    if dir < zero(dir)
        mul!(ampl, F.forward_plan, ampl)
    elseif dir > zero(dir)
        mul!(ampl, F.backward_plan, ampl)
    else
        throw(ArgumentError("FFT direction has undefinite sign"))
    end
    return F
end
