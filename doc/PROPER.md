# Random notes about `PROPOER` library

## FFT

In IDL 2-D forward FFT is given by `fft_x = FFT(x, -1)` and is defined by:

```
fft_x[k1,k2] = 1/(n1*n2) * sum_{j1,j2} x[j1,j2]*exp(-2iπ*(j1*k1/n1 + j2*k2/n2))
```

the backward (inverse) FFT is given by `x = FFT(fft_x, +1)` and is defined by:

```
x[j1,j2] = sum_{k1,k2} fft_x[k1,k2]*exp(+2iπ*(j1*k1/n1 + j2*k2/n2))
```

In `prop_stw.pro` (with some simplifications to make clear what is done, `dz`
is the propagation distance):

``` idl
direct = dz ge 0
if direct then begin
    a.wavefront =  fft(a.wavefront, -1)*n
endif else begin
    a.wavefront =  fft(a.wavefront, +1)/n
endelse
```

Since the grid is square, `n1 = n2 = n`, what is done is:

```
w_out[k1,k2] = 1/sqrt(n1*n2) * sum_{j1,j2} w_in[j1,j2]*exp(-2*sign(dz)*i*π*(j1*k1/n1 + j2*k2/n2))
```

i.e. the orthonormal transform is applied.

In FFTW: `fft` (the forward transform) multiplies by `exp(-2iπ...)` while
`bfft` (the backward transform) multiplies by `exp(-2iπ...)`, both with no
scaling.


## Propagation

Notations:

- `w0` = minimum possible waist radius for current beam
- `z_w0` = location of minimum beam waist
- `z_Rayleigh = π*w0^2/λ` = Rayleigh distance from minimum beam waist
- `w_at_surface` = waist radius at current surface
- `R_beam` = beam radius of curvature

Structure members used for propagation:

- `propagator_type` (string):
  - is the concatenation of the input beam type (`INSIDE` or `OUTSIDE`), of
    `_TO_`, and of the output beam type (`INSIDE` or `OUTSIDE`);
  - is initialized to `INSIDE_TO_INSIDE` in `prop_begin`;
  - is changed by `prop_select_propagator` and `prop_lens`;
  - is used by `prop_propagate` and `prop_lens` to modify the wavefront.

- `reference_surface` (string):
  - is the shape of the reference surface;
  - is `PLANAR` or `SPHERI`;
  - is initialized to `PLANAR` in `prop_begin`;
  - is changed by `prop_ptp` (`PLANAR -> PLANAR`), `prop_stw` (`SPHERI ->
    PLANAR`), `prop_wts` (`PLANAR -> SPHERI`), and `prop_lens` (to `PLANAR` if
    output beam type `beam_type_new` is `INSIDE`, to `SPHERI` otherwise);

- `beam_type_old` (string):
  - is the output beam type after propagation or insertion of a lens;
  - is `INSIDE` or `OUTSIDE`.
  - is initialized to `PLANAR` in `prop_begin`;
  - is changed by `prop_ptp` (`PLANAR -> PLANAR`), `prop_stw` (`SPHERI ->
    PLANAR`), `prop_wts` (`PLANAR -> SPHERI`), and `prop_lens` (to `PLANAR` if
    output beam type `beam_type_new` is `INSIADE`, to `SPHERI` otherwise);

- `R_beam` (float)
 - is the beam radius of curvature;
 - is initialized to `0` in `prop_begin`;
 - is updated in `prop_select_propagator` and `prop_lens`.

- `R_beam_inf` (boolean)
 - indicates whether the beam radius of curvature is infinite;
 - is initialized to `true` in `prop_begin`;
 - is updated in `prop_select_propagator` and `prop_lens` (`true` if output
   beam is spherical, `false` if it is planar);

In the following, `PROPER` methods for propagation are given in Julia
pseudo-code (omitting messages and assuming all parameters are stored in
wavefront instance `A`):

- plane-to-plane (PTP) propagator (2 FFTs):

  ```julia
  @assert A.reference_surface == PLANAR
  A.reference_surface = PLANAR
  A.z += dz
  A.wavefront = idl_fft(A.wavefront, -1) * n
  x = (indgen(n) - n/2)/(n*A.dx)
  y = x'
  rhosqr = x^2 + y^2
  A.wavefront *= exp_i(-π*λ*dz*rhosqr)
  A.wavefront = idl_fft(A.wavefront, +1) / n
  maybe_apply_phase_shift!(a, exp_i(2π*dz/λ))
  ```

- spherical-to-waist (STW) propagator:

  ```julia
  dz = A.z_w0 - A.z # by default
  @assert A.reference_surface == SHERICAL
  A.z += dz
  A.dx = λ * abs(dz) / (n * A.dx)
  direct = dz ≥ 0
  if direct # forward transform (-1 for IDL FFT)
      A.wavefront = idl_fft(A.wavefront, -1) * n
  else # backward transform (+1 for IDL FFT)
      A.wavefront = idl_fft(A.wavefront, +1) / n
  end
  _qphase!(a, dz)
  maybe_apply_phase_shift!(a, exp_i(2π*dz/λ))
  A.reference_surface = PLANAR
  ```

- waist-to-spherical (WTS) propagator:

  ```julia
  @assert A.reference_surface == PLANAR
  A.reference_surface = SPHERICAL
  iszero(dz) && return
  λ = get_wavelength(A)
  direct = dz ≥ 0
  A.z += dz
  _qphase!(a, dz)
  if direct # forward transform (-1 for IDL FFT)
      A.wavefront = idl_fft(A.wavefront, -1) * n
  else # backward transform (+1 for IDL FFT)
      A.wavefront = idl_fft(A.wavefront, +1) / n
  end
  maybe_apply_phase_shift!(a, exp_i(2π*dz/λ))
  A.dx = λ * abs(dz) / (n * A.dx)
  ```

- `prop_select_propagator`:

  ```julia
  dzw = A.z_w0 - A.z
  newz = A.z + dz
  if abs(A.z_w0 - newz) < A.rayleigh_factor*A.z_Rayleigh
      beam_type_new = INSIDE
  else
      beam_type_new = OUTSIDE
  end
  A.propagator_type = A.beam_type_old : beam_type_new
  A.beam_type_old = beam_type_new
  return dzw
  ```

  Notes:

  1. This subroutine is only called by `prop_propagate`.
  2. Computed value `dzw`, the distance to new focus position from new position,
     is never used.
  3. The only purpose of this subroutine is to determine whether propagation is
     to *near* or *far* field (respectively called `INSIDE` and `OUTSIDE` in
     the code).
  4. Side effects: internal properties `propagator_type` and `beam_type_old`
     are modified but could be managed in a more simple way.

  The following Julia function does the job with no side effects:

  ```julia
  @enum FieldRange::UInt32; NEAR = 0; FAR = 1; end
  field_range(A::Wavefront, dz::Length) =
      abs((A.z + dz) - A.z_w0) < A.rayleigh_factor*A.z_Rayleigh ? NEAR : FAR
  ```

- `prop_propagate`:
  ```julia
  dzw = _select_propagator!(A, dz)
  z1 = A.z
  z2 = z1 + dz
  if to_plane
      # Force output beam type to be PLANAR.
      A.propagator_type = first(A.propagator_type) : INSIDE
  end
  if A.propagator_type === INSIDE_TO_INSIDE
      apply_planar_to_planar_propagator!(A, dz)            # PLANAR    -> PLANAR
  elseif A.propagator_type === INSIDE_TO_OUTSIDE
      apply_planar_to_planar_propagator!(A, A.z_w0 - z1)   # PLANAR    -> PLANAR
      apply_waist_to_spherical_propagator!(A, z2 - A.z_w0) # PLANAR    -> SPHERICAL
   elseif A.propagator_type === OUTSIDE_TO_INSIDE
      apply_spherical_to_waist_propagator!(A, A.z_w0 - z1) # SPHERICAL -> PLANAR
      apply_planar_to_planar_propagator!(A, z2 - A.z_w0)   # PLANAR    -> PLANAR
  elseif A.propagator_type === OUTSIDE_TO_OUTSIDE
      apply_spherical_to_waist_propagator!(A, A.z_w0 - z1) # SPHERICAL -> PLANAR
      apply_waist_to_spherical_propagator!(A, z2 - A.z_w0) # PLANAR    -> SPHERICAL
  else
      error("invalid propagator type: `$(A.propagator_type)`")
  end
  ```

  Notes:
  1. `...TO_INSIDE` (that is propagation to near field) always results in a
     `PLANAR` reference surface. Conversely `...TO_OUTSIDE` (that is
     propagation to far field) always result in a `SPHERICAL` reference
     surface. There is therefore some redundancy in the properties
     `beam_type_old`, `reference_surface`, and `propagator_type`.
  2. Simplifications:
     ```julia
     dz1 = A.z_w0 - z1 = A.z_w0 - A.z
     dz2 = z2 - A.z_w0 = A.z_w0 - A.z + dz - A.z_w0 = dz - A.z
     ```

  Can be simplified to:
  ```julia
  @enum Surface::UInt32 begin; PLANAR = 0; SPHERICAL = 1; end
  from = A.current_field_range
  to = to_planar ? NEAR : field_range(A, dz)
  dz1 = A.z_w0 - A.z # 1st propagation distance
  dz2 = dz - A.z     # 2nd propagation distance
  if A.reference_surface === PLANAR
      # Current reference surface is planar, previous propagation was to near field.
      if to === NEAR
          # Apply near-to-near propagator.
          apply_planar_to_planar_propagator!(A, dz)
      else # to === FAR
          # Apply near-to-far propagator.
          apply_planar_to_planar_propagator!(A, dz1)
          apply_waist_to_spherical_propagator!(A, dz2)
      end
  else
      # Current reference surface is spherical, previous propagation was to far field.
      if to === NEAR
          # Apply far-to-near propagator.
          apply_spherical_to_waist_propagator!(A, dz1)
          apply_planar_to_planar_propagator!(A, dz2)
      else # to === FAR
          # Apply far-to-far propagator.
          apply_spherical_to_waist_propagator!(A, dz1)
          apply_waist_to_spherical_propagator!(A, dz2)
      end
  end
  @assert A.reference_surface === (to === NEAR ? PLANAR : SPHERICAL)
  ```
