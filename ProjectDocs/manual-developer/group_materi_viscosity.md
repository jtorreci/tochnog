# group_materi_viscosity (and heatgeneration / user)

## Files and functions

- `viscosit.cc` — `viscous_stress()` (line 22): computes the viscous stress.
  - Reads `group_materi_viscosity` (`nu`) and `group_materi_viscosity_user`;
    if user viscosity is `-yes`, calls `user_viscosity()` to override `nu`.
  - Builds `D[idim,jdim] = 0.5*(du_i/dx_j + du_j/dx_i)` from
    `grad_unknowns[...vel_indx...]`.
  - `array_multiply(D, new_sig, 2.*viscosity, MDIM*MDIM)` — the viscous stress
    is ADDED to `new_sig` in-place.
  - If `group_materi_viscosity_heatgeneration -yes`, sets
    `viscosity_heatgeneration = 2*viscosity * (D:D)`.
- Entry point called from `stress.cc:686-689` (old and new configuration)
  inside the `materi_stress` block.
- `user.cc:167` — `user_viscosity(user_data, new_unknowns, &visc)`. The shipped
  body is a stub: it prints "Error: routine user_viscosity not programmed" and
  exits. The header comment documents the expected interface (e.g.
  `visc = user_data[0]*(user_data[1]-temp)` using `new_unknowns[6]` when
  `condif_temperature` is initialized after `materi_strain_total`).
- `database.cc:2749-2761` — registrations:
  - `GROUP_MATERI_VISCOSITY`: `DOUBLE_PRECISION`, length 1.
  - `GROUP_MATERI_VISCOSITY_HEATGENERATION`: `INTEGER`, length 1.
  - `GROUP_MATERI_VISCOSITY_USER`: `INTEGER`, length 1.
- Enums mirrored in `tochnog.h` and `tochnog-mod.h`.

## Implementation details

- `viscous_stress` only fills the deviatoric viscous part from the symmetric
  velocity gradient; the volumetric (compressibility) behaviour comes from the
  `group_materi_elasti_compressibility` record used by the surrounding stress
  update.
- `user_data` passed to `user_viscosity` comes from the calling context in
  `stress.cc` (the element's `user_data`, configured via the `user_data` input
  mechanism).
- `ndim` limits the loop over dimensions; the `D` tensor is stored in `MDIM*MDIM`
  arrays, so in 2D the out-of-plane entries stay 0.

## External dependencies

- Core `db()`/`get_group_data()` accessors; globals `ndim`, `vel_indx`,
  `nuknwn`, `nder`, `MDIM`, `new_sig`.

## Hardcoded parameters / pending refactorings

- `user_viscosity()` is a hardcoded stub — the only way to use
  `group_materi_viscosity_user` is to edit `user.cc` and recompile (no runtime
  plugin mechanism).
- The heat-generation flag is read with plain `db(...GET_IF_EXISTS)` (not
  `get_group_data`); both work for an INTEGER group record, but the two read
  styles are inconsistent within the same routine.
- Viscosity data is element-group indexed only; there is no `no_index` variant
  to apply a single viscosity to the whole model.
