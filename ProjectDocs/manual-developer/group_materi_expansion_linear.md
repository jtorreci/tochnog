# group_materi_expansion_linear (and group_materi_expansion_volume)

## Files and functions

- `stress.cc` — `group_materi_expansion_linear` is read inside the
  `materi_stress` block under `if ( condif_temperature )` (lines 437-448):
  ```c
  get_group_data( GROUP_MATERI_EXPANSION_LINEAR, gr, element, new_unknowns,
    &materi_expansion_linear, ldum, GET_IF_EXISTS );
  for ( idim=0; idim<MDIM; idim++ ) {
    inc_temperature_strain[idim*MDIM+idim] = -materi_expansion_linear *
      ( new_unknowns[temp_indx] - old_unknowns[temp_indx] );
    new_temperature_strain[idim*MDIM+idim] = -materi_expansion_linear *
      new_unknowns[temp_indx];
  }
  ```
  These temperature strains are added to the strain tensor before the stress
  update (`inc_epe`/`new_epe`).
- `materi.cc` — `group_materi_expansion_volume` is read at lines 166-169
  (same `condif_temperature` guard) and applied to the density for gravity
  forces at line 170: `dens = (1.-materi_expansion_volume*temp)*dens;`.
- `database.cc:2346` — `GROUP_MATERI_EXPANSION_LINEAR` (`DOUBLE_PRECISION`,
  length 1). `database.cc:2352` — `GROUP_MATERI_EXPANSION_VOLUME`
  (`DOUBLE_PRECISION`, length 1).
- Enums mirrored in `tochnog.h` and `tochnog-mod.h`.

## Implementation details

- Both records are read with `GET_IF_EXISTS`, so the temperature strain is
  zero and the density is unmodified when the record is missing.
- The MINUS sign on the temperature strain is a sign convention of the strain
  tensor as stored (positive = compressive-like in this codebase's strain
  definition); it must NOT be "fixed" to a plus without re-validating
  `expans1.dat`/`expans2.dat`, which encode the reference behaviour.
- `expans1.dat` covers the linear expansion; `expans2.dat` is a second
  expansion regression; `bimet1.dat` combines linear expansion with
  `condif_temperature` (and is currently failing in `mesh_delete_small` for an
  unrelated test-input reason, not the model).

## External dependencies

- Core `db()`/`get_group_data()` accessors; globals `condif_temperature`,
  `temp_indx`, `MDIM`, `new_unknowns`, `old_unknowns`.
- Density infrastructure `get_materi_density()` and `force_gravity_calculate()`
  in `materi.cc`.

## Hardcoded parameters / pending refactorings

- The expansion is applied only when `condif_temperature` is active; there is
  no "temperature-free" expansion record.
- Both coefficients are element-group indexed; no `no_index` global variant.
- The volumetric expansion only changes the density (gravity forces); it does
  not add a volumetric strain to the constitutive update (that is the job of
  the linear coefficient). Documented here to avoid confusion with the
  professional manual's wording.
