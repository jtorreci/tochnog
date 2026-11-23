# strain_volume_* — developer

## Files and functions

- `database.cc` — `STRAIN_VOLUME_RELATIVE_TIME` / `_ABSOLUTE_TIME`
  (DOUBLE, variable length), `_ELEMENT` / `_ELEMENT_GROUP` (INTEGER,
  variable length), `_GEOMETRY` (INTEGER, 2); internal per-element
  `ELEMENT_STRAIN_VOLUME`. Enums in `tochnog.h` / `tochnog-mod.h`.
- `materi.cc`:
  - `strain_volume_setup()` (called at the end of `input()`): rejects
    both tables in one index and stores `POST_STRAIN_VOLUME_INITIAL` of
    every index, integrating the element volumes with `pol()` on the
    input coordinates (the absolute tables need it from the first step).
  - `strain_volume_in_set(ind, element, group)`: intersection of the
    given selectors, the element list, the group list and the geometry
    (`element_in_geometry -all`); no selector = every element.
  - `strain_volume_eigen(element, group, old, new)`: eps_v at
    `time_current` and `time_current + dtime`, last index wins.
  - in `materi()`, after `set_deften_etc`: copies of `inc_ept` / `new_ept`
    minus (new-old)/3 and new/3 on the diagonal go to `set_stress`; the
    actual strains stay for the strain dofs. Each integration point adds
    `volume * (tr F - 3)` to `ELEMENT_STRAIN_VOLUME`.
- `top.cc` — preallocates `ELEMENT_STRAIN_VOLUME` before the parallel
  element loop when a `strain_volume_` record exists.
- `post.cc` — `post_count_data_volume()`: `POST_STRAIN_VOLUME_ABSOLUTE`
  = sum of `ELEMENT_STRAIN_VOLUME` over the set, `_RELATIVE` =
  100 * absolute / initial.

## Implementation details

Professional 02-08-2026 probes (hex8 / quad4, E 1e6, nu 0.2,
`-total_linear`): fully fixed hex8 3 % -> -16666.67 in xx, yy, zz; free
hex8 absolute 0.1234 -> post 0.1234 (the small strain trace, not the exact
deformed volume 0.1247); plane strain two quad4 -> 0.0454000102547 at t = 2;
geometry needs all nodes inside (a centroid inside is not enough). GNU
reproduces each value to 1e-8.

## History

- F5 T2 V9 (2026-11-04): `post_strain_volume_initial` only.
- F5 T5 U8 (2026-11-22): the volume strain, the remaining records and the
  post results.
- 2026-11-23: selectors combined as an intersection (Professional
  02-08-2026: element 1 with groups 0 1 -> element 1, initial 4; element 2
  with a geometry holding element 1 -> empty set, all results 0). Tests
  `strain_volume_selector_and1/2`.

## Hardcoded parameters / pending refactorings

- Large deformation memories are not rejected (the manual restricts the
  option to small deformation).
- The element set is not recomputed after mesh refinement.
