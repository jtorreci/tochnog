# control_print_vtk_dof_calcul

## Implementación

- Keyword `control_print_vtk_dof_calcul` (data_class CONTROL, type
  INTEGER, data_length DATA_ITEM_SIZE, fixed_length 0) registered in
  `database.cc` with `data_required = CONTROL_PRINT_VTK`.
- New enum `CONTROL_PRINT_VTK_DOF_CALCUL` in `tochnog.h` /
  `tochnog-mod.h` (kept in sync).
- Read in `print_vtk()` (`print_vt.cc`) like `control_print_vtk_dof`
  (array + `nvtk_dof_calcul` count, `GET_IF_EXISTS`).
- A per-field flag array `print_post_field[icalcul]` is computed once
  after reading `POST_CALCUL_SCAL_VEC_MAT` and gates BOTH post loops
  (scalars/vectors block and component block), while `icalcul`/`idim`
  always advance (same pattern as the dof filter).

## Matching rule (GNU adaptation, documented)

- The Professional manual says "See `post_calcul_label` for the allowed
  names". `post_calcul_label` does NOT exist in the GNU; the post field
  labels are the `post_calcul_names` strings built in `calcul.cc`.
- The names the user can type are resolved by `db_number` (standard
  parser). A post field `icalcul` is written when ANY listed name
  matches:
  1. **exact** on the underlying unknown
     `post_calcul_unknown_operat[icalcul*2+0]` (initia names like
     `-materi_stress` select every operator of that data), or
  2. **substring** of the field label `post_calcul_names[icalcul]`
     (e.g. `-sigyy` matches `sigyy` and `tosigyy`).
- `-none` (`vtk_dof_calcul[0]==-NONE`) sets every flag to 0.
- Without the record every flag is 1 (default behavior unchanged).

## Detalles / gotchas

- The vtk post block writes TWO blocks per field: the aggregated
  scalars/vectors (header `db_name(unknown)_db_name(operat)`, e.g.
  `materi_strain_total_average`) and the component scalars (header
  `post_calcul_names[icalcul]`, e.g. `aept`). Both are gated by the
  same flag.
- The vector `idim` extension counter advances even for filtered-out
  vector fields (keeps the `_0/_1/_2` numbering consistent).
- `post_calcul` requires an operator in the GNU (there is no bare
  `-materi_stress` default); the filter therefore works per
  unknown+operator pair.
- Verification: `vtk_dofcalc1` test — no filter (both post headers),
  `-none` (no post headers, primary fields intact), filter
  `-materi_strain_total` (only `materi_strain_total_average`).
