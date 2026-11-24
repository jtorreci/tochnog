# group_groundflow_permeability

## Where implemented

- `database.cc` (db_initialize): record metadata. `data_length = ndim`
  (one slot per space direction), `fixed_length = 0` — variable length
  1..ndim: the Professional (manual 6.618) allows a single value that is
  then used in EVERY direction.
- `groundda.cc` (`groundflow_data`): the read now uses a dedicated
  length variable (`permeability_length`); when exactly ONE value was
  stored it is replicated over the ndim directions (`pe[idim] = pe[0]`).
- `groundfl.cc`: presence gate (`db_active_index`) unchanged; the
  velocity/flux assembly consumes the replicated `pe[]`.

## Implementation details

Record storage: the database layer marks variable-length records with a
`DBL_MAX` sentinel at the first unused slot; `db_len()`/`db()` GET
return the ACTUAL stored length (1..ndim). Fixed-length inputs (exactly
ndim values, e.g. `kx ky` in 2D) behave exactly as before — the
variable-length parse stops at the next keyword, which is how every
data record terminates.

Consumers of the raw length semantics (dependency diagrams over
permeability, e.g. `ground11_nonsaturated`) are unchanged: the
dependency branch in `get_group_data` derives the value count from the
diagram size, independent of `fixed_length`.

## Verification

- `ground19_water_under_dam` (2D seepage under a dam with a concrete
  cut-off wall, single isotropic value 1.01937e-4): post_node_result
  flux out of the right top edge = -0.00557863 vs target
  -0.00557268 ± 1.e-5 (the Professional itself gives -0.00557268).
- `large2` / `large3` (3D consolidation bricks, single value in 3D):
  run to completion (rc=0; no targets, solver/memory tests).
- Existing 1D single-value and 2D two-value tests (ground1/2/3/6,
  matrix8, matrix4...) unchanged (all rc=0).

## Pending

- `group_groundflow_permeability_nonlinear_method` /
  `_nonlinear_parameters` (nonlinear permeability iteration) and
  `group_groundflow_permeability_vertical_stress` remain as before
  (the latter implemented, non-linear method records registered as
  partials).
