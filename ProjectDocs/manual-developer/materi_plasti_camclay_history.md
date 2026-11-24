# materi_plasti_camclay_history

## Where

- Initia parsing: `input.cc` (`materi_history_variables` branch). Sets
  `materi_plasti_camclay_history`, `materi_history_variables = 2` and
  creates the history dof (type `MATERI_PLASTI_CAMCLAY_HISTORY`, scalar
  x2) at `hisv_indx`, so the generic history machinery
  (`materi.cc` element_rhside update, `general.cc` inertia/conv_part)
  integrates it.
- Enum: `MATERI_PLASTI_CAMCLAY_HISTORY` in `tochnog.h` +
  `tochnog-mod.h` (kept in sync; full clean build required).
- Names: `database.cc` — name entry + `db_initialize()` basename branch
  `cchis` + counter (cchis0, cchis1).
- Law: `plasti.cc` Cam Clay branch reads `e = old_hisv[0]` and
  `p0 = old_hisv[1]` (stored state). Legacy fallback: if `p0<=0` and
  the group record carries a 4th value `N`, `p0` is derived from the
  normal compression line as before. `stress.cc` unchanged (bulk
  modulus K = (1+e)p/kappa from `old_hisv[0]`).
- `database.cc`: `group_materi_plasti_camclay` is now variable-length
  3..4 (Professional layout `M kappa lambda`).
- `check.cc`: camclay groups accept either `materi_history_variables`
  or `materi_plasti_camclay_history`.

## Design notes

The GNU legacy Cam Clay derived `p0` from the current pressure `p`,
the void ratio `e` and a fourth material parameter `N`
(p0 = exp((N − kappa ln p − (1+e))/(lambda − kappa)), the wall
equation of the NCL). The Professional model (manual 2.2.x) stores
`e` and `p0` as explicit histories and derives `N` from the initial
state; its group record has only `M kappa lambda`. The wall equation
through the current state (e,p) with the SAME `N` is exactly what
relates `e` and `p0` in the Professional — so with `p0` tracked
explicitly both models coincide (verified: final `p0` in the
validation_13 .dbs matches the Professional to the printed precision).

## Gotchas

- Reset indices 1..3 in validation_13 are timeless control steps that
  run before the timestep at index 20; they pre-load the initial state.
- `p0<=0` (uninitialized `cchis1`) with a 3-value group record fails
  with an explicit error, guiding the user to the correct initia.
