# group_materi_plasti_diprisco_density (developer)

## Where

- `database.cc`: `GROUP_MATERI_PLASTI_DIPRISCO_DENSITY` (DOUBLE, 22,
  MATERI); enums in `tochnog.h` / `tochnog-mod.h` after
  `GROUP_MATERI_PLASTI_DIPRISCO_RT`; `check.cc` requires `materi_stress`,
  `materi_history_variables` or `materi_plasti_diprisco_history`, and
  `materi_strain_plasti`.
- `plasti.cc` `plasti_rule()`: before the di Prisco block. When the record
  exists, the ten data are replaced by `(1-w) loose + w dense`,
  `w = (el - e_old)/(el - ed)` clamped to [0,1], `e_old = old_hisv[11]`,
  and `new_hisv[11] = e_old + (1 + e_old) tr(inc_ept)` (total strain
  increment, tochnog sign). The di Prisco block then runs unchanged
  (`prisco_record` selects the record for `db_error`); the model keeps the
  plasti type `GROUP_MATERI_PLASTI_DIPRISCO`.
- `materi.cc`: the per-model dof `materi_strain_plasti_diprisco` is also
  filled for density groups.
- `conspr.cc`: contact spring friction approximation, same 1 degree as the
  plain di Prisco model.

## Implementation details

- Errors: `el <= ed` (Professional message) and fewer than 12 history
  variables.
- `plasti_rule()` runs at least once per step (yield test), so the void
  ratio is updated in elastic steps too; with `plasti_apply -no` it is not
  updated.
- The void ratio update is incremental (explicit); the weight uses the
  value at the start of the step.

## Verification

diprisc1 rig (isochoric triaxial, t 0.004): e 0.75 gives exactly the plain
model with the averaged parameters; e 0.9 the loose set; e 0.5 (clamped) the
dense set; oedometric compression gives e = 1.75 (1 - 1.e-5)^400 - 1
exactly. The Professional uses the loose set only (DIVERGENCES
DIPRISCO-DENSITY).

## Differences with the Professional version

- The Professional 02-08-2026 uses the loose set only (DIVERGENCES
  DIPRISCO-DENSITY, class D, owner decision 2026-11-07).

## Pending

- The plain di Prisco model itself differs from the Professional by about
  0.2 % on this rig (BACKLOG C119).
