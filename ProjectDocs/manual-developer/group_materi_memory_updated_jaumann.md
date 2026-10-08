# group_materi_memory -updated_jaumann

## Implementation (F3, 2026-10-07, BACKLOG C20)

- **Enum**: `UPDATED_JAUMANN` appended right before `LAST_DUMMY` in BOTH
  `tochnog.h` and `tochnog-mod.h` (same position, so existing values stay
  stable). Keyword `updated_jaumann` registered in `database.cc`
  (`db_initialize`, next to `updated_linear`).
- **materi.cc, `set_deften_etc()`**:
  - old/new rotation matrices: identity (same branch as
    `-UPDATED_LINEAR`/`-TOTAL_LINEAR`);
  - incremental rotation: `inc_rot = I + 0.5 (inc_deften - inc_deften^T)`
    (manual 2.2.1, "R = 0.5(F - F^T) + I"); every other memory keeps
    `set_deften_u_rot(inc_deften, inc_u, inc_rot)` (polar rotation);
  - strain: linear engineering increment (same list as `-UPDATED_LINEAR`).
- **materi.cc, `materi()`**: the old stress (and rho, epi, maxwell
  stresses, damage effective stress) is taken unrotated, as for the other
  updated memories; after `set_stress()` the new stress (old + extra) is
  rotated with `inc_rot` (same list as `-UPDATED`), i.e. the Jaumann
  first-order update.
- **stress.cc, `set_stress()`**: compressibility block treats it as the
  linear memories (increment from `inc_epe`).
- **tendon.cc**: the assert of the no-rotation tendon branch accepts it.
- Not touched: the explicit-memory `materi_displacement` restriction of
  stress.cc only covers `-UPDATED`/`-UPDATED_WITHOUT_ROTATION`, so
  `-updated_jaumann` with `materi_displacement` is accepted, as in the
  Professional; elem.cc updates the coordinates for it (every memory but
  `-TOTAL_LINEAR`).

## Verification

Professional 02-08-2026 vs GNU (simple shear gamma 0.5, single quad4,
node 4):

| case | sigxx | sigxy | sigyy |
|---|---|---|---|
| 1 step, Pro | 0.125 | 0.234375 | -0.125 |
| 1 step, GNU | 0.125 | 0.234375 | -0.125 |
| 100 steps, Pro | 0.06183352822 | 0.2396333009 | -0.06183352822 |
| 100 steps, GNU | 0.06183352822 | 0.2396333009 | -0.06183352822 |

Identical for `materi_velocity_integrated` with `-fixed_in_space` and
`-follow_material`, and for `materi_displacement`. Two-quad4 column
sheared by 1 in one step: `-fixed_in_space` GNU sigxx -0.02164855 / Pro
-0.02164855; `materi_displacement` GNU 0.0207070 / Pro 0.0207068;
`-follow_material` GNU 0.0207 / Pro -0.0216 (in-step geometry
divergence, see DIVERGENCES-PROFESSIONAL.md).

Tests: `validation-suite/test-2026/materi_memory_updated_jaumann_shear1`
and `_shear2` (rc=1 before, rc=0 after), registered in
`scripts/build_safe.sh`.

## Hardcoded parameters / pending

- The first-order rotation is not orthogonal for finite increments (by
  definition of the Professional formula); no re-orthogonalisation.

## Differences with the Professional version

- **explicit `-follow_material`: in-step geometry** (known difference).
  GNU end-of-step geometry vs Professional start-of-step geometry; see
  the user page. Status: Open, LOW. Index:
  [DIVERGENCES-PROFESSIONAL.md](../DIVERGENCES-PROFESSIONAL.md).
