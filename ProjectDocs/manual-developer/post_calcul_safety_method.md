# post_calcul_safety_method — developer

## Files and functions

- `tochnog.h` / `tochnog-mod.h` — enum block of the safety family:
  `POST_CALCUL_SAFETY_METHOD` (-vertical default / -prival / -global),
  `POST_CALCUL_SAFETY_MAXIMUM` (cap, DOUBLE), operator values
  `SAFETY_PIPING`/`SAFETY_LIFTING`, selectors `VERTICAL`/`GLOBAL`, and the
  12 pure-name item entries `SAFETY_PIPING_PRIVAL_0..2`,
  `SAFETY_LIFTING_PRIVAL_0..2`, `SAFETY_PIPING_GLOBAL_X/Y/Z`,
  `SAFETY_LIFTING_GLOBAL_X/Y/Z`.
- `database.cc` — registration block after the
  `POST_CALCUL_STATIC_PRESSURE_HEIGHT` entries: the method/maximum records
  (INTEGER/DOUBLE, no_index 1, class POST), the operator/selector pure
  values and the item pure names (parse + target resolution, same pattern as
  the `*_SIG`/`TO_PRES` item names).
- `calcul.cc` — `calculate()`: naming branch for
  `unknown == -MATERI_STRESS` with operator `-safety_piping`/`-safety_lifting`
  (generates the item names per method); `calculate_operat()`: value branch
  (below).

## Implementation details

- Formulas (measured against the Professional .dbs on ground15/16):
  `safety_lifting = (sigma_i + p_total)/p_total`,
  `safety_piping = (sigma_i + p_dynamic)/p_dynamic` where
  `p_dynamic = p_total - p_static` and `p_total`/`p_static` come from
  `groundflow_phreatic_coord()` (so `p_dynamic == pres_dof` when a level
  applies). Division by a pressure with `|p| < TINY` yields 0
  (post_calcul_safety_default: eps very small, value 0); the optional
  `post_calcul_safety_maximum` caps the result.
- `sigma_i` selection: `-vertical` reads the diagonal entry of the
  `(ndim-1)` axis of the sigma matrix (sigma_yy 2D / sigma_zz 3D /
  sigma_xx 1D); `-global` reads the diagonal entries sigma_xx/sigma_yy/
  sigma_zz (slots `_x`, `_y`, `_z`); `-prival` uses the principal values
  sorted from the MOST compressive. GOTCHA: the GNU `sort()` (math.cc)
  orders the principal values DESCENDING (most positive first), while the
  Professional lists them from the most compressive (`prival_0` = the
  minimum, verified on ground15 where only the sigma_yy = -10 slot yields
  2.3333); the safety branch therefore reads `prival[MDIM-1-j]`.
- Result length: 1 value (vertical) or MDIM (prival/global), matching the
  number of generated item names per operator.

## Verification

- ground13 (no safety) passes with `-topres`; the safety pipeline runs on
  ground14/15/16 (vertical/prival/global). The factors converge to the
  Professional targets (piping 2.3333, lifting 2.0) once the coupled
  consolidation state reaches the drained steady pore-pressure field; within
  the 1 s window of the corpus tests the GNU u-p state is still transient
  (see SEGUIMIENTO: the drained targets are reproduced exactly at ~2000 s).
