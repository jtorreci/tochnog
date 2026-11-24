# to_pres

## Where

- Pure name entries: `database.cc` TO_PRES/ST_PRES/DY_PRES (`"to_pres"`,
  `"st_pres"`, `"dy_pres"`) next to the `*_sig` family (same pattern:
  no record, only names).
- Item name generation: `calcul.cc` groundflow TOTAL/STATIC/DYNAMIC
  branches generate `to_pres`, `st_pres`, `dy_pres` (Professional
  `post_calcul_label` naming; the GNU used to generate `topres`).
- Slot resolution: `miscel.cc` `exit_tn()` target resolution loop —
  extended from NODE_DOF_CALCUL to the POST_POINT/LINE/QUADRILATERAL
  DOF_CALCUL records: the slot is the position of the item name inside
  `post_calcul_names` (calcul.cc calculate()).

## Design notes

The GNU `db_number()` alias `"topres" -> GROUNDFLOW_PRESSURE` (bounda
dof names, e.g. `bounda_dof ... -topres`) is unrelated and untouched.

## Pending

ground8 (multiple phreatic levels) still fails on the VALUE: the GNU
total pressure at a point below a phreatic level of its element group
is −76.9 vs the Professional −10. Diagnosis: with the multiple
phreatic levels the Professional total pressure at (1,−60)/(1,−20) is
exactly −10 = pres + g·dens·(y_ref−y) with pres = +50/+10 and
y_ref = 0 — the GNU's pres field and/or the static reference differ.
A dedicated A/B against the Professional phreatic-multiple machinery
(groundfl.cc groundflow_phreatic_coord) is needed.
