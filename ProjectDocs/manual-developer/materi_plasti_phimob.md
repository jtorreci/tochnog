# materi_plasti_phimob

## Files and functions

- `input.cc`: initialization word `materi_plasti_phimob` → flag
  `materi_plasti_phimob`, index `phimob_indx`, one SCALAR dof of type
  `-MATERI_PLASTI_PHIMOB`; requires `materi_stress`.
- `tochnog.h` / `tochnog-mod.h`: enum (F5 T3 C3 block), flag and index
  externs; `initia.cc` definitions.
- `database.cc`: name; dof label `phimob`. The label reuses the number of
  the GNU `post_calcul` operator `PHIMOB` (same word, which `db_number()`
  finds first), so `-phimob` names the dof in `target_item` / print
  records and still is the operator in `post_calcul`.
- `general.cc`: inertia 1 (state dof, like `materi_strainenergy`);
  `check.cc`: check_unknown.
- `materi.cc`: `materi_plasti_phimob_value()` and the fill next to
  `materi_strainenergy` (`volume * h_lump * (new - old) / dtime`).

## Implementation details

- c from `group_materi_plasti_mohr_coul`, else
  `group_materi_plasti_mohr_coul_direct` (the record value, not the wall
  reduced one); no record → 0.
- Root: `b sin(phi) - c cos(phi) = r sin(phi - delta)`, `r = sqrt(b^2+c^2)`,
  `delta = atan2(c, b)`, roots `delta + asin(-a/r)` and
  `delta + pi - asin(-a/r)`; the root in [-pi/2, pi/2] with the smallest
  absolute value; 90 degrees when `|a| > r`; 0 when `r = 0`.

## Verification

Professional 02-08-2026 probes (t3/c3) on the `plasti_apply1` square:
yield 28.6479 (GNU 28.6479), elastic c 1 5.39561 (GNU 5.39561),
`mohr_coul_direct` 5.39561, hardening/softening and Drucker-Prager 0. With
c 0 in tension the stress states of the two codes differ (pre-existing
Mohr-Coulomb difference); on the Professional stress the GNU formula gives
its -18.6186.

## Pending

- `group_materi_plasti_mohr_coul_hardening_softening` gives 0, as the
  Professional, although it is a Mohr-Coulomb law.
