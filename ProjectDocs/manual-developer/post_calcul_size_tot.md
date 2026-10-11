# post_calcul_size_tot

## Where

`calcul.cc`, `calculate_operat()`: branches `SIZETOT` and `SIZEDEV`
(names `-sizetot`/`-size_tot`, `-sizedev`/`-size_dev`; the Professional
spellings are translated in `db_number()`).

## Implementation

- `SIZETOT`: `sqrt( array_inproduct(A, A, MDIM*MDIM) )` for a matrix,
  `sqrt( array_inproduct(a, a, ndim) )` for a vector.
- `SIZEDEV`: `sqrt( array_inproduct(B, B, MDIM*MDIM) )`, B = A minus the
  average on the diagonal.
- F4 T6a (2026-10-11): the factor 0.5 inside both square roots was
  removed (manual Professional `post_calcul`; Professional 02-08-2026
  `calcul1`: 1 and 0.816497; DIVERGENCES K12; the validation_3 /
  validation_14 calcul ratio 1/sqrt(2)).
- `refine_l.cc` `element_residue_norm_set()` also uses `-SIZETOT` for
  the refinement indicator; the selection compares each element with
  the others and with `EPS_NORM` times their average, so the uniform
  scale change does not change which elements are refined.

## Pending refactorings

- `-prival` ordering (largest first in the GNU, smallest first in the
  Professional) is a separate item (K12).
