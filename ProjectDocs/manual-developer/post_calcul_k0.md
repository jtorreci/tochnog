# post_calcul -k0 operator

## Implementation

- Enums appended in `tochnog.h`/`tochnog-mod.h` (same enum, same
  order): `K0` (the `-k0` operator token, registered in `database.cc`
  as an INTEGER pure name entry) and `K0_SIG` (the output item name
  `k0_sig`, resolved by `target_item`).
- Slot generation in `calculate()` (`calcul.cc`): a branch for
  `unknown==-MATERI_STRESS && labs(calcul_operat)==K0` creates ONE
  scalar post_calcul slot named `k0_sig` (matching the Professional
  `post_calcul_label`, measured on `k0.dat`).
- Value computation in `calculate_operat()` (`calcul.cc`): new branch
  for `labs(calcul_operat)==K0` on a matrix unknown:
  - 2D: `0.5*(uv[0]+uv[8])/uv[4]`  (sigxx, sigzz over sigyy)
  - 3D: `0.5*(uv[0]+uv[4])/uv[8]`  (sigxx, sigyy over sigzz)
  - 1D fallback: `0.5*(uv[4]+uv[8])/uv[0]`
  - zero vertical stress gives 0 (guarded division).
- Verified against the Professional binary 25-10-2023: `k0.dat` gives
  `node_dof_calcul = post_point_dof_calcul = 0.1111111111` EXACT
  (nu/(1-nu) of the laterally confined column).
- Blast radius: only activated when a `post_calcul` line carries `-k0`
  (no previous corpus test used it); the apparent
  young/poisson branches and all other operators untouched.

## Pending

- The 1D branch formula is a documented fallback (the Professional
  manual only defines 2D/3D).
