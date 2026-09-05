# group_truss_initial_force

## Implementation

- Keyword `group_truss_initial_force` (DOUBLE_PRECISION, one value,
  class TRUSS) registered in `database.cc`; enum
  `GROUP_TRUSS_INITIAL_FORCE` appended in `tochnog.h`/`tochnog-mod.h`
  (same enum, same order).
- Read in `truss()` (`truss.cc`) with the other `group_truss_*` records
  (GET_IF_EXISTS, default 0).
- Seeding: when the stored state `ELEMENT_TRUSS_FORCE` does not exist
  yet (the element has just come to life, i.e. its first evaluation —
  also covers elements born later in phased analyses), the "old" force
  of the incremental update is initialized with
  `group_truss_initial_force` instead of 0:
  `F_new = initial_force + (E*A/L) * incremental_length` for the birth
  step. The state written at the end of the step carries the force, so
  the following steps update from it normally.
- No initial strain: the record only seeds the force state; the
  mechanical strain develops from the actual deformation (matches the
  Professional: `element_truss_strain = 0` with a stationary
  initial-force truss).
- Verified against the Professional binary 25-10-2023: `truss12.dat`
  gives the nodal right-hand-sides +1/-1 in vely and force +1 EXACT.
  Note: the Professional renumbers the truss element 2 -> 3 in its
  `.dbs` output (its distribute machinery); the GNU keeps the input
  numbering — the physics is identical.
- Blast radius: only elements of truss groups that specify the record
  are affected (none of the previously passing truss tests use it);
  truss1/4/5/6/8/14/15 stay rc=0.

## Pending

- The initial force enters the state at the element's first evaluation;
  a truss group whose elements never participate keeps no force record
  (harmless).
