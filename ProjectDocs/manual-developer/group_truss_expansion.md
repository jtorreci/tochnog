# group_truss_expansion

## Implementation

- Keyword `group_truss_expansion` (DOUBLE_PRECISION, one value, class
  TRUSS) registered in `database.cc`; enum `GROUP_TRUSS_EXPANSION`
  appended in `tochnog.h`/`tochnog-mod.h` (same enum, same order).
- Read in `truss()` (`truss.cc`) together with the other `group_truss_*`
  records (GET_IF_EXISTS, defaults to 0 = no thermal expansion).
- The thermal correction is subtracted from the incremental length
  BEFORE the force update, for every memory model (the code path is
  shared by `-updated` and the total models):
  `incremental_length -= alpha * dT * initial_length`, where `dT` is the
  average of the nodal temperature increments of the step
  (`new_dof - old_dof` at the `temp_indx` slot, guarded by
  `condif_temperature`).
- The resulting force follows the existing incremental update
  `F_new = F_old + (E*A/L) * incremental_length` and the plastic/rope
  caps, so no other part of the formulation changed.
- Verified against the Professional binary 25-10-2023: `truss11.dat`
  (1D, truss of E=A=alpha=1 with fixed ends, temperature ramp to T=2)
  gives `element_truss_force 3 = -2` EXACT, identical `.dbs` value.
- Blast radius: none of the previously passing truss tests sets
  `group_truss_expansion` (default 0 keeps the code path unchanged);
  truss1/4/5/6/8/14/15 stay rc=0.

## Pending

- The temperature is read per element node and averaged; a temperature
  gradient along the truss uses the mean increment (matches the
  uniform-field corpus tests).
