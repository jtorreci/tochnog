# group_interface_groundflow_capacity

## Files and functions

- `interface.cc` — in `interface_element()` (lines 368-459). When
  `groundflow_pressure` is active and `GROUP_INTERFACE_GROUNDFLOW_CAPACITY` is
  present, a lumped storage term is added on the pressure dofs of the
  interface nodes:
  ```c
  if ( has_C ) {
    for ( inol=0; inol<nnol; inol++ ) {
      long int jndx = inol*npuknwn + pres_indx/nder;
      double dpres = ( new_dof[inol*nuknwn+pres_indx] -
        old_dof[inol*nuknwn+pres_indx] ) / dtime;
      element_rhside[jndx] -= C_iface * dpres;
    }
  }
  ```
- `database.cc` — keyword registration (alphabetical, between
  `GROUP_INTERFACE_GAP` and `GROUP_INTERFACE_MATERI_ELASTI_STIFFNESS`): type
  `DOUBLE_PRECISION`, `data_length = 1`, `data_class = GROUNDFLOW`,
  `data_required = GROUP_INTERFACE`.
- Enum `GROUP_INTERFACE_GROUNDFLOW_CAPACITY` in `tochnog.h`.

## Implementation details

- Only evaluated when `groundflow_pressure` is initialised; the block guards
  on `has_pe || has_C || has_total_pressure_tension`.
- Lumped (diagonal) storage on the pressure dofs of the interface nodes,
  following the same `dof` indexing as the rest of `interface_element()`.
- `old_dof`/`new_dof` are the time-stepped nodal dofs (per-node arrays).

## External dependencies

- Core `db()` accessor; globals `groundflow_pressure`, `pres_indx`, `dtime`.

## Hardcoded parameters / pending refactorings

- Lumped storage; a consistent (mass-matrix) formulation is not implemented.
