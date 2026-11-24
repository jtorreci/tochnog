# group_interface_groundflow_total_pressure_tension

## Files and functions

- `interface.cc` — in `interface_element()` (lines 411-455). When
  `groundflow_pressure` is active and
  `GROUP_INTERFACE_GROUNDFLOW_TOTAL_PRESSURE_TENSION` is present, the static
  water pressure from `water_height` is forced on the pressure dofs of the
  interface nodes when the accumulated normal strain exceeds
  `strain_normal_minimum`:
  ```c
  if ( strain_normal > gitpt[0] && dens>0. ) {
    for ( inol=0; inol<nnol; inol++ ) {
      long int jndx = inol*npuknwn + pres_indx/nder;
      double static_pres = force_gravity[ndim-1] * dens * gitpt[1];
      double pres_n = new_dof[inol*nuknwn+pres_indx];
      if ( scalar_dabs(static_pres) > scalar_dabs(pres_n) )
        element_rhside[jndx] += static_pres - pres_n;
    }
  }
  ```
- `database.cc` — keyword registration (alphabetical, between
  `GROUP_INTERFACE_GROUNDFLOW_PERMEABILITY` and
  `GROUP_INTERFACE_MATERI_ELASTI_STIFFNESS`): type `DOUBLE_PRECISION`,
  `data_length = 2`, `data_class = GROUNDFLOW`,
  `data_required = GROUP_INTERFACE`.
- Enum `GROUP_INTERFACE_GROUNDFLOW_TOTAL_PRESSURE_TENSION` in `tochnog.h`.

## Implementation details

- `strain_normal` is the ACCUMULATED interface normal strain (history,
  `ELEMENT_INTERFACE_STRAIN_NORMAL`), read before the current increment is
  applied — same value used by gap/tension/Mohr-Coulomb.
- The correction adds `static_pres - pres_n` to the pressure rhs (a pseudo
  forcing), effectively clamping the crack pressure at the static value.
- Uses the water density (`GROUNDFLOW_DENSITY`) and gravity.

## External dependencies

- Core `db()` accessor; globals `groundflow_pressure`, `pres_indx`, `strain_normal`.

## Hardcoded parameters / pending refactorings

- A strain-based threshold differs from the volume-element counterpart
  (eigenvalue-based); this matches the manual wording for interfaces.
