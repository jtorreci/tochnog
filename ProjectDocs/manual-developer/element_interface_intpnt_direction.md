# element_interface_intpnt_direction + element_interface_intpnt_gap_status (developer)

## Implementation

- `tochnog.h` / `tochnog-mod.h`: enums `ELEMENT_INTERFACE_INTPNT_DIRECTION`,
  `ELEMENT_INTERFACE_INTPNT_GAP_STATUS` before `LAST_DUMMY`.
- `database.cc`: DIRECTION DOUBLE, `data_length` 4*9; GAP_STATUS INTEGER
  (negative enums `-OPENED` / `-CLOSED`), `data_length` 4; both
  `version_all`, `print_only`, class ELEMENT, not fixed length.
- `top.cc`: lengths widened to `ns1_max` / `ns1_max*9` (as the tension
  status) and both records allocated VERSION_NORMAL and VERSION_NEW
  before the parallel element loop.
- `interface.cc`, `interface_element()`, output block: per interface
  integration point (`n_intpnt = nnol/2`) the element `normal`, `tangent`
  and `tangent2` (3D) or `normal x tangent` (2D, `array_outproduct_3D`);
  gap status from `strain_eff_ip <= gap` (the same rule as the stiffness
  and the tension status; an opening by `group_interface_materi_tension`
  limit is not a gap opening). PUT VERSION_NORMAL and VERSION_NEW.
- The directions are those of the element (one normal per interface
  element), equal on every integration point.

## Tests

- `test-2026/element_interface_intpnt_direction1` (2D),
  `element_interface_intpnt_direction_3d1` (hex8),
  `element_interface_intpnt_gap_status1`.

## Pending refactorings

- `element_interface_intpnt_materi_tension_status` uses the gap rule
  too (it does not show an opening by the tension limit).
