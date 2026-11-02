# control_print_vtk_coord

## Implementation

- Keyword `control_print_vtk_coord` (data_class CONTROL, type INTEGER,
  data_length 1) registered in `database.cc` with
  `data_required = CONTROL_PRINT_VTK`.
- Enum `CONTROL_PRINT_VTK_COORD` in `tochnog.h` / `tochnog-mod.h`
  (kept in sync), right after `CONTROL_PRINT_VTK`.
- Read in `print_vtk()` (`print_vt.cc`) through `print_option(
  PRINT_VTK_COORD, CONTROL_PRINT_VTK_COORD, icontrol, -YES )`: the control
  record overrules the global `print_vtk_coord`.
- `POINTS` is always written. When the switch is not `-NO`, `POINT_DATA`
  starts with `VECTORS coord double`: the `NODE` coordinates (VERSION_PRINT),
  `0.0` for the components beyond `ndim`. `POINT_DATA` is now also written
  when there is no primary unknown but the coord field is on.

## Decision (C98, owner 2026-11-02 "dale, seguí a Pro")

- Sprint 11 read manual 6.340 ("the coordinates of nodes is plotted in
  vtk") as "omit the `POINTS` block", which gave files a post-processor
  cannot open.
- The Professional 02-08-2026 (F5 W2 V2 probes `vt2`, `cvtk*`) always
  writes `POINTS` and a point field `VECTORS coord` first in `POINT_DATA`;
  `-no` drops only that field. With `-node_deformed_mesh` its `POINTS` are
  deformed and `coord` holds the stored coordinates (`cvtknt`: POINTS x
  0.998714285714, coord x 1). The GNU follows it.
- Effect: every GNU vtk file gains the `coord` field (about `max_node+1`
  lines).

## Details

- The Professional writes the field as `float`, the GNU as `double` like
  its other fields.
- Verification: test-2014 `vtk_coord1` (`vtk_coord1_030.vtk` with the field,
  `vtk_coord1_031.vtk` `-no` without it, `POINTS` kept in both); test-2026
  `print_vtk_coord1/2` (global `-no`, control `-yes` overrule, node 4
  coord `1 1 0.0`).

## Pending

- Nothing for the switch. The gmsh-style coordinate data fields
  (`node_0`/`node_1`/`node_2`) of the gmsh writer remain unimplemented.
