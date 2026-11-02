# element_interface_intpnt_direction + element_interface_intpnt_gap_status

## Description

Output records of an interface element (manual Professional 6.430 and
6.431), filled after every time step:

- `element_interface_intpnt_direction <element> n_x n_y n_z t1_x t1_y t1_z
  t2_x t2_y t2_z ...`: per integration point of the interface the normal
  direction, the first tangential direction and the second tangential
  direction (normal x first tangential). In 2D the z components of the
  normal and the first tangential are 0 and the second tangential is the
  z axis with its sign.
- `element_interface_intpnt_gap_status <element> status ...`: per
  integration point `-closed` when the accumulated normal strain is at
  most `group_interface_gap` (default 0), else `-opened`.

The interface has one integration point per pair of facing nodes
(quad4: 2, quad6: 3, prism6: 3, hex8: 4).

## Example

`element_interface_intpnt_direction1.dat` (validation-suite/test-2026),
a 2D quad4 interface along (0.6, 0.8):

```
element_interface_intpnt_direction 2 -0.8 0.6 0 0.6 0.8 0 0 0 -1 -0.8 0.6 0 0.6 0.8 0 0 0 -1
```

`element_interface_intpnt_gap_status1.dat`: an interface with
`group_interface_gap 0 -1.e-4` opened by 1e-3 gives `-opened -opened`; a
compressed one gives `-closed -closed`.

## Differences with the Professional version

- The Professional 02-08-2026 writes neither record (also with an opened
  gap). It writes `element_interface_direction` instead (not a GNU
  record), with the same vectors in 2D, and the gap opening in
  `element_interface_intpnt_materi_tension_status`. The GNU follows the
  manual (owner decision 2026-11-02; DIVERGENCES
  INTERFACE-INTPNT-DIRECTION-GAP, class A).
