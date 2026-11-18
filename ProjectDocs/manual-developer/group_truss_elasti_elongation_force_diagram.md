# group_truss_elasti_elongation_force_diagram — developer

## Files and functions

- `database.cc` — `GROUP_TRUSS_ELASTI_ELONGATION_FORCE_DIAGRAM`
  (DOUBLE_PRECISION, DATA_ITEM_SIZE, class TRUSS, requires GROUP_TYPE);
  enum after `GROUP_TRUSS_YOUNG` in `tochnog.h` / `tochnog-mod.h`.
- `truss.cc` — `truss()`: after the elastic update
  `new_truss_force = old + k * incremental_length`, the diagram replaces
  the force with `table_xy(l)`, l = (new_length - initial_length) /
  initial_length, and the stiffness with the segment slope /
  initial_length (0 outside the diagram). The plasticity and rope limits
  follow unchanged.

## Implementation details

`new_length` is the current length of each memory branch: the real
length for `-updated` and `-total`, the length projected on the initial
axis for the linear memories. The Professional gives the same force with
`materi_displacement` and with `materi_velocity_integrated`, after one or
two steps, and after unloading (truss_force_diagram2/3): the diagram works
on the total strain, not on increments.

## History

- F5 T5 U6 (2026-11-19): added.

## Hardcoded parameters / pending refactorings

- Thermal expansion (`group_truss_expansion`) is not subtracted from l
  (not measured).
- Outside the diagram the stiffness is 0; a free truss end there has no
  stiffness from the truss.
