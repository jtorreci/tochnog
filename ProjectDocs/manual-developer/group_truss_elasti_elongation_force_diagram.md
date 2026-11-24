# group_truss_elasti_elongation_force_diagram — developer

## Files and functions

- `database.cc` — `GROUP_TRUSS_ELASTI_ELONGATION_FORCE_DIAGRAM`
  (DOUBLE_PRECISION, DATA_ITEM_SIZE, class TRUSS, requires GROUP_TYPE);
  enum after `GROUP_TRUSS_YOUNG` in `tochnog.h` / `tochnog-mod.h`.
- `truss.cc` — `truss()`: after the elastic update
  `new_truss_force = old + k * incremental_length`, the diagram replaces
  the force with `table_xy(l)`, l = (new_length - initial_length) /
  initial_length, and the stiffness with the segment slope /
  initial_length (0 outside the diagram). The accumulated thermal strain
  (`ELEMENT_TRUSS_STRAIN_TEMPERATURE` of the last step plus
  `group_truss_expansion * inc_temperature`) is subtracted from l first.
  The plasticity and rope limits follow unchanged.

## Implementation details

`new_length` is the current length of each memory branch: the real
length for `-updated` and `-total`, the length projected on the initial
axis for the linear memories. The Professional gives the same force with
`materi_displacement` and with `materi_velocity_integrated`, after one or
two steps, and after unloading (truss_force_diagram2/3): the diagram works
on the total strain, not on increments.

## History

- F5 T5 U6 (2026-11-19): added.
- 2026-11-23: the thermal strain is subtracted (TRUSS-DIAGRAM-THERMAL).

## Hardcoded parameters / pending refactorings

- Outside the diagram the stiffness is 0; a free truss end there has no
  stiffness from the truss.

## Differences with the Professional version

- **TRUSS-DIAGRAM-THERMAL (deliberate, owner 2026-11-23).** The GNU enters
  the diagram with the mechanical strain, the total strain minus alpha dT.
  The Professional 02-08-2026 uses the total strain, so a free truss that is
  only heated gets a force. Probe: length 2, total strain 0.014, alpha 1e-3.
  The Professional gives force 12 for dT +10, -10 and 0. The GNU gives 4 for
  dT +10 (`truss_force_diagram_thermal1`) and 14.5 for dT -5 in four steps
  (`truss_force_diagram_thermal2`). Without the diagram, both subtract the
  thermal strain (3.4 = 3 + 100 * 0.004). See `DIVERGENCES-PROFESSIONAL.md`.

