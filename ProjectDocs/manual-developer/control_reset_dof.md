# control_reset_dof

## Implementación

- **Logic**: applied in `data()` in `data.cc`, after the `change_dataitem`
  block. For each active `CONTROL_RESET_DOF` record, iterate the nodes and
  modify the target dof in `NODE_DOF`.
- **Keywords** (data_class CONTROL) registered in `database.cc`:
  - `control_reset_dof` (INTEGER, DATA_ITEM_SIZE) — the dof to reset.
  - `control_reset_value_constant` (DOUBLE, length 1, required
    CONTROL_RESET_DOF).
  - `control_reset_value_dof` (INTEGER, length 1, required
    CONTROL_RESET_DOF) — the driving dof.
  - `control_reset_value_dof_diagram` (DOUBLE, DATA_ITEM_SIZE, required
    CONTROL_RESET_VALUE_DOF) — table z_i value_i.
  - `control_reset_value_method` (INTEGER, length 1, required
    CONTROL_RESET_DOF) — `-use`/`-add`/`-multiply`.
- **New enums**: `CONTROL_RESET_DOF`, `CONTROL_RESET_VALUE_CONSTANT`,
  `CONTROL_RESET_VALUE_DOF`, `CONTROL_RESET_VALUE_DOF_DIAGRAM`,
  `CONTROL_RESET_VALUE_METHOD` in `tochnog.h` / `tochnog-mod.h`, and a
  generic `MULTIPLY` enum (keyword `multiply`) for the `-multiply` method.
- **Dof resolution**: negative dof names (e.g. `-hisv0`) are resolved with
  `array_member(dof_label, ...)`; if `db_len(NODE_DOF, inod) == npuknwn`
  the index is divided by `nder` (matrix dofs).

## Diseño / decisiones

- Follows the `change_dataitem` / `control_distribute` pattern: iterate
  nodes, read `NODE_DOF`, modify, and write back implicitly (the pointer
  from `db_dbl` is in-place; no explicit PUT needed).
- The diagram is evaluated with `table_xy` (linear interpolation), same as
  `change_dataitem_time`.

## Detalles

- `data_length[NODE_DOF] = nuknwn`; the dof index from `array_member` is
  used directly on `node_dof` (no `*nder` unless the length is npuknwn).
- `db_error(CONTROL_RESET_DOF, ireset)` on out-of-range index.

## Pendiente

- Only `_value_constant`, `_value_dof` (+`_diagram`) and `_value_method`
  are implemented. The other `control_reset_value_*` distributions
  (exponent, linear, logarithmic, power, square root, multi_linear) are
  not.
- The reset applies on every call to `data()` (every time step); there is
  no time gate (Professional may apply it once).
