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

## Sprint 8 additions (2026-08-24): region filters + interface resets

- Enums CONTROL_RESET_ELEMENT_DOF/_ELEMENT_GROUP/_GEOMETRY/_INTERFACE/
  _INTERFACE_STRAIN/_NODE (headers in sync, 1059).
- data.cc: the node filter is computed ONCE per control_reset_dof record
  before the value variants (constant/dof-diagram/spatial): a marker
  array reset_dof_node_filter[1+max_node] built from the geometry test
  (elements fully inside) OR the node list (elements with all nodes
  listed), optionally restricted by element groups; every node loop of
  the value variants applies `filter[inod]` (2 loops with 12-space
  indent + 1 with 11 — the file mixes indentation styles; match
  exactly when editing).
- Interface resets: separate block in data() BEFORE the control_reset_dof
  loop; geometry test per element (all nodes inside), then PUT 0. to
  ELEMENT_INTERFACE_STRAIN_NORMAL (both variants) and additionally
  ELEMENT_INTERFACE_FORCE_TANG/_TANG2 (full reset). Normal PUTs here are
  safe: data() runs in step_close, and these records are only written
  (never read back with GET_IF_EXISTS gating behaviour).
- control_reset_element_dof: registered but the -yes branch (only
  element_dof) is NOT wired — the GNU reset path writes NODE_DOF;
  routing to ELEMENT_DOF needs the per-integration-point storage.
  Documented as partial.

## Verification

creset_geom (hisv0 0.5 -> 0 left column only; A/B without the filter
fails), creset_iface (strain history 0.0 exact). Suite 87/87.
