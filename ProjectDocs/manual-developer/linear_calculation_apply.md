# linear_calculation_apply — developer

- Record `LINEAR_CALCULATION_APPLY` (INTEGER, length 1, class CONTROL, no
  index; `database.cc` before `control_materi_dynamic`, enum after
  `INTERFACE_GAP_APPLY`).
- `linear_calculation_apply_transform()` (`general.cc`), called once in
  `top()` right after `groundflow_apply_transform()`, before the start
  version is stored. Order: the material substitutions (they read the
  records the plasticity deletion removes), plasticity, `OPTIONS_MESH`,
  the `group_*_memory` records (found by name: prefix `group_`, suffix
  `_memory`, INTEGER), the dependencies (`DEPENDENCY_ITEM[0]` a `group_*`
  item, `[2]` a dof label), then the deletions.
- Helpers: `general_delete_prefix()` (deletes every record whose name
  starts with a prefix), `table_xy()` (`math.cc`) on a copy of the spring
  table (it sorts in place).
- Hardcoded: Poisson 0.2 for wolffersdorff (manual); the value positions
  E50ref = `hardsoil[0]`, E0 = `young_power[0]`, hs = `wolffersdorff[1]`.
- Pending: the four manual records the GNU does not have (see the user
  page); records added later by `control_input` or the runtime file are
  not rewritten.

## Differences with the Professional version

`LINCALC-WOLFF`, `LINCALC-CONTACT` (class D), see the user page.
