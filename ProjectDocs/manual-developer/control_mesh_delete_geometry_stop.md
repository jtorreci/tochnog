# control_mesh_delete_geometry_stop

## Files and functions

- `tochnog.h` / `tochnog-mod.h` — two new enum entries appended before
  `LAST_DUMMY` (values stay stable):
  `CONTROL_MESH_DELETE_GEOMETRY_STOP`,
  `CONTROL_MESH_DELETE_GEOMETRY_STOP_GEOMETRY`.
- `database.cc` — name/type/class registration of
  `control_mesh_delete_geometry_stop` (INTEGER, 1) and
  `control_mesh_delete_geometry_stop_geometry` (INTEGER, 2,
  `data_required` = `_stop`).
- `delete.cc` — `void delete_geom( double time_current )`: the stop block
  runs before the `CONTROL_MESH_DELETE_GEOMETRY` gate, so the stop index
  does not need to carry the geometry record.
- `elem.cc` — `elem()`, the `ELEMENT_INTPNT_DOF` write.

## Implementation details

### The stop (delete.cc)

The block reads `CONTROL_MESH_DELETE_GEOMETRY_STOP` at the current
`icontrol` (the record may live at a different index than the geometry and
factor records). When the switch is `-yes`:

1. Optionally read `CONTROL_MESH_DELETE_GEOMETRY_STOP_GEOMETRY`
   (`item, index`). The membership test is the same as
   `control_mesh_delete_geometry`: an element is part of the geometry when
   ALL its nodes lie inside it. The node membership uses `geometry()` with
   `CONTROL_MESH_DELETE_GEOMETRY` as the projection type (filled interior
   of area geometries).
2. For every element with an active `ELEMENT_DELETE_FACTOR` (restricted by
   the stop geometry when given), `db_delete_index` removes
   `ELEMENT_DELETE_FACTOR` and `ELEMENT_RHSIDE_DELETE`.

Elements removed completely (`delete_element`, factor end value above
`1 - 1.e-10`, manual 6.170) have no element record and are not reactivated.

`delete_geom()` is called from `step_start()` for every control index
(with or without `control_timestep`), so a stop index that carries no
timestep still runs (`top.cc` calls `step_start(NO, ...)` in that case).

### Freezing the deleted element's integration-point state (elem.cc)

Manual 6.174 says the reactivated element takes its strains/stresses from
the moment just before being deleted, and `global_element_dof_apply`
defaults to `-yes` (manual 6.553): the element remembers its own
integration-point values instead of the averaged nodal values.

The GNU `elem()` recomputed `ELEMENT_INTPNT_DOF` on every iteration from the
current nodal dofs, so a partially deleted element's stored stress drifted
away from the deletion moment (measured on `delete2.dat`: the element
recomputed to `1.4` while the Professional's `element_intpnt_dof_average`
stays `1.0`). The fix skips the `ELEMENT_INTPNT_DOF` write once the released load
has been captured (`ELEMENT_RHSIDE_DELETE` active). On the first delete
iteration the write still runs (storing the pre-deletion value); from then
on `ELEMENT_INTPNT_DOF VERSION_NEW` keeps the value copied from
`VERSION_NORMAL` at the step start, and the step-end version copy promotes
it. The deleted element therefore keeps the pre-deletion state until the
stop reactivates it.

### Measured status (delete2.dat, Professional 02-08-2026)

`delete2.dat` (1D, two `bar2` in series, `control_reset_dof -sigxx 1`,
delete the right bar 40 %, stop, continue) now parses and runs. The
element integration-point stress matches the Professional
(`element_dof = 0.8` in both elements after the stop), and the `node_dof`
target (`node 2 sigxx = 0.8`) is CLOSED: the nodal stress recovery was
aligned with the Professional's absolute weighted average in the
2026-09-15 work unit (see `nodal_stress_recovery.md`). Before that change
the GNU recovery accumulated the element stress INCREMENTS, so at a node
shared by two elements whose increments cancel the reported stress kept
its pre-step value (`0.6`); the Professional assigns the absolute average
(`node_rhside 1.6 / node_lhside 2.0 = 0.8`). `delete2.dat` is now `rc=0`,
`delete3.dat` is unchanged (`sigyy = 0.6`) and the full `node_dof` field of
the delete family matches the Professional exactly.
