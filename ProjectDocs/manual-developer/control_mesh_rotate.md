# control_mesh_rotate

## Files and functions

- `mesh.cc` — `mesh_rotate_3d(long int nrot)`.
- `top.cc` — dispatch in `step_start` when `CONTROL_MESH_ROTATE` is active.
- `database.cc` — keyword `CONTROL_MESH_ROTATE` (INTEGER, `data_length=1`).
- `delete.cc` — `delete_element()` used to remove the 2D source elements.

## Implementation details

`mesh_rotate_3d(nrot)` sweeps a flat 2D mesh to 3D around the original y-axis
(which becomes the new z). Steps:

1. Duplicate every node with the rotation `(x,y,z) -> (z,y,-x)` (a 90° sweep).
   For each duplicated node, all node-class data items are copied with a
   generic loop:
   ```c
   for ( idat=0; idat<MDAT; idat++ )
     if ( idat!=NODE && db_data_class(idat)==NODE &&
          db_active_index(idat, inod, VERSION_NORMAL) ) {
       long int len = db_len(idat, inod, VERSION_NORMAL);
       if ( db_type(idat)==DOUBLE_PRECISION ) { double *d=db_dbl(...); db(...,PUT); }
       else { long int *i=db_int(...); db(...,PUT); }
     }
   ```
   This copies `node_dof`, `node_dof_start_refined`, `node_start_refined`, etc.,
   avoiding the per-item fragility seen in the earlier `mesh_mirror` attempt.
   `NODE_START_REFINED` additionally gets the rotated coordinates.
2. Duplicate each 2D element: `-tria3` → `-prism6`, `-quad4` → `-hex8`, using
   `create_element()` (which copies all element data: `element_dof`,
   `element_dof_initialised`, `nonlocal`, ...). The new element connects the
   original nodes (bottom) with the rotated copies (top).
3. Delete the 2D source elements with `delete_element()` (they are not valid
   in 3D), matching the manual's "all data not valid in 3D will be deleted".

Dispatch in `top.cc`:
```c
if ( db_active_index( CONTROL_MESH_ROTATE, icontrol, VERSION_NORMAL ) ) {
  long int nrot=1;
  db( CONTROL_MESH_ROTATE, icontrol, &nrot, ddum, ldum, VERSION_NORMAL, GET );
  mesh_rotate_3d( nrot );
}
```

## External dependencies

- `create_element()` in `create.cc` (copies element data).
- `delete_element()` in `delete.cc`.

## Hardcoded parameters / pending refactorings

- Only `nrot = 1` (one rotational segment) is implemented. Multi-segment sweep
  (`n > 1`, repeating the rotation with `angle = 360/n`) is pending.
- The sweep rotation is fixed to 90° (`(x,y,z) -> (z,y,-x)`). The
  `CONTROL_MESH_ROTATE_ANGLE` is handled separately by `mesh_rotate_2d`
  (flat in-plane rotation), not combined with the 3D sweep.
- Requires `number_of_space_dimensions 3` and sufficient
  `number_of_integration_points` (8 for hex8, 6 for prism6).
- The generic node-data copy loop could be factored out and reused by
  `mesh_mirror`/`mesh_copy`/`mesh_rotate_3d` (they currently each have their
  own copy logic).

## Professional comparison (F4 T8, BACKLOG C46c, parked)

Valid probe: the Professional example `test/other/rotate1.dat` (2D
axisymmetric quad4 `1 3 2 4` with its side on the y axis,
`control_mesh_rotate 10 10`, `control_mesh_rotate_angle 10 180.`).
Professional: 14 nodes (axis nodes kept once, ring nodes numbered
`old + 4 k`), 5 `-prism6`, y -> z. GNU: rc 1, "Length too small of : node
... 2 ... at least 3" (`NODE` data length is `ndim`; `mesh_rotate_3d` writes
3 coordinates in a 2D database). Missing: run-time switch of
`number_of_space_dimensions` and of every ndim-sized record and dof list,
the n-segment sweep with the angle, axis-node collapse to prism6, the
data transfer/deletion of manual 6.238 and the `control_input` re-read.
Large; parked with this diagnosis (BACKLOG C46c, not scheduled: owner
decision 2026-10-10). `control_mesh_rotate_angle` is a
parameter of this sweep in the Professional, not a 2D rotation.
