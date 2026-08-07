# control_mesh_extrude

## Files and functions

- `mesh.cc` — `mesh_extrude(double z_layer[], long int n_layer)`.
- `extrude.cc` — dispatch in `extrude()`, which reads `CONTROL_MESH_EXTRUDE`
  and calls `mesh_extrude()`. Called from `top.cc`.
- `database.cc` — keywords `CONTROL_MESH_EXTRUDE` (DOUBLE_PRECISION,
  `data_length=DATA_ITEM_SIZE`, `fixed_length=0`) and
  `CONTROL_MESH_EXTRUDE_N` (INTEGER, `data_length=DATA_ITEM_SIZE`,
  `fixed_length=0`, `data_required=CONTROL_MESH_EXTRUDE`).
- `tochnog.h` / `tochnog-mod.h` — enum entries `CONTROL_MESH_EXTRUDE`,
  `CONTROL_MESH_EXTRUDE_N`.

## Implementation details

`mesh_extrude(z_layer, n_layer)` sweeps a flat 2D mesh (z=0) to 3D along the
z-axis. Steps:

1. Compute the copy base `nbase = max_node + 1`. The copy of node `inod` in
   layer `layer` is at `new_node = nbase + layer*nbase + inod`
   (`inod + (max_node+1)*(layer+1)`), so all copies start after the original
   node block.
2. Create the extruded node copies for every layer boundary. For each copy the
   z-coordinate is set from `z_layer[layer]` and all other node-class data
   items are copied with a generic loop:
   ```c
   for ( int idat=0; idat<MDAT; idat++ )
     if ( idat!=NODE && db_data_class(idat)==NODE &&
          db_active_index( idat, inod, VERSION_NORMAL ) ) {
       long int ndata_len = db_len( idat, inod, VERSION_NORMAL );
       if ( db_type(idat)==DOUBLE_PRECISION ) { double *dold=db_dbl(...); db(...,PUT); }
       else { long int *iold=db_int(...); db(...,PUT); }
     }
   ```
   `NODE_START_REFINED` additionally gets the extruded coordinate.
3. Create one 3D element per 2D element per layer with `create_element()`:
   `-tria3` → `-prism6` (bottom layer nodes `layer*nbase`, top layer nodes
   `(layer+1)*nbase`), `-quad4` → `-hex8`, using the same node ordering.
4. Delete the 2D source elements with `delete_element()` (not valid in 3D) and
   signal `mesh_has_changed()`.

Dispatch in `extrude.cc`:
```c
db( ICONTROL, 0, &icontrol, ddum, ldum, VERSION_NORMAL, GET );
if ( db_active_index( CONTROL_MESH_EXTRUDE, icontrol, VERSION_NORMAL ) ) {
  db( CONTROL_MESH_EXTRUDE, icontrol, idum_e, ext_z, ext_length,
    VERSION_NORMAL, GET );
  n_layer = ext_length;
  if ( n_layer>0 ) mesh_extrude( ext_z, n_layer );
}
```

The original `extrude.cc` (which aborted with "extrude is not available") was
rewritten to read `CONTROL_MESH_EXTRUDE` and delegate to `mesh_extrude`.

## External dependencies

- `create_element()` in `create.cc` (copies element data: `element_dof`,
  `element_dof_initialised`, `nonlocal`, ...).
- `delete_element()` in `delete.cc`.
- `mesh_has_changed()` in `mesh.cc` / `top.cc`.

## Hardcoded parameters / pending refactorings

- `control_mesh_extrude_n` (number of elements per layer) is defined in the
  database and enum but not used; only 1 element per layer is created.
- The source 2D mesh is assumed to lie at z=0; the extrusion is always along
  the z-axis.
- Only `-tria3` and `-quad4` are handled; other 2D element types are ignored
  (their nodes are still copied).
- The generic node-data copy loop could be factored out and reused by
  `mesh_mirror`/`mesh_copy`/`mesh_rotate_3d` (they currently each have their
  own copy logic).
- Requires `number_of_space_dimensions 3` and sufficient
  `number_of_integration_points` (8 for hex8, 6 for prism6).
