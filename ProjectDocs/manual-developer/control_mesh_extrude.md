# control_mesh_extrude

## Files and functions

- `mesh.cc` — `mesh_extrude(double z_layer[], long int n_layer, long int quad9_layers)`.
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

## quad9 -> hex27 branch (shared inner-plane nodes)

`mesh_extrude()`, branch `el[0]==-QUAD9 && nnol==9` (third argument
`quad9_layers` = first value of `control_mesh_extrude_n`, read in
`extrude.cc`; `nq = max(quad9_layers, 1)` hex27 layers over
`z_layer[0]..z_layer[n_layer-1]`). Plane numbering, with
`next_node = nbase + n_layer*nbase` (first index after the boundary
copies) and `nbase = max_node + 1`:

- plane 0 (z0): the source nodes `nodes[inol]`;
- plane q = 1..2nq-1 (z ascending; odd q = mid-planes at the layer
  mid-height, even q = layer tops): node `next_node + (q-1)*nbase +
  nodes[inol]`;
- plane 2nq (z1): the existing boundary copy `nodes[inol] + n_layer*nbase`.

Layer `lay` uses `mid0 = next_node + 2*lay*nbase` and, below the last
layer, `top_abs = next_node + (2*lay+1)*nbase`; its base is plane 0 for
`lay==0` and the previous layer's `top_abs` (`last_top`) otherwise. A
plane node is created (`db(NODE, ..., PUT)` with the z coordinate, plus
every other node-class record copied from the source node with
`mesh_copy_node_item()`) only when `db_active_index(NODE, nn)` is false,
so the first element that reaches it creates it and its neighbours reuse
it. The hex27 connectivity is base `base+nodes[inol]`, mid `mid0+nodes[inol]`,
top `top_abs+nodes[inol]` (or the boundary copy) for inol = 0..8 in the
lexicographic quad9 order.

History: before F4 T12 (2026-10-11, BACKLOG C47) `top_abs` and `mid0`
were fresh per-element blocks of 9 nodes (`next_node++`), so the hex27
elements were disconnected at every inner plane (Professional `force11`:
15600 nodes instead of 12000, planes 2400/3600/3600/3600/2400 instead of
2400 each, ring 7.9x too soft, `mom1y_sig` 5.80223 vs 5.76460 because
the section moment is equilibrium-based). Test:
`validation-suite/test-2026/extrude_quad9_shared_nodes1.dat`
(`build_safe.sh`: rc 0 + 75 node records in the .dbs).

## External dependencies

- `create_element()` in `create.cc` (copies element data: `element_dof`,
  `element_dof_initialised`, `nonlocal`, ...).
- `delete_element()` in `delete.cc`.
- `mesh_has_changed()` in `mesh.cc` / `top.cc`.

## Hardcoded parameters / pending refactorings

- `control_mesh_extrude_n`: only its first value is used, and only by
  the quad9 -> hex27 branch (number of hex27 layers); the linear branches
  create one element per z segment.
- The source 2D mesh is assumed to lie at z=0; the extrusion is always along
  the z-axis.
- Only `-tria3`, `-bar2`, `-quad4` and `-quad9` are handled; other 2D
  element types are ignored (their nodes are still copied).
- quad9 -> hex27 node numbering has gaps (`nbase` blocks) and is not the
  Professional's compact plane-major numbering; renumbering would also
  move the boundary-copy block shared with the linear branches.
- hex27 nodal stress recovery differs from the Professional at stress
  peaks (BACKLOG C57).
- The generic node-data copy loop could be factored out and reused by
  `mesh_mirror`/`mesh_copy`/`mesh_rotate_3d` (they currently each have their
  own copy logic).
- Requires `number_of_space_dimensions 3` and sufficient
  `number_of_integration_points` (8 for hex8, 6 for prism6).
