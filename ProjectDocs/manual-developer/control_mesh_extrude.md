# control_mesh_extrude

## Files and functions

- `extrude.cc` — `extrude()`: called from `top.cc` `step_start()` on every
  step (task YES and NO); it runs `mesh_extrude(icontrol)` once (static
  flag), when `ICONTROL` is an index with a `CONTROL_MESH_EXTRUDE` record.
- `mesh.cc` — `mesh_extrude(long int icontrol)`: reads every record at that
  index and builds the 3D mesh.
- `input.cc` — the integer reader accepts a negative number in
  `CONTROL_MESH_EXTRUDE_ELEMENT_GROUP_NEW` (`-1` would otherwise be read
  as a word).
- `database.cc` — `CONTROL_MESH_EXTRUDE` (DOUBLE, `DATA_ITEM_SIZE`),
  `CONTROL_MESH_EXTRUDE_N`, `_DIRECTION` (INTEGER, 1, words `-y` / `-z`),
  `_ELEMENT_GROUP_NEW`, `_CONTACTSPRING_ELEMENT_GROUP` (requires `_NEW`),
  `_CONTACTSPRING_ELEMENT_GROUP_NEW` (INTEGER, `DATA_ITEM_SIZE`); all but
  the first require `CONTROL_MESH_EXTRUDE`. The names of the contact
  spring records are the Professional spelling
  (`control_mesh_extrude_contact_spring_element_group(_new)`).
- `tochnog.h` / `tochnog-mod.h` — the enums after `CONTROL_MESH_EXTRUDE_N`.

## Implementation details

`mesh_extrude(icontrol)`:

1. Reads `zv[0..nz-1]`; one value means `0 zv[0]` (GNU extension). `n_sub[i]`
   from `_n` (default 1, at least `nz-1` values, each >= 1). `axis` = 2, or 1
   for `-y`.
2. `max_node` / `max_elem` from `db_highest_index` (`db_max_index` is the
   allocated range: on a 9-node mesh it gave 10 and shifted the plane
   numbers). `nnod` / `nelem` = highest number, plus 1 when number 0 is
   active. `attached[]` marks the nodes of the `-bar2`, `-tria3`, `-quad4`,
   `-quad9`, `-tria6` elements.
3. Layers: `nl = sum n_sub`, `lay_int[l]` = interval of layer `l`. Planes:
   `ppl` = 2 with a `-quad9` or a `-tria6` -> `-prism18` in the mesh (a
   mid plane per layer), else 1;
   `np = nl*ppl + 1` planes at `zp[]`, evenly spaced inside each interval.
4. `_element_group_new`: `ng` values, `nold = ng/(1+n_int)` (error when not
   divisible); group of layer `l` of an element of old group `k` =
   `gtab[nold + k*n_int + lay_int[l]]`.
5. Contact spring groups (`sgroup[inod]`, index into `_new`), computed
   before the elements change: one new value = every attached node;
   otherwise the lowest index `k` of the listed old groups of the
   elements attached to the node.
6. Planes: plane 0 is the 2D node, its coordinate `axis` set to `zp[0]`;
   plane `p` is node `p*nnod + inod`, a copy of every node-class record
   (`mesh_copy_node_item`) with `coords[axis] = zp[p]`. For `-y` the 2D y
   moves to z first. `NODE_START_REFINED` gets the same coordinates.
   Planes are written from the top down, so plane 0 is changed last.
7. Elements: layer `l` of element `e` is `l*nelem + e`, created from `e`
   with `create_element()` (copies the element data), layers from the top
   down so layer 0 overwrites `e` last; a new group is PUT after the
   create; a negative group skips the layer (layer 0: `delete_element(e)`).
   Base / top planes `l*ppl`, `(l+1)*ppl`, mid plane `l*ppl+1`.
   Node orders: quad4 `a b a' b' c d c' d'` (Professional; the 2D eta
   direction becomes the hex8 zeta direction, so the two sides of an
   interface quad4 stay the two hex8 faces), prism6 `a b c a' b' c'`,
   quad4 from bar2 `a b a' b'`, hex27 base / mid / top in the quad9 order.
8. Springs: `-CONTACTSPRING` `inod` -> `(np-1)*nnod + inod`, numbered after
   `max(highest element, nl*nelem-1)`.
9. Nodes of the planes that no element uses are deleted (`delete_node`);
   `mesh_has_changed()`.

Measured on the Professional 02-08-2026 (same `.dbs` mesh: node numbers
and coordinates, element numbers, node orders, groups): `other/extrude1`
(45 elements, 119 nodes), `other/interface_bar2_hex8` (also the nodal
dofs to 7.6e-12, after `interface_convert`; only the node order of the
two converted interface hex8 differs, `interface_convert`), and the probes
of the tests below. `force11`, `extrude_quad9_shared_nodes1` and
`post_point_extrude_hex27_1`: the same coordinates and nodal results as
before U5a (by coordinates, 7e-9 of a scale of 819), new numbers.
`other/extrude2` still stops in both binaries (GNU: `point_el` has no
`-prism6` for its `post_point`, BACKLOG C125). `other/safety2` and
`validation_14_mesh` cannot run in the GNU for other reasons
(`safety_slip`, `change_dataitem_time_until_data`).

Tests (`validation-suite/test-2026`, `build_safe.sh` `ext_check`: the
`.dbs` element / group / node lines equal the Professional mesh):
`extrude_layers1`, `extrude_direction1`, `extrude_element_group_new1`,
`extrude_contact_spring1`, `extrude_contact_spring2`, `extrude_single1`;
input errors `extrude_error1` (group table length), `extrude_error2`
(direction), `extrude_error4` (`_n` too short).

## History

- F5 T4 U5a (2026-11-08): rewritten from `mesh_extrude(z_layer, n_layer,
  quad9_layers)`. The old linear branches made one layer per z value from
  the untouched z = 0 nodes (a collapsed first layer for `0. 1.`), ignored
  `_n`, kept the source `-bar2`, used the hex8 order `a b c d a' b' c'
  d'` and `nbase = max_node+1` blocks with gaps; the quad9 branch took
  only `_n[0]` layers over `z0..z_last`. `extrude()` applied the first
  record of any index at the first timestep.
- F4 T12 (2026-10-11, BACKLOG C47): shared hex27 inner-plane nodes.
- F4 T6f C57 (2026-10-13): `NODE_START_REFINED` gets the plane coordinate.

## Hardcoded parameters / pending refactorings

- The extruded element types are fixed (`-bar2`, `-tria3`, `-quad4`,
  `-quad9`, `-tria6`); other elements are left as they are and their nodes
  are not copied. `-tria6` gives `-prism18` (`control_mesh_extrude_element`,
  see its page) or `-prism12` without mid planes (F5 T4 U9c).
- `-y` is the only other direction (manual 6.180).
- The node-data copy loop could be shared with `mesh_mirror`,
  `mesh_copy` and `mesh_rotate_3d`.
