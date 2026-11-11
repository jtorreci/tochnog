# control_mesh_map

## Files and functions

- `mesh_map.cc` — `mesh_map(long int first_of_control)`: called from
  `top.cc` `step_start()` after `extrude()` (so after
  `mesh_convert_hex20` and the extrusion) on every step; it acts on the
  first step of a control index (`first_of_control`) whose
  `CONTROL_MESH_MAP` record is `-yes`, also an index without
  `control_timestep`, before its timestep.
  - `mesh_map_linear()`: quadratic name / node count -> linear name, node
    count and corner positions (`-hex27` 1 3 7 9 19 21 25 27, `-tet10`
    1 3 6 10, `-prism15` 1-6).
  - `mesh_map_weights()`: the linear shape functions at the natural
    coordinates of a quadratic node (hex27 lattice product of 1, 1/2, 0;
    edge mid node 1/2 + 1/2).
  - `mesh_map_element_put()`: gives an element its new node list in place.
- `database.cc` — `CONTROL_MESH_MAP` (INTEGER 1, CONTROL),
  `ELEMENT_MAP` (INTEGER `1+MNOL`, ELEMENT: the quadratic element, the
  Professional `.dbs` record) and `ELEMENT_MAP_NODE` (DOUBLE
  `2*MNOL*MDIM`, ELEMENT, GNU: `NODE` then `NODE_START_REFINED` of the
  quadratic nodes).
- `tochnog.h` / `tochnog-mod.h` — the three enums after `TRANSLATE`, the
  prototype after `extrude()`.
- `makefile`, `makefile - Windows-Mingw64` — `mesh_map` object.

## Implementation details

`mesh_map()`:

1. Direction: to linear when any element is a `-hex27` / `-tet10` /
   `-prism15`, else back (the elements with an `ELEMENT_MAP`).
2. `used[]`: the nodes the current elements use (sized up to the highest
   node of any `ELEMENT_MAP`).
3. To linear, per quadratic element: `ELEMENT_MAP` = the element,
   `ELEMENT_MAP_NODE` = the coordinates of its nodes, element = the linear
   element of the corners. `mesh_has_changed()` then runs `nod_nod`, which
   deletes the mid nodes (no element uses them).
4. Back, per linear element whose `ELEMENT_MAP` corners are its nodes:
   the corner `NODE_DOF` / `NODE_DOF_START_REFINED`, then for every node
   of the quadratic element not in `used[]` the interpolated dofs; a
   missing node is PUT again (`NODE`, `NODE_START_REFINED` from
   `ELEMENT_MAP_NODE`, `NODE_DOF`, `NODE_DOF_START_REFINED`). The
   element gets the quadratic node list. A shared mid node is set once.
5. `mesh_has_changed()` when an element changed.

`mesh_map_element_put()` deletes the records tied to the old node list or
integration points (`ELEMENT_DOF`, `element_intpnt_*`, `element_matrix_*`,
`ELEMENT_RHSIDE_DELETE`, the nonlocal records), deletes `ELEMENT` before
the PUT (`db_len` finds a record length by bisection on the `LONG_MIN` end
mark, so a shorter PUT over a longer record gives a wrong length), PUTs
the new list, a zero `ELEMENT_INTPNT_DOF` and `ELEMENT_DOF_INITIALISED`
0, so the element initialises its integration point dofs from the nodes
(what `create_element` does for a remeshed element). The other records
(group, geometry, activity, ...) stay. `create_element` / `create_node`
are not used in place: `create_element` copies every element record
through a `DATA_ITEM_SIZE` buffer, too short for the `element_matrix_*`
records of a solved element (heap corruption seen), and `create_node`
copies through `db_dbl` pointers that a PUT beyond a record's highest
index reallocates.

Interpolation at the natural coordinates (not `point_el` on the
coordinates): exact, no iteration, and the same as the Professional on
the curved ring of `other/map1` (mid node displacements to 1.6e-11 of a
scale of 1.9e-5 after the map back).

Measured on the Professional 02-08-2026 (`other/map1` and probes, `.dbs`
of every stage with `control_print_database -separate_index`): element
lists after each map, `element_map` kept after the map back, a further
map toggles again, the mid node dofs after the map back = trilinear
interpolation of the corner dofs (all 12 dofs, 4.5e-12), the mid node
dofs unchanged by the linear calculation, a `-hex20` -> `-hex8 1-8`,
`-prism18` / `-prism12` / 2D `-quad9` / `-tria6` unchanged, a map at an
index with a timestep runs before it, `global_element_dof_apply` not
checked. GNU = Professional on the hex27 column: corner dofs after the
linear calculation 1.3e-12, all dofs after the map back 2.5e-12, final
3.7e-9 (bicg), relative to the largest dof.

Tests (`validation-suite/test-2026`, `build_safe.sh`): `map_mesh_check`
(element / element_map lines and node count), `map_dof_check` (node dofs
against the Professional values, `tol*(1+|value|)`), `map_mid_check`
(mid node dofs = mean of the edge corners): `control_mesh_map1` -
`_map3`, `_mixed1`, `_toggle1`, `_no1`, `_keep1`, `_keep2`, `_tet10_1`,
`_prism15_1`, `_hex20_1`, `_tet10_2`, `_prism15_2`.

## History

- F5 T4 U9d (2026-11-11): new (owner decision, implement manual 6.210).

## Hardcoded parameters / pending refactorings

- Mapped types fixed: `-hex27`, `-tet10`, `-prism15` (and `-hex20`
  through `mesh_convert_hex20`). `-prism18` / `-prism12` / 2D elements
  are left as they are (= Professional; a mesh mixing them with mapped
  elements is not conforming while linear).
- The node coordinates of a mapped element are the ones of the map to
  linear; a mesh moved in between (updated formulations) is not followed
  by the recreated mid nodes.
- With `-separate_index` the GNU `control_print_database` still writes one
  `.dbs` (pre-existing), so the tests stop the run at the stage they check.

## Differences with the Professional version

See the user page: MESH-MAP-MID-NODES, MESH-MAP-TET10, MESH-MAP-PRISM15,
MESH-MAP-HEX20; BACKLOG C134-C136.
