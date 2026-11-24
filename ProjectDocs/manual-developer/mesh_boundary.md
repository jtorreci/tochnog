# mesh_boundary + element_boundary + node_boundary (developer)

## Where

- `area.cc`, `void mesh_boundary_set( void )` (renamed from
  `contact_node_boundary_set`, F5 T2 V3): one scan for both callers.
  - `mesh_boundary -yes` and no `ELEMENT_BOUNDARY` record yet: fills
    `NODE_BOUNDARY` and `ELEMENT_BOUNDARY` (`-yes`).
  - otherwise, `CONTACT_TARGET_ELEMENT_GROUP 0` active and no
    `NODE_BOUNDARY` record: fills `NODE_BOUNDARY` only (unchanged K05
    behavior, see `contact_target_element_group.md`).
- `top.cc`, `iteration_start()`: calls it after `element_normal_set()`,
  serially, before `element_loop()` and `parallel_contact()`.
- `mesh.cc`, `mesh_has_changed()`: deletes `ELEMENT_BOUNDARY`
  (VERSION_NORMAL); `NODE_BOUNDARY` goes with the NODE-class delete. The
  next `iteration_start()` regenerates both for the new mesh.
- `database.cc`: `MESH_BOUNDARY` (INTEGER, length 1, no index) and
  `ELEMENT_BOUNDARY` (INTEGER, length 1, class ELEMENT, required ELEMENT).
  `NODE_BOUNDARY` existed already.
- Enums appended before `LAST_DUMMY` in `tochnog.h` / `tochnog-mod.h`.

## Implementation details

- Side tables: the static `element_normal_sides()` of area.cc. A node ->
  side-bearing element compressed list gives the elements around the
  first node of a side; a side is free when no other element holds all
  its nodes.
- `element_boundary` is written for the owner of every free side
  (`length = 1` PUT inside the side loop; writing twice is harmless).
- Hand-written `node_boundary` records are kept: with `mesh_boundary -yes`
  the computed nodes are added to them; for the contact path they disable
  the generation.
- The `-yes` gate is "no ELEMENT_BOUNDARY record": a mesh with no
  side-bearing element rescans every iteration (cheap, no records).

## Hardcoded parameters / pending

- Element types without a side table (prism6/15, tet10, hex20, quad8,
  tria15, ...) are ignored. Extending `element_normal_sides()` covers both
  `element_normal` and this scan.
- Divergences: MESH-BOUNDARY-DEFAULT (class A, default -no = manual) and
  ELEMENT-BOUNDARY-ALL (class D, Pro marks every element).

## Verification

Professional 02-08-2026 probes (3 x 3 and 4 x 4 quad4, with and without
`mesh_boundary`, `mesh_boundary -no`): `node_boundary` = the outer ring in
both programs. Tests `mesh_boundary1`, `mesh_boundary_default1`
(test-2026, checks in `scripts/build_safe.sh`).
